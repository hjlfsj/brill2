#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include <TChain.h>
#include <TMinuit.h>
#include <TString.h>

#include "external/cxxopts.hpp"
#include "include/config.h"
#include "include/energy_calculator/lost_energy_calculator.h"
#include "include/event/t0/t0_event.h"
#include "include/utils.h"

namespace {

constexpr int kLayerCount = 5;

struct FitSample {
	int calculator = -1;
	int de_layer = -1;
	int e_layer = -1;
	double de_raw = 0.0;
	double e_raw = 0.0;
};

struct CalculatorKey {
	int charge = 0;
	int mass = 0;
	int detector = -1;

	bool operator<(const CalculatorKey &other) const {
		if (charge != other.charge) return charge < other.charge;
		if (mass != other.mass) return mass < other.mass;
		return detector < other.detector;
	}
};

std::vector<FitSample> *g_samples = nullptr;
std::vector<std::unique_ptr<brill::LostEnergyCalculator>> *g_calculators = nullptr;

void PrintUsage(const cxxopts::Options &options) {
	std::cout << options.help() << "\n";
}

bool MatchParticle(
	const brill::T0Event &event,
	int index,
	int layer,
	int charge,
	int mass
) {
	return
		event.layer[index] == layer
		&& event.charge[index] == charge
		&& event.mass[index] == mass;
}

bool BuildSample(
	const brill::T0Event &event,
	int index,
	FitSample &sample
) {
	if (MatchParticle(event, index, 2, 2, 4)) {
		sample.de_layer = 0;
		sample.e_layer = 1;
	} else if (
		MatchParticle(event, index, 3, 2, 4)
		|| MatchParticle(event, index, 3, 4, 7)
		|| MatchParticle(event, index, 3, 6, 12)
	) {
		sample.de_layer = 1;
		sample.e_layer = 2;
	} else if (
		MatchParticle(event, index, 4, 2, 4)
		|| MatchParticle(event, index, 4, 4, 7)
		|| MatchParticle(event, index, 4, 6, 12)
	) {
		sample.de_layer = 2;
		sample.e_layer = 3;
	} else if (
		MatchParticle(event, index, 5, 2, 4)
		|| MatchParticle(event, index, 5, 4, 7)
	) {
		sample.de_layer = 3;
		sample.e_layer = 4;
	} else {
		return false;
	}

	sample.de_raw = event.energy[index][sample.de_layer];
	sample.e_raw = event.energy[index][sample.e_layer];
	return sample.de_raw > 0.0 && sample.e_raw > 0.0;
}

std::string DetectorNameFromLayer(int layer) {
	switch (layer) {
		case 0: return "t0d1";
		case 1: return "t0d2";
		case 2: return "t0d3";
		case 3: return "t0d4";
		case 4: return "t0s";
		default: return "";
	}
}

std::string CalculatorCachePath(
	const brill::AppConfig &config,
	int charge,
	int mass,
	const brill::SquareDetectorConfig &detector
) {
	return TString::Format(
		"%s/si_%.0fum_z%d_a%d.root",
		brill::JoinPath(config.workspace, config.paths.energy_calculator).c_str(),
		detector.thickness_um,
		charge,
		mass
	).Data();
}

int GetCalculatorIndex(
	const brill::AppConfig &config,
	int charge,
	int mass,
	int detector_layer,
	const brill::SquareDetectorConfig &detector,
	std::map<CalculatorKey, int> &indices,
	std::vector<std::unique_ptr<brill::LostEnergyCalculator>> &calculators
) {
	CalculatorKey key{charge, mass, detector_layer};
	auto iter = indices.find(key);
	if (iter != indices.end()) return iter->second;

	int index = int(calculators.size());
	calculators.push_back(std::make_unique<brill::LostEnergyCalculator>(
		charge,
		mass,
		brill::SiliconMaterial(),
		detector.thickness_um,
		CalculatorCachePath(config, charge, mass, detector)
	));
	indices[key] = index;
	return index;
}

void CalibrationFcn(
	int &,
	double *,
	double &value,
	double *parameters,
	int
) {
	value = 0.0;
	if (!g_samples || !g_calculators) return;

	for (const auto &sample : *g_samples) {
		double de_cal = parameters[sample.de_layer * 2] + parameters[sample.de_layer * 2 + 1] * sample.de_raw;
		double e_cal = parameters[sample.e_layer * 2] + parameters[sample.e_layer * 2 + 1] * sample.e_raw;
		if (de_cal <= 0.0 || e_cal <= 0.0) {
			value += 1e12;
			continue;
		}
		double de_phys = (*g_calculators)[sample.calculator]->EnergyLossFromResidual(e_cal);
		double residual = de_cal - de_phys;
		value += residual * residual;
	}
}

int WriteCalibrationParameters(
	const std::string &path,
	const double *parameters
) {
	std::filesystem::path file_path(path);
	if (!file_path.parent_path().empty()) {
		std::filesystem::create_directories(file_path.parent_path());
	}
	std::ofstream fout(path);
	if (!fout.good()) {
		std::cerr << "Error: Open output parameter file " << path << " failed.\n";
		return -1;
	}

	fout << "# layer p0 p1\n";
	for (int layer = 0; layer < kLayerCount; ++layer) {
		fout << layer << " " << parameters[layer * 2] << " " << parameters[layer * 2 + 1] << "\n";
	}
	return 0;
}

} // namespace

int main(int argc, char **argv) {
	cxxopts::Options options("calibrate_t0", "Calibrate T0 energies from tracked events.");
	options.add_options()
		("h,help", "Print help information.")
		("r,run", "Start run number.", cxxopts::value<int>(), "run")
		("e,end-run", "End run number.", cxxopts::value<int>(), "run")
		("t,trigger", "Trigger type.", cxxopts::value<std::string>(), "trigger")
		(
			"c,config",
			"Config file path.",
			cxxopts::value<std::string>()->default_value("config.toml"),
			"file"
		);

	auto result = options.parse(argc, argv);
	if (result.count("help")) {
		PrintUsage(options);
		return 0;
	}
	if (!result.count("run")) {
		std::cerr << "Error: Missing required option --run.\n";
		PrintUsage(options);
		return 1;
	}

	brill::AppConfig config;
	if (brill::LoadConfig(result["config"].as<std::string>(), config)) {
		return 1;
	}
	if (result.count("trigger")) {
		config.trigger = result["trigger"].as<std::string>();
	}

	const int run = result["run"].as<int>();
	const int end_run = result.count("end-run") ? result["end-run"].as<int>() : run;
	if (end_run < run) {
		std::cerr << "Error: end run " << end_run << " is smaller than run " << run << ".\n";
		return 1;
	}

	const brill::SquareDetectorConfig *detectors[kLayerCount] = {
		brill::FindDetectorConfig(config, "t0d1"),
		brill::FindDetectorConfig(config, "t0d2"),
		brill::FindDetectorConfig(config, "t0d3"),
		brill::FindDetectorConfig(config, "t0d4"),
		brill::FindDetectorConfig(config, "t0s")
	};
	for (int i = 0; i < kLayerCount; ++i) {
		if (!detectors[i]) {
			std::cerr << "Error: Missing detector config for " << DetectorNameFromLayer(i) << ".\n";
			return 1;
		}
	}

	const std::string trigger_infix = brill::TriggerInfix(config.trigger);
	const std::string track_dir = brill::JoinPath(config.workspace, config.paths.track);

	TChain chain("tree");
	int added_files = 0;
	for (int current_run = run; current_run <= end_run; ++current_run) {
		if (brill::IsJumpRun(config, current_run)) continue;
		std::string path = TString::Format(
			"%s/t0_%s%04d.root",
			track_dir.c_str(),
			trigger_infix.c_str(),
			current_run
		).Data();
		if (!std::filesystem::exists(path)) {
			std::cerr << "Warning: Skip missing tracked file " << path << ".\n";
			continue;
		}
		chain.Add(path.c_str());
		++added_files;
	}
	if (added_files == 0) {
		std::cout << "No tracked files to process.\n";
		return 0;
	}

	brill::T0Event event;
	brill::SetupInput(&chain, event);

	std::vector<FitSample> samples;
	std::vector<std::unique_ptr<brill::LostEnergyCalculator>> calculators;
	std::map<CalculatorKey, int> calculator_indices;

	const long long total = chain.GetEntries();
	long long last_percentage = -1;
	std::printf("Collecting T0 calibration samples   0%%");
	std::fflush(stdout);
	for (long long entry = 0; entry < total; ++entry) {
		long long percentage = total > 0 ? entry * 100ll / total : 100ll;
		if (percentage > last_percentage) {
			last_percentage = percentage;
			std::printf("\b\b\b\b%3lld%%", percentage);
			std::fflush(stdout);
		}
		chain.GetEntry(entry);
		for (int i = 0; i < event.num; ++i) {
			FitSample sample;
			if (!BuildSample(event, i, sample)) continue;
			sample.calculator = GetCalculatorIndex(
				config,
				event.charge[i],
				event.mass[i],
				sample.de_layer,
				*detectors[sample.de_layer],
				calculator_indices,
				calculators
			);
			samples.push_back(sample);
		}
	}
	std::printf("\b\b\b\b100%%\n");

	if (samples.empty()) {
		std::cerr << "Error: No calibration samples selected from tracked events.\n";
		return 1;
	}

	g_samples = &samples;
	g_calculators = &calculators;

	TMinuit minuit(kLayerCount * 2);
	minuit.SetPrintLevel(-1);
	minuit.SetFCN(CalibrationFcn);

	int ierflg = 0;
	double arglist[2] = {1.0, 0.0};
	minuit.mnexcm("SET ERR", arglist, 1, ierflg);

	const char *layer_names[kLayerCount] = {"d1", "d2", "d3", "d4", "s"};
	for (int layer = 0; layer < kLayerCount; ++layer) {
		minuit.DefineParameter(
			layer * 2,
			(std::string("p0_") + layer_names[layer]).c_str(),
			0.0,
			1.0,
			0.0,
			0.0
		);
		minuit.DefineParameter(
			layer * 2 + 1,
			(std::string("p1_") + layer_names[layer]).c_str(),
			0.006,
			1e-5,
			1e-8,
			1.0
		);
	}

	arglist[0] = 2000.0;
	arglist[1] = 1e-6;
	minuit.mnexcm("MIGRAD", arglist, 2, ierflg);
	if (ierflg != 0) {
		std::cerr << "Error: MIGRAD failed with code " << ierflg << ".\n";
		return 1;
	}

	double parameters[kLayerCount * 2] = {0.0};
	double errors[kLayerCount * 2] = {0.0};
	for (int index = 0; index < kLayerCount * 2; ++index) {
		minuit.GetParameter(index, parameters[index], errors[index]);
	}

	std::cout << "Calibration parameters:\n";
	for (int layer = 0; layer < kLayerCount; ++layer) {
		std::cout
			<< "  " << layer_names[layer]
			<< ": p0 = " << parameters[layer * 2]
			<< ", p1 = " << parameters[layer * 2 + 1]
			<< "\n";
	}

	std::string output_path = TString::Format(
		"%s/t0_%s%04d_%04d.txt",
		brill::JoinPath(config.workspace, config.paths.calibration).c_str(),
		trigger_infix.c_str(),
		run,
		end_run
	).Data();
	if (run == end_run) {
		output_path = TString::Format(
			"%s/t0_%s%04d.txt",
			brill::JoinPath(config.workspace, config.paths.calibration).c_str(),
			trigger_infix.c_str(),
			run
		).Data();
	}

	return WriteCalibrationParameters(output_path, parameters) == 0 ? 0 : 1;
}
