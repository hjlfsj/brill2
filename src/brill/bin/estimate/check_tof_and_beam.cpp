#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

#include <TCanvas.h>
#include <TCutG.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TString.h>
#include <TTree.h>

#include <catima/catima.h>

#include "external/cxxopts.hpp"
#include "include/config.h"
#include "include/energy_calculator/range_energy_calculator.h"
#include "include/event/beam/beam_sort.h"
#include "include/event/ingot/beam_event.h"
#include "include/event/t0/dssd_match_event.h"
#include "include/utils.h"

namespace {

constexpr int kCalibrationLayers = 5;
constexpr int kMaxEvents = 1000000;

constexpr double k14OMass = 13044.84;
constexpr double kLightSpeed = 299.792458;
constexpr double kTofDistance = 17000.0;

constexpr double kTargetDensity = 0.79867;

enum class TargetType { CD2, CH2 };

struct TargetConfig {
	TargetType type;
	double thickness_um;
	const char *name;
	const char *cache_suffix;
};

TargetConfig GetTargetConfig(int run) {
	if (run >= 57 && run <= 156) {
		return {TargetType::CD2, 121.7, "CD2", "cd2"};
	} else if (run >= 159 && run <= 183) {
		return {TargetType::CH2, 180.0, "CH2", "ch2"};
	} else {
		std::fprintf(
			stderr,
			"Warning: run %d outside known target ranges [57,156] and [159,183], "
			"defaulting to CD2 121.7 um.\n",
			run
		);
		return {TargetType::CD2, 121.7, "CD2", "cd2"};
	}
}

struct CalibrationParameters {
	double p0[kCalibrationLayers] = {0.0, 0.0, 0.0, 0.0, 0.0};
	double p1[kCalibrationLayers] = {1.0, 1.0, 1.0, 1.0, 1.0};
};

void PrintUsage(const cxxopts::Options &options) {
	std::cout << options.help() << "\n";
}

int ReadCalibrationParameters(
	const std::string &path,
	CalibrationParameters &parameters
) {
	std::ifstream fin(path);
	if (!fin.good()) {
		std::cerr << "Error: Open calibration parameter file " << path << " failed.\n";
		return -1;
	}

	std::string header;
	std::getline(fin, header);
	int index = -1;
	double p0 = 0.0;
	double p1 = 1.0;
	while (fin >> index >> p0 >> p1) {
		if (index < 0 || index >= kCalibrationLayers) continue;
		parameters.p0[index] = p0;
		parameters.p1[index] = p1;
	}
	return 0;
}

double CalibrateEnergy(
	const CalibrationParameters &parameters,
	int layer,
	double raw_energy
) {
	return parameters.p0[layer] + parameters.p1[layer] * raw_energy;
}

catima::Material CreateCD2Material() {
	catima::Material material;
	material.add_element(12.0, 6.0, 1.0);
	material.add_element(2.0, 1.0, 2.0);
	material.density(kTargetDensity);
	return material;
}

catima::Material CreateCH2Material() {
	catima::Material material;
	material.add_element(12.0, 6.0, 1.0);
	material.add_element(1.0, 1.0, 2.0);
	material.density(kTargetDensity);
	return material;
}

double CorrectBeamEnergy(
	double E0,
	brill::RangeEnergyCalculator &target_calc,
	double target_thickness_um
) {
	double range_target = target_calc.Range(E0);
	double range_at_target = range_target + target_thickness_um;
	return target_calc.Energy(range_at_target);
}

int OpenRunFiles(
	const std::string &match_dir,
	const std::string &trigger_infix,
	int run,
	TFile *&d1_file,
	TFile *&d2_file,
	TFile *&d3_file,
	TTree *&d1_tree,
	TTree *&d2_tree,
	TTree *&d3_tree,
	brill::DssdMatchEvent &d1_event,
	brill::DssdMatchEvent &d2_event,
	brill::DssdMatchEvent &d3_event
) {
	TString d1_path = TString::Format(
		"%s/t0d1_%s%04d.root",
		match_dir.c_str(),
		trigger_infix.c_str(),
		run
	);
	TString d2_path = TString::Format(
		"%s/t0d2_%s%04d.root",
		match_dir.c_str(),
		trigger_infix.c_str(),
		run
	);
	TString d3_path = TString::Format(
		"%s/t0d3_%s%04d.root",
		match_dir.c_str(),
		trigger_infix.c_str(),
		run
	);

	if (!std::filesystem::exists(d1_path.Data())) {
		std::cerr << "Error: d1 match file not found: " << d1_path << "\n";
		return 1;
	}
	if (!std::filesystem::exists(d2_path.Data())) {
		std::cerr << "Error: d2 match file not found: " << d2_path << "\n";
		return 1;
	}
	if (!std::filesystem::exists(d3_path.Data())) {
		std::cerr << "Error: d3 match file not found: " << d3_path << "\n";
		return 1;
	}

	d1_file = new TFile(d1_path, "read");
	d2_file = new TFile(d2_path, "read");
	d3_file = new TFile(d3_path, "read");

	d1_tree = static_cast<TTree*>(d1_file->Get("tree"));
	if (!d1_tree) {
		std::cerr << "Error: Get tree from " << d1_path << " failed.\n";
		return 1;
	}
	d2_tree = static_cast<TTree*>(d2_file->Get("tree"));
	if (!d2_tree) {
		std::cerr << "Error: Get tree from " << d2_path << " failed.\n";
		return 1;
	}
	d3_tree = static_cast<TTree*>(d3_file->Get("tree"));
	if (!d3_tree) {
		std::cerr << "Error: Get tree from " << d3_path << " failed.\n";
		return 1;
	}

	brill::SetupInput(d1_tree, d1_event);
	brill::SetupInput(d2_tree, d2_event);
	brill::SetupInput(d3_tree, d3_event);
	return 0;
}

TCutG *LoadBeamCut() {
	TString cut_path = "/home/ribll2026/ribll2026_www/github_code/brill2/src/brill/Cut/cal_d2_d3_14O_beam_cut.C";
	gROOT->ProcessLine(TString::Format(".x %s", cut_path.Data()));
	TCutG *beam_cut = static_cast<TCutG*>(gROOT->FindObject("cal_d2_d3_14O_beam_cut"));
	if (!beam_cut) {
		std::cerr << "Error: Load TCutG from " << cut_path << " failed.\n";
	}
	return beam_cut;
}

} // namespace

int main(int argc, char **argv) {
	cxxopts::Options options("check_tof_and_beam", "Check TOF and beam energy correlation.");
	options.add_options()
		("h,help", "Print help information.")
		(
			"c,config",
			"Config file path.",
			cxxopts::value<std::string>()->default_value("config.toml"),
			"file"
		)
		("r,run", "Run number.", cxxopts::value<int>(), "run")
		("L,distance",
			"TOF flight distance in mm (default: 17000).",
			cxxopts::value<double>()->default_value("17000.0"),
			"mm"
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
	const int run = result["run"].as<int>();
	config.trigger = "main";

	const double tof_distance = result["distance"].as<double>();
	const double light_speed = 299.792458;

	const std::string trigger_infix = brill::TriggerInfix(config.trigger);
	const std::string match_dir = brill::JoinPath(config.workspace, config.paths.match);
	const std::string ingot_dir = brill::JoinPath(config.workspace, config.paths.ingot);
	const std::string output_dir = brill::JoinPath(config.workspace, config.paths.beam);
	std::filesystem::create_directories(output_dir);

	TCutG *beam_cut = LoadBeamCut();
	if (!beam_cut) return 1;

	if (brill::IsJumpRun(config, run)) {
		std::cout << "Skipping jump run " << run << ".\n";
		return 0;
	}

	TFile *d1_file = nullptr;
	TFile *d2_file = nullptr;
	TFile *d3_file = nullptr;
	TTree *d1_tree = nullptr;
	TTree *d2_tree = nullptr;
	TTree *d3_tree = nullptr;
	brill::DssdMatchEvent d1_event;
	brill::DssdMatchEvent d2_event;
	brill::DssdMatchEvent d3_event;

	if (OpenRunFiles(match_dir, trigger_infix, run,
		d1_file, d2_file, d3_file,
		d1_tree, d2_tree, d3_tree,
		d1_event, d2_event, d3_event)) {
		return 1;
	}

	TString beam_path = TString::Format(
		"%s/beam_%s%04d.root",
		ingot_dir.c_str(),
		trigger_infix.c_str(),
		run
	);
	if (!std::filesystem::exists(beam_path.Data())) {
		std::cerr << "Error: beam file not found: " << beam_path << "\n";
		d1_file->Close();
		d2_file->Close();
		d3_file->Close();
		delete d1_file;
		delete d2_file;
		delete d3_file;
		return 1;
	}
	TFile beam_file(beam_path, "read");
	TTree *beam_tree = static_cast<TTree*>(beam_file.Get("tree"));
	if (!beam_tree) {
		std::cerr << "Error: Get tree from " << beam_path << " failed.\n";
		d1_file->Close();
		d2_file->Close();
		d3_file->Close();
		delete d1_file;
		delete d2_file;
		delete d3_file;
		return 1;
	}
	brill::BeamEvent beam_event;
	brill::SetupInput(beam_tree, beam_event);

	TString sorted_beam_path = TString::Format(
		"%s/beam_%s%04d.root",
		output_dir.c_str(),
		trigger_infix.c_str(),
		run
	);
	if (!std::filesystem::exists(sorted_beam_path.Data())) {
		std::cerr << "Error: sorted beam file not found: " << sorted_beam_path << "\n";
		d1_file->Close();
		d2_file->Close();
		d3_file->Close();
		delete d1_file;
		delete d2_file;
		delete d3_file;
		return 1;
	}
	TFile sorted_beam_file(sorted_beam_path, "read");
	TTree *sorted_beam_tree = static_cast<TTree*>(sorted_beam_file.Get("tree"));
	if (!sorted_beam_tree) {
		std::cerr << "Error: Get tree from " << sorted_beam_path << " failed.\n";
		d1_file->Close();
		d2_file->Close();
		d3_file->Close();
		delete d1_file;
		delete d2_file;
		delete d3_file;
		return 1;
	}
	bool v_14O = false;
	bool v_13N = false;
	bool v_12C = false;
	brill::SetupInputSortBeamTree(sorted_beam_tree, v_14O, v_13N, v_12C);

	int calib_run = brill::GetT0CalibrationRun(config, run);
	const std::string calibration_path = TString::Format(
		"%s/t0_%04d.txt",
		brill::JoinPath(config.workspace, config.paths.calibration).c_str(),
		calib_run
	).Data();

	CalibrationParameters calibration;
	if (ReadCalibrationParameters(calibration_path, calibration)) return 1;
	std::printf(
		"Calibration (t0_%04d.txt): "
		"d1(p0=%+.6f,p1=%.6f) d2(p0=%+.6f,p1=%.6f) d3(p0=%+.6f,p1=%.6f)\n",
		calib_run,
		calibration.p0[0], calibration.p1[0],
		calibration.p0[1], calibration.p1[1],
		calibration.p0[2], calibration.p1[2]
	);

	std::printf("Run %d: collecting up to %d events\n", run, kMaxEvents);
	std::printf("TOF distance: %.1f mm\n", tof_distance);
	std::printf("14O mass: %.2f MeV/c^2\n", k14OMass);

	const TargetConfig target_config = GetTargetConfig(run);
	std::printf(
		"Target: %s, thickness = %.1f um, density = %.5f g/cm^3\n",
		target_config.name, target_config.thickness_um, kTargetDensity
	);

	const double mass = k14OMass;

	const std::string cache_dir = brill::JoinPath(config.workspace, config.paths.energy_calculator);
	TString target_cache_path = TString::Format(
		"%s/target_z8_a14_%s.root",
		cache_dir.c_str(),
		target_config.cache_suffix
	);

	catima::Material target_material;
	if (target_config.type == TargetType::CH2) {
		target_material = CreateCH2Material();
	} else {
		target_material = CreateCD2Material();
	}

	brill::RangeEnergyCalculator target_calc(8, 14, target_material, target_cache_path.Data());
	std::printf(
		"Energy loss correction: %s target = %.1f um\n",
		target_config.name, target_config.thickness_um
	);

	TH1D *h_beam_E = new TH1D("h_beam_E", "Beam Energy;E (MeV);Counts",
		1000, 0.0, 1000.0);
	h_beam_E->SetDirectory(0);

	TH1D *h_tof = new TH1D("h_tof", "TOF;TOF (ns);Counts",
		2000, 0.0, 200.0);
	h_tof->SetDirectory(0);

	TH2D *h_E_vs_tof = new TH2D("h_E_vs_tof", "E vs TOF;TOF (ns);E (MeV)",
		400, 50.0, 70.0,
		500, 450.0, 550.0);
	h_E_vs_tof->SetDirectory(0);

	std::vector<double> E_vec;
	std::vector<double> tof_vec;
	E_vec.reserve(kMaxEvents);
	tof_vec.reserve(kMaxEvents);

	Long64_t total = d1_tree->GetEntries();
	Long64_t passed = 0;
	Long64_t last_percentage = -1;

	std::printf("Processing events   0%%");
	std::fflush(stdout);

	for (Long64_t entry = 0; entry < total && passed < kMaxEvents; ++entry) {
		Long64_t percentage = total > 0 ? entry * 100ll / total : 100ll;
		if (percentage > last_percentage) {
			last_percentage = percentage;
			std::printf("\b\b\b\b%3lld%%", percentage);
			std::fflush(stdout);
		}

		if (d1_tree->GetEntry(entry) <= 0
			|| d2_tree->GetEntry(entry) <= 0
			|| d3_tree->GetEntry(entry) <= 0) continue;

		if (d1_event.num != 1 || d2_event.num != 1 || d3_event.num != 1) continue;

		double dx = d3_event.x[0] - d1_event.x[0];
		double dy = d3_event.y[0] - d1_event.y[0];
		if (dx * dx + dy * dy >= 4.0) continue;

		double d1_cal_e = CalibrateEnergy(calibration, 0, d1_event.energy[0]);
		double d2_cal_e = CalibrateEnergy(calibration, 1, d2_event.energy[0]);
		double d3_cal_e = CalibrateEnergy(calibration, 2, d3_event.energy[0]);

		if (!beam_cut->IsInside(d3_cal_e, d2_cal_e)) continue;

		if (entry < beam_tree->GetEntries()) {
			beam_tree->GetEntry(entry);
		}
		if (entry < sorted_beam_tree->GetEntries()) {
			sorted_beam_tree->GetEntry(entry);
		}

		if (!beam_event.valid || !v_14O) continue;

		double E_meas = d1_cal_e + d2_cal_e + d3_cal_e;
		double E = CorrectBeamEnergy(E_meas, target_calc, target_config.thickness_um);
		double tof = beam_event.tof;

		h_beam_E->Fill(E);
		h_tof->Fill(tof);
		h_E_vs_tof->Fill(tof, E);
		E_vec.push_back(E);
		tof_vec.push_back(tof);

		++passed;
	}
	std::printf("\b\b\b\b100%%\n");
	std::printf("Passed events: %lld\n", passed);

	d1_file->Close();
	d2_file->Close();
	d3_file->Close();
	delete d1_file;
	delete d2_file;
	delete d3_file;
	beam_file.Close();
	sorted_beam_file.Close();

	std::printf("\nFitting beam energy spectrum...\n");
	int max_bin_E = h_beam_E->GetMaximumBin();
	double peak_pos_E = h_beam_E->GetBinCenter(max_bin_E);
	double fit_low_E = peak_pos_E - 6.0;
	double fit_high_E = peak_pos_E + 6.0;
	std::printf("  E peak bin = %d, center = %.2f MeV, fit range [%.2f, %.2f]\n",
		max_bin_E, peak_pos_E, fit_low_E, fit_high_E);
	int fit_status_E = h_beam_E->Fit("gaus", "Q", "", fit_low_E, fit_high_E);
	double E0 = 0.0;
	double E0_err = 0.0;
	double sigma_E = 0.0;
	if (fit_status_E == 0) {
		TF1 *fit_E = h_beam_E->GetFunction("gaus");
		if (fit_E) {
			E0 = fit_E->GetParameter(1);
			E0_err = fit_E->GetParError(1);
			sigma_E = std::abs(fit_E->GetParameter(2));
			std::printf("  E0 = %.4f +/- %.4f MeV, sigma = %.4f MeV\n",
				E0, E0_err, sigma_E);
		}
	} else {
		std::printf("  Beam energy fit failed (status=%d)\n", fit_status_E);
	}

	std::printf("\nFitting TOF spectrum...\n");
	int max_bin_tof = h_tof->GetMaximumBin();
	double peak_pos_tof = h_tof->GetBinCenter(max_bin_tof);
	double fit_low_tof = peak_pos_tof - 1.0;
	double fit_high_tof = peak_pos_tof + 1.0;
	std::printf("  TOF peak bin = %d, center = %.2f ns, fit range [%.2f, %.2f]\n",
		max_bin_tof, peak_pos_tof, fit_low_tof, fit_high_tof);
	int fit_status_tof = h_tof->Fit("gaus", "Q", "", fit_low_tof, fit_high_tof);
	double t0 = 0.0;
	double t0_err = 0.0;
	double sigma_tof = 0.0;
	if (fit_status_tof == 0) {
		TF1 *fit_tof = h_tof->GetFunction("gaus");
		if (fit_tof) {
			t0 = fit_tof->GetParameter(1);
			t0_err = fit_tof->GetParError(1);
			sigma_tof = std::abs(fit_tof->GetParameter(2));
			std::printf("  t0 = %.4f +/- %.4f ns, sigma = %.4f ns\n",
				t0, t0_err, sigma_tof);
		}
	} else {
		std::printf("  TOF fit failed (status=%d)\n", fit_status_tof);
	}

	double dt = 0.0;
	if (E0 > 0.0 && t0 > 0.0) {
		double gamma = 1.0 + E0 / mass;
		double beta = std::sqrt(1.0 - 1.0 / (gamma * gamma));
		double v = beta * light_speed;
		double t_true = tof_distance / v;
		dt = t0 - t_true;
		std::printf("\nRelativistic kinetics (14O, m=%.2f MeV/c^2, L=%.1f mm):\n",
			mass, tof_distance);
		std::printf("  E0 = %.4f MeV -> gamma = %.6f, beta = %.6f\n",
			E0, gamma, beta);
		std::printf("  True TOF = %.4f ns, measured TOF = %.4f ns\n",
			t_true, t0);
		std::printf("  TOF offset dt = t_measured - t_true = %.4f ns\n", dt);
	}

	// ==================== Gated histograms ====================

	TH1D *h_beam_E_tof_gated = new TH1D("h_beam_E_tof_gated",
		"Beam E (TOF-gated);E (MeV);Counts",
		1000, 0.0, 1000.0);
	h_beam_E_tof_gated->SetDirectory(0);

	TH1D *h_tof_E_gated = new TH1D("h_tof_E_gated",
		"TOF (E-gated);TOF (ns);Counts",
		2000, 0.0, 200.0);
	h_tof_E_gated->SetDirectory(0);

	if (sigma_tof > 0.0) {
		double tof_low = t0 - 0.25 * sigma_tof;
		double tof_high = t0 + 0.25 * sigma_tof;
		std::printf("\nTOF-gated E: selecting tof in [%.4f, %.4f] ns (t0 +/- 0.25*sigma_tof)\n",
			tof_low, tof_high);
		for (size_t i = 0; i < E_vec.size(); ++i) {
			if (tof_vec[i] >= tof_low && tof_vec[i] <= tof_high) {
				h_beam_E_tof_gated->Fill(E_vec[i]);
			}
		}
		std::printf("  Entries in TOF gate: %.0f\n", h_beam_E_tof_gated->GetEntries());

		int max_bin_E_gated = h_beam_E_tof_gated->GetMaximumBin();
		double peak_pos_E_gated = h_beam_E_tof_gated->GetBinCenter(max_bin_E_gated);
		double fit_low_E_gated = peak_pos_E_gated - 6.0;
		double fit_high_E_gated = peak_pos_E_gated + 6.0;
		int fit_status_E_gated = h_beam_E_tof_gated->Fit("gaus", "Q", "",
			fit_low_E_gated, fit_high_E_gated);
		if (fit_status_E_gated == 0) {
			TF1 *fit_E_gated = h_beam_E_tof_gated->GetFunction("gaus");
			if (fit_E_gated) {
				double E0_gated = fit_E_gated->GetParameter(1);
				double E0_gated_err = fit_E_gated->GetParError(1);
				double sigma_E_gated = std::abs(fit_E_gated->GetParameter(2));
				std::printf("  TOF-gated E0 = %.4f +/- %.4f MeV, sigma = %.4f MeV\n",
					E0_gated, E0_gated_err, sigma_E_gated);
			}
		} else {
			std::printf("  TOF-gated E fit failed (status=%d)\n", fit_status_E_gated);
		}
	}

	if (sigma_E > 0.0) {
		double E_low = E0 - 0.25 * sigma_E;
		double E_high = E0 + 0.25 * sigma_E;
		std::printf("\nE-gated TOF: selecting E in [%.4f, %.4f] MeV (E0 +/- 0.25*sigma_E)\n",
			E_low, E_high);
		for (size_t i = 0; i < tof_vec.size(); ++i) {
			if (E_vec[i] >= E_low && E_vec[i] <= E_high) {
				h_tof_E_gated->Fill(tof_vec[i]);
			}
		}
		std::printf("  Entries in E gate: %.0f\n", h_tof_E_gated->GetEntries());

		int max_bin_tof_gated = h_tof_E_gated->GetMaximumBin();
		double peak_pos_tof_gated = h_tof_E_gated->GetBinCenter(max_bin_tof_gated);
		double fit_low_tof_gated = peak_pos_tof_gated - 1.0;
		double fit_high_tof_gated = peak_pos_tof_gated + 1.0;
		int fit_status_tof_gated = h_tof_E_gated->Fit("gaus", "Q", "",
			fit_low_tof_gated, fit_high_tof_gated);
		if (fit_status_tof_gated == 0) {
			TF1 *fit_tof_gated = h_tof_E_gated->GetFunction("gaus");
			if (fit_tof_gated) {
				double sigma_tof_gated = std::abs(fit_tof_gated->GetParameter(2));
				std::printf("  E-gated TOF sigma = %.4f ns\n", sigma_tof_gated);
			}
		} else {
			std::printf("  E-gated TOF fit failed (status=%d)\n", fit_status_tof_gated);
		}
	}

	TString output_path = TString::Format(
		"%s/check_tof_and_beam_%04d.root",
		output_dir.c_str(),
		run
	);
	TFile output_file(output_path, "recreate");

	TCanvas c_E("c_E", "Beam Energy Spectrum", 800, 600);
	h_beam_E->Draw();
	output_file.cd();
	c_E.Write();
	h_beam_E->Write();

	TCanvas c_tof("c_tof", "TOF Spectrum", 800, 600);
	h_tof->Draw();
	output_file.cd();
	c_tof.Write();
	h_tof->Write();

	TCanvas c_2d("c_2d", "E vs TOF", 800, 600);
	h_E_vs_tof->Draw("colz");

	if (dt != 0.0) {
		auto f_theory = [=](double *x, double *) -> double {
			double t_meas = x[0];
			double t_corr = t_meas - dt;
			if (t_corr <= 0.0) return 0.0;
			double v = tof_distance / t_corr;
			double beta = v / light_speed;
			double beta2 = beta * beta;
			if (beta2 >= 1.0) return 0.0;
			double gamma = 1.0 / std::sqrt(1.0 - beta2);
			return mass * (gamma - 1.0);
		};

		TF1 *theory_curve = new TF1("theory_curve", f_theory, 50.0, 70.0, 0);
		theory_curve->SetLineColor(kRed);
		theory_curve->SetLineWidth(2);
		theory_curve->Draw("same");

		TPaveText *pt = new TPaveText(0.15, 0.75, 0.45, 0.88, "NDC");
		pt->SetFillStyle(0);
		pt->SetBorderSize(1);
		pt->AddText(TString::Format("E_{0} = %.2f #pm %.2f MeV", E0, E0_err));
		pt->AddText(TString::Format("t_{0} = %.2f #pm %.2f ns", t0, t0_err));
		pt->AddText(TString::Format("#Deltat = %.3f ns", dt));
		pt->Draw();
	}

	output_file.cd();
	c_2d.Write();
	h_E_vs_tof->Write();

	TCanvas c_E_tof_gated("c_E_tof_gated", "Beam E (TOF-gated)", 800, 600);
	h_beam_E_tof_gated->Draw();
	output_file.cd();
	c_E_tof_gated.Write();
	h_beam_E_tof_gated->Write();

	TCanvas c_tof_E_gated("c_tof_E_gated", "TOF (E-gated)", 800, 600);
	h_tof_E_gated->Draw();
	output_file.cd();
	c_tof_E_gated.Write();
	h_tof_E_gated->Write();

	output_file.Close();

	std::printf("\nOutput saved to: %s\n", output_path.Data());

	TString txt_path = TString::Format(
		"%s/check_tof_and_beam_%04d.txt",
		output_dir.c_str(),
		run
	);
	std::ofstream ofs(txt_path.Data());
	if (ofs.is_open()) {
		ofs << "# E0        sigma_E0   t0         sigma_t0   dt\n";
		ofs << TString::Format("%.4f   %.4f     %.4f     %.4f      %.6f\n",
			E0, sigma_E, t0, sigma_tof, dt);
		ofs.close();
		std::printf("Parameters saved to: %s\n", txt_path.Data());
	}

	return 0;
}