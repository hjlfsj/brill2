#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TCutG.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH2F.h>
#include <TROOT.h>
#include <TString.h>

#include "external/cxxopts.hpp"
#include "include/config.h"
#include "include/energy_calculator/delta_energy_calculator.h"
#include "include/energy_calculator/range_energy_calculator.h"
#include "include/utils.h"

namespace {

constexpr int kLayerCount = 5;

struct LayerPairInfo {
	const char *name;
	int layer;
	const char *graph_name;
};

const LayerPairInfo kLayerPairs[] = {
	{"d1d2", 0, "g_d1d2"},
	{"d2d3", 1, "g_d2d3"},
	{"d3d4", 2, "g_d3d4"},
	{"d4s",  3, "g_d4s"},
};

struct ParticlePidInfo {
	int layer;
	int charge;
	int mass;
	double left;
	double right;
	double offset;
};

const std::vector<ParticlePidInfo> pid_info {
	{0, 2,  4,  2000.0, 11000.0,   	   0},  // d1d2 4He
	// {0, 3,  6,  4500.0, 15000.0,   12000+1000},  // d1d2 6Li
	{1, 2,  4,  3900.0,  8000.0,   29000+1000},  // d2d3 4He
	{1, 4,  7, 11500.0, 19000.0,   38000+1000},  // d2d3 7Be
	{1, 6, 12, 22500.0, 46000.0,   58000+1000},  // d2d3 12C
	{2, 2,  4,  3400.0,  7000.0,  105000+1000},  // d3d4 4He
	{2, 4,  7,  9500.0, 19000.0,  113000+1000},  // d3d4 7Be
	{2, 6, 12, 22000.0, 40000.0,  133000+1000},  // d3d4 12C
	{3, 1,  1,  1000.0,  2300.0,  174000+1000},  // d4s  1H
	{3, 2,  4,  4000.0, 10000.0,  178000+1000},  // d4s  4He
	{3, 4,  7, 12000.0, 28000.0,  190000+1000},  // d4s  7Be
};

class PidFitFunc {
public:
	PidFitFunc(
		const std::vector<std::pair<int, int>> &projectiles,
		const brill::AppConfig &config
	) {
		for (const auto &p : projectiles) {
			calculators_.insert(
				std::make_pair(
					p.first * 100 + p.second,
					std::make_unique<brill::DeltaEnergyCalculator>(
						config,
						p.first,
						p.second
					)
				)
			);
		}
	}

	double operator()(double *x, double *par) const {
		for (const auto &info : pid_info) {
			if (
				x[0] > info.left + info.offset
				&& x[0] < info.right + info.offset
			) {
				double de =
					par[info.layer*2]
					+ par[info.layer*2+1] * (x[0] - info.offset);
				auto calculator = calculators_.at(info.charge*100+info.mass);
				double e = calculator->Energy(info.layer, de);
				return (e - par[info.layer*2+2]) / par[info.layer*2+3];
			}
		}
		return 0.0;
	}
private:
	std::map<int, std::shared_ptr<brill::DeltaEnergyCalculator>> calculators_;
};

std::string layer_names[kLayerCount] = {
	"t0d1", "t0d2", "t0d3", "t0d4", "t0s"
};

struct ParticleIdentity {
	int charge;
	int mass;
};

ParticleIdentity ParseParticleName(const std::string &name) {
	static const std::map<std::string, ParticleIdentity> kParticleMap = {
		{"1H",  {1, 1}},
		{"4He", {2, 4}},
		{"6Li", {3, 6}},
		{"7Be", {4, 7}},
		{"12C", {6, 12}},
	};
	auto it = kParticleMap.find(name);
	if (it != kParticleMap.end()) return it->second;
	return {0, 0};
}

void PrintUsage(const cxxopts::Options &options) {
	std::cout << options.help() << "\n";
}

void WriteCalibrationParameters(
	const char *path,
	const double *parameters
) {
	std::filesystem::path file_path(path);
	if (!file_path.parent_path().empty()) {
		std::filesystem::create_directories(file_path.parent_path());
	}
	std::ofstream fout(path);
	if (!fout.good()) {
		std::cerr << "Error: Open output parameter file " << path << " failed.\n";
		return;
	}

	fout << "# layer p0 p1\n";
	for (int layer = 0; layer < kLayerCount; ++layer) {
		fout << layer << " " << parameters[layer * 2] << " " << parameters[layer * 2 + 1] << "\n";
	}
}

struct LoadedCut {
	const LayerPairInfo *pair;
	ParticleIdentity particle;
	std::unique_ptr<TCutG> cut;
};

std::string FindPreCalibrationFile(
	const std::string &estimate_dir,
	const std::string &prefix,
	int run
) {
	for (const auto &entry : std::filesystem::directory_iterator(estimate_dir)) {
		if (!entry.is_regular_file()) continue;
		std::string name = entry.path().filename().string();
		char expected[256];
		std::snprintf(expected, sizeof(expected), "pre_calibration_%s%04d_", prefix.c_str(), run);
		if (name.find(expected) == 0 && name.find(".root") != std::string::npos) {
			return entry.path().string();
		}
	}
	return "";
}

int LoadCuts(
	const std::string &cut_dir,
	const std::string &cut_suffix,
	std::vector<LoadedCut> &cuts
) {
	for (const auto &pair : kLayerPairs) {
		for (const auto &entry : std::filesystem::directory_iterator(cut_dir)) {
			if (!entry.is_regular_file()) continue;
			std::string name = entry.path().filename().string();
			if (name.size() <= 2 || name.substr(name.size()-2) != ".C") continue;
			std::string basename = name.substr(0, name.size()-2);
			std::string expected_prefix = std::string(pair.name) + "_";
			if (basename.find(expected_prefix) != 0) continue;
			if (basename.find(cut_suffix) == std::string::npos) continue;

			std::string particle_name = basename.substr(
				expected_prefix.size(),
				basename.size() - expected_prefix.size() - cut_suffix.size()
			);
			ParticleIdentity pid = ParseParticleName(particle_name);
			if (pid.charge == 0 && pid.mass == 0) {
				std::cerr << "Warning: Unknown particle '" << particle_name
					<< "' in cut file " << name << ", skipping.\n";
				continue;
			}

			TString cut_object_name = TString::Format("%s_%s%s",
				pair.name, particle_name.c_str(), cut_suffix.c_str());
			std::string macro_path = entry.path().string();
			gROOT->ProcessLine(TString::Format(".x %s", macro_path.c_str()));
			TCutG *cutg = dynamic_cast<TCutG*>(gROOT->FindObject(cut_object_name));
			if (!cutg) {
				std::cerr << "Error: Load TCutG " << cut_object_name
					<< " from " << macro_path << " failed.\n";
				return -1;
			}

			LoadedCut loaded;
			loaded.pair = &pair;
			loaded.particle = pid;
			loaded.cut.reset(dynamic_cast<TCutG*>(cutg->Clone()));
			cuts.push_back(std::move(loaded));
			std::cout << "Loaded cut: " << name << " (Z=" << pid.charge
				<< ", A=" << pid.mass << ")\n";
		}
	}
	return cuts.empty() ? -1 : 0;
}

const ParticlePidInfo *GetPidInfo(int layer, int charge, int mass) {
	for (const auto &info : pid_info) {
		if (info.layer == layer && info.charge == charge
			&& info.mass == mass) {
			return &info;
		}
	}
	return nullptr;
}

const std::vector<ParticleIdentity> kCommonParticles = {
	{1, 1}, {1, 2}, {1, 3},
	{2, 3}, {2, 4}, {2, 6},
	{3, 6}, {3, 7}, {3, 8}, {3, 9},
	{4, 7}, {4, 9}, {4, 10},
	{5, 10}, {5, 11}, {5, 12},
	{6, 10}, {6, 11}, {6, 12}, {6, 13},
	{7, 12}, {7, 13}, {7, 14},
	{8, 13}, {8, 14},
};

TGraph *GenerateTheoryCurve(
	const brill::RangeEnergyCalculator &calculator,
	double t_shallow,
	double t_deep
) {
	auto *curve = new TGraph();
	curve->SetLineColor(kRed);
	curve->SetLineWidth(1);

	double E_min = calculator.Energy(t_shallow);
	double E_max = calculator.Energy(t_shallow + t_deep);

	constexpr double kStep = 0.1;
	for (double E = E_min; E <= E_max; E += kStep) {
		double r_residual = calculator.Range(E) - t_shallow;
		if (r_residual <= 0.0) continue;
		if (r_residual > t_deep) break;

		double E_residual = calculator.Energy(r_residual);
		double delta_E = E - E_residual;
		if (delta_E < 0.0 || E_residual < 0.0) continue;

		curve->SetPoint(curve->GetN(), E_residual, delta_E);
	}

	return curve;
}

} // namespace

int main(int argc, char **argv) {
	cxxopts::Options options("calibrate_t0",
		"Calibrate T0 energies using pre_calibration TGraph and TCutG cuts.");
	options.add_options()
		("h,help", "Print help information.")
		("r,run", "Start run number.", cxxopts::value<int>(), "run")
		("t,trigger", "Trigger type.", cxxopts::value<std::string>(), "trigger")
		("c,config", "Config file path.",
			cxxopts::value<std::string>()->default_value("config.toml"), "file");

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

	const brill::SquareDetectorConfig *detectors[kLayerCount] = {
		brill::FindDetectorConfig(config, "t0d1"),
		brill::FindDetectorConfig(config, "t0d2"),
		brill::FindDetectorConfig(config, "t0d3"),
		brill::FindDetectorConfig(config, "t0d4"),
		brill::FindDetectorConfig(config, "t0s")
	};
	for (int i = 0; i < kLayerCount; ++i) {
		if (!detectors[i]) {
			std::cerr << "Error: Missing detector config for "
				<< layer_names[i] << ".\n";
			return 1;
		}
	}

	const std::string trigger_infix = brill::TriggerInfix(config.trigger);

	const std::string estimate_dir = brill::JoinPath(
		config.workspace, config.paths.estimate);
	std::string precal_path = FindPreCalibrationFile(
		estimate_dir, trigger_infix, run);
	if (precal_path.empty()) {
		std::cerr << "Error: No pre_calibration file found for run "
			<< run << " in " << estimate_dir << ".\n";
		return 1;
	}
	std::cout << "Reading pre_calibration: " << precal_path << "\n";

	TFile precal_file(precal_path.c_str(), "read");
	if (precal_file.IsZombie()) {
		std::cerr << "Error: Open pre_calibration file "
			<< precal_path << " failed.\n";
		return 1;
	}

	TGraph *graphs[4] = {nullptr};
	const char *graph_names[4] = {"g_d1d2", "g_d2d3", "g_d3d4", "g_d4s"};
	for (int i = 0; i < 4; ++i) {
		graphs[i] = dynamic_cast<TGraph*>(precal_file.Get(graph_names[i]));
		if (!graphs[i]) {
			std::cerr << "Error: TGraph " << graph_names[i]
				<< " not found in " << precal_path << ".\n";
			return 1;
		}
		std::cout << "  " << graph_names[i] << ": "
			<< graphs[i]->GetN() << " points\n";
	}

	const std::string cut_dir = "src/brill/Cut";
	if (!std::filesystem::exists(cut_dir)) {
		std::cerr << "Error: Cut directory " << cut_dir << " not found.\n";
		return 1;
	}
	std::cout << "\nLoading cuts from: " << cut_dir << "\n";

	std::vector<LoadedCut> cuts;
	if (LoadCuts(cut_dir, "_stop", cuts)) {
		std::cerr << "Error: No cuts loaded from " << cut_dir << ".\n";
		return 1;
	}

	TString output_path = TString::Format(
		"%s/t0_%s%04d.root",
		brill::JoinPath(config.workspace, config.paths.calibration).c_str(),
		trigger_infix.c_str(),
		run
	);
	TFile opf(output_path, "recreate");
	std::cout << "Output ROOT file: " << output_path << "\n";
	TGraph gcali;

	std::set<std::pair<int, int>> projectile_set;
	for (const auto &cut : cuts) {
		int layer = cut.pair->layer;
		const ParticlePidInfo *info = GetPidInfo(
			layer, cut.particle.charge, cut.particle.mass);
		if (!info) {
			std::cerr << "Error: No pid_info for (layer=" << layer
				<< ", Z=" << cut.particle.charge
				<< ", A=" << cut.particle.mass << ").\n";
			continue;
		}

		const TGraph *g = graphs[layer];
		int hit_count = 0;
		for (int pt = 0; pt < g->GetN(); ++pt) {
			double deep = g->GetPointX(pt);
			double shallow = g->GetPointY(pt);
			if (cut.cut->IsInside(deep, shallow)
				&& shallow > info->left
				&& shallow < info->right) {
				gcali.AddPoint(shallow + info->offset, deep);
				++hit_count;
			}
		}
		std::cout << "  " << cut.pair->name << " " << cut.particle.charge
			<< "/" << cut.particle.mass << ": " << hit_count
			<< " points selected\n";
		projectile_set.insert({cut.particle.charge, cut.particle.mass});
	}

	if (gcali.GetN() == 0) {
		std::cerr << "Error: No calibration points after cut selection.\n";
		return 1;
	}
	std::cout << "Total calibration points: " << gcali.GetN() << "\n";

	std::vector<std::pair<int, int>> projectiles(
		projectile_set.begin(), projectile_set.end());

	double initial_calibration_parameters[10] = {
		0.0, 0.002,
		0.0, 0.006,
		0.0, 0.006,
		0.0, 0.003,
		0.0, 0.003
	};

	PidFitFunc pid_fit(projectiles, config);
	TF1 fcali("fcali", pid_fit, 0.0, 330'000.0, 10);
	fcali.SetNpx(10000);
	for (size_t i = 0; i < 10; ++i) {
		fcali.SetParameter(i, initial_calibration_parameters[i]);
	}
	for (int i = 0; i < 10; ++i) {
		fcali.SetParLimits(i, 0.0, (i % 2 == 0) ? 10.0 : 1.0);
	}
	std::cout << "Fitting with \"R S\" option...\n";
	gcali.Fit(&fcali, "R S");

	std::cout << "Fit result:\n"
		<< "  chi2     = " << fcali.GetChisquare() << "\n"
		<< "  ndf      = " << fcali.GetNDF() << "\n"
		<< "  chi2/ndf = " << fcali.GetChisquare() / fcali.GetNDF() << "\n"
		<< "  prob     = " << fcali.GetProb() << "\n";

	double *parameters = fcali.GetParameters();
	std::cout << "Calibration parameters:\n";
	for (int layer = 0; layer < kLayerCount; ++layer) {
		std::cout
			<< "  " << layer_names[layer]
			<< ": p0 = " << parameters[layer * 2]
			<< ", p1 = " << parameters[layer * 2 + 1]
			<< "\n";
	}

	TString txt_path = TString::Format(
		"%s/t0_%04d.txt",
		brill::JoinPath(config.workspace, config.paths.calibration).c_str(),
		run
	);
	WriteCalibrationParameters(txt_path.Data(), parameters);
	std::cout << "Calibration written to: " << txt_path << "\n";

	double thickness[kLayerCount];
	for (int i = 0; i < kLayerCount; ++i) {
		thickness[i] = detectors[i]->thickness_um;
	}

	const char *th2_names[4] = {
		"d1d2", "d2d3", "d3d4", "d4s"};
	const char *th2_cal_names[4] = {
		"d1d2_cal", "d2d3_cal", "d3d4_cal", "d4s_cal"};

	std::cout << "\nCalibrating TH2Fs and drawing theory curves...\n";

	for (int i = 0; i < 4; ++i) {
		TH2F *h2_raw = dynamic_cast<TH2F*>(
			precal_file.Get(th2_names[i]));
		if (!h2_raw) {
			std::cerr << "Warning: TH2F " << th2_names[i]
				<< " not found, skipping.\n";
			continue;
		}

		int deep_layer = i + 1;
		int shallow_layer = i;

		double p0_deep = parameters[deep_layer * 2];
		double p1_deep = parameters[deep_layer * 2 + 1];
		double p0_shallow = parameters[shallow_layer * 2];
		double p1_shallow = parameters[shallow_layer * 2 + 1];

		double x_min_raw = h2_raw->GetXaxis()->GetXmin();
		double x_max_raw = h2_raw->GetXaxis()->GetXmax();
		double y_min_raw = h2_raw->GetYaxis()->GetXmin();
		double y_max_raw = h2_raw->GetYaxis()->GetXmax();

		double x_min_cal = p0_deep + p1_deep * x_min_raw;
		double x_max_cal = p0_deep + p1_deep * x_max_raw;
		double y_min_cal =
			p0_shallow + p1_shallow * y_min_raw;
		double y_max_cal =
			p0_shallow + p1_shallow * y_max_raw;

		int n_bins_x = h2_raw->GetNbinsX();
		int n_bins_y = h2_raw->GetNbinsY();

		auto *h2_cal = new TH2F(
			th2_cal_names[i],
			TString::Format("%s Calibrated;%s [MeV];%s [MeV]",
				kLayerPairs[i].name,
				layer_names[deep_layer].c_str(),
				layer_names[shallow_layer].c_str()),
			n_bins_x, x_min_cal, x_max_cal,
			n_bins_y, y_min_cal, y_max_cal);

		for (int bx = 1; bx <= n_bins_x; ++bx) {
			for (int by = 1; by <= n_bins_y; ++by) {
				double content =
					h2_raw->GetBinContent(bx, by);
				if (content > 0.0) {
					h2_cal->SetBinContent(bx, by, content);
				}
			}
		}
		h2_cal->SetEntries(h2_raw->GetEntries());

		opf.cd();
		h2_cal->Write();

		TCanvas *c = new TCanvas(
			TString::Format("c_%s", th2_names[i]),
			TString::Format("%s Calibrated with Theory",
				kLayerPairs[i].name),
			1200, 900);
		c->cd();
		h2_cal->Draw("colz");

		std::vector<TGraph*> theory_curves;
		for (const auto &p : kCommonParticles) {
			std::string cache_path = TString::Format(
				"%s/si_z%d_a%d.root",
				brill::JoinPath(
					config.workspace,
					config.paths.energy_calculator).c_str(),
				p.charge, p.mass).Data();

			try {
				brill::RangeEnergyCalculator rec(
					p.charge, p.mass,
					brill::SiliconMaterial(),
					cache_path);

				TGraph *curve = GenerateTheoryCurve(
					rec,
					thickness[shallow_layer],
					thickness[deep_layer]);
				if (curve->GetN() > 0) {
					curve->SetLineColor(kRed);
					curve->SetLineWidth(1);
					curve->Draw("l same");
					theory_curves.push_back(curve);
				} else {
					delete curve;
				}
			} catch (...) {
				continue;
			}
		}

		opf.cd();
		c->Write();

		delete c;
	}

	opf.cd();
	gcali.Write("gcali");
	fcali.Write("fcali");
	opf.Close();
	return 0;
}