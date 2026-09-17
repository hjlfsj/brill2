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
#include "include/t0/calibrate_t0_utils.h"
#include "include/utils.h"

namespace {

void PrintUsage(const cxxopts::Options &options) {
	std::cout << options.help() << "\n";
}

class PidFitFuncStage1 {
public:
	PidFitFuncStage1(
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
		for (const auto &info : brill::kT0PidInfo) {
			if (info.layer == 0) continue;
			if (
				x[0] > info.left + info.offset
				&& x[0] < info.right + info.offset
			) {
				int upper_idx = (info.layer - 1) * 2;
				int lower_idx = info.layer * 2;
				double de =
					par[upper_idx]
					+ par[upper_idx + 1]
						* (x[0] - info.offset);
				auto it = calculators_.find(
					info.charge * 100 + info.mass);
				if (it == calculators_.end()) return 0.0;
				double e = it->second->Energy(
					info.layer, de);
				return (e - par[lower_idx])
					/ par[lower_idx + 1];
			}
		}
		return 0.0;
	}

private:
	std::map<int,
		std::shared_ptr<brill::DeltaEnergyCalculator>>
		calculators_;
};

class PidFitFuncStage2 {
public:
	PidFitFuncStage2(
		const std::vector<std::pair<int, int>> &projectiles,
		const brill::AppConfig &config,
		double d2_p0,
		double d2_p1
	) : d2_p0_(d2_p0), d2_p1_(d2_p1) {
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
		for (const auto &info : brill::kT0PidInfo) {
			if (info.layer != 0) continue;
			if (
				x[0] > info.left + info.offset
				&& x[0] < info.right + info.offset
			) {
				double de =
					par[0]
					+ par[1] * (x[0] - info.offset);
				auto it = calculators_.find(
					info.charge * 100 + info.mass);
				if (it == calculators_.end()) return 0.0;
				double e = it->second->Energy(
					info.layer, de);
				return (e - d2_p0_) / d2_p1_;
			}
		}
		return 0.0;
	}

private:
	double d2_p0_;
	double d2_p1_;
	std::map<int,
		std::shared_ptr<brill::DeltaEnergyCalculator>>
		calculators_;
};

int RebuildTheoryCurves(const brill::AppConfig &config) {
	std::string ec_dir = brill::JoinPath(
		config.workspace,
		config.paths.energy_calculator);
	std::cout << "Removing cached delta-E theory curves in "
		<< ec_dir << "...\n";

	int removed = 0;
	for (const auto &entry :
		std::filesystem::directory_iterator(ec_dir)) {
		if (!entry.is_regular_file()) continue;
		std::string name = entry.path().filename().string();
		if (name.find("t0_delta_z") == 0
			&& name.find(".root") != std::string::npos) {
			std::filesystem::remove(entry.path());
			std::cout << "  Deleted: " << name << "\n";
			++removed;
		}
	}

	if (removed == 0) {
		std::cout << "  No cached files found.\n";
	} else {
		std::cout << "  Removed " << removed
			<< " cached file(s). "
			"They will be regenerated on next use.\n";
	}

	return 0;
}

bool PromptRebuildTheoryCurves(
	const brill::SquareDetectorConfig *detectors[],
	int layer_count,
	const brill::AppConfig &config
) {
	std::cout << "\nCurrent detector thicknesses:\n";
	for (int i = 0; i < layer_count; ++i) {
		if (detectors[i]) {
			std::cout << "  " << brill::kT0LayerNames[i]
				<< ": " << detectors[i]->thickness_um
				<< " um\n";
		}
	}
	std::cout
		<< "Recalculate delta-E theory curves? [y/N]: ";
	std::cout.flush();

	std::string answer;
	std::getline(std::cin, answer);

	if (answer == "y" || answer == "Y") {
		RebuildTheoryCurves(config);
		return true;
	}
	return false;
}

int FillGcaliForLayers(
	const std::vector<brill::LoadedCut> &cuts,
	TGraph *graphs[4],
	const std::set<int> &allowed_layers,
	TGraph &gcali,
	std::set<std::pair<int, int>> &projectile_set
) {
	for (const auto &cut : cuts) {
		int layer = cut.pair->layer;
		if (allowed_layers.find(layer)
			== allowed_layers.end()) continue;

		const brill::T0ParticlePidInfo *info =
			brill::GetPidInfo(
				layer,
				cut.particle.charge,
				cut.particle.mass);
		if (!info) {
			std::cerr << "Error: No pid_info for (layer="
				<< layer
				<< ", Z=" << cut.particle.charge
				<< ", A=" << cut.particle.mass
				<< ").\n";
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
				gcali.AddPoint(
					shallow + info->offset,
					deep);
				++hit_count;
			}
		}
		std::cout << "  " << cut.pair->name
			<< " Z=" << cut.particle.charge
			<< " A=" << cut.particle.mass
			<< ": " << hit_count
			<< " points selected\n";
		projectile_set.insert(
			{cut.particle.charge,
			cut.particle.mass});
	}
	return 0;
}

void CreateCalibratedTH2(
	TH2F *h2_raw,
	const char *name,
	const char *title,
	double p0_deep,
	double p1_deep,
	double p0_shallow,
	double p1_shallow,
	TFile &opf
) {
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
		name,
		title,
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
}

void DrawTheoryCurves(
	TH2F *h2_cal,
	const brill::SquareDetectorConfig *detectors[],
	int shallow_layer,
	int deep_layer,
	const brill::AppConfig &config,
	TFile &opf,
	const char *canvas_name,
	const char *canvas_title
) {
	TCanvas *c = new TCanvas(
		canvas_name, canvas_title, 1200, 900);
	c->cd();
	h2_cal->Draw("colz");

	double t_shallow = detectors[shallow_layer]->thickness_um;
	double t_deep = detectors[deep_layer]->thickness_um;

	std::vector<TGraph*> theory_curves;
	for (const auto &p : brill::kCommonParticles) {
		std::string cache_path = TString::Format(
			"%s/si_z%d_a%d.root",
			brill::JoinPath(
				config.workspace,
				config.paths.energy_calculator).c_str(),
			p.charge,
			p.mass).Data();

		try {
			brill::RangeEnergyCalculator rec(
				p.charge, p.mass,
				brill::SiliconMaterial(),
				cache_path);

			TGraph *curve = brill::GenerateTheoryCurve(
				rec,
				t_shallow,
				t_deep);
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

} // namespace

int main(int argc, char **argv) {
	cxxopts::Options options(
		"calibrate_t0_v1",
		"T0 calibration v1: two-stage fit "
		"(stage1: d2-d3/d3-d4/d4-s1, "
		"stage2: d1-d2 with fixed d2).");
	options.add_options()
		("h,help", "Print help information.")
		("r,run", "Start run number.",
			cxxopts::value<int>(), "run")
		("t,trigger", "Trigger type.",
			cxxopts::value<std::string>(), "trigger")
		("c,config", "Config file path.",
			cxxopts::value<std::string>()
				->default_value("config.toml"),
			"file");

	auto result = options.parse(argc, argv);
	if (result.count("help")) {
		PrintUsage(options);
		return 0;
	}
	if (!result.count("run")) {
		std::cerr
			<< "Error: Missing required option --run.\n";
		PrintUsage(options);
		return 1;
	}

	brill::AppConfig config;
	if (brill::LoadConfig(
		result["config"].as<std::string>(),
		config)) {
		return 1;
	}
	if (result.count("trigger")) {
		config.trigger =
			result["trigger"].as<std::string>();
	}

	const int run = result["run"].as<int>();

	const brill::SquareDetectorConfig *detectors[
		brill::kT0LayerCount] = {
		brill::FindDetectorConfig(config, "t0d1"),
		brill::FindDetectorConfig(config, "t0d2"),
		brill::FindDetectorConfig(config, "t0d3"),
		brill::FindDetectorConfig(config, "t0d4"),
		brill::FindDetectorConfig(config, "t0s")
	};
	for (int i = 0; i < brill::kT0LayerCount; ++i) {
		if (!detectors[i]) {
			std::cerr
				<< "Error: Missing detector config for "
				<< brill::kT0LayerNames[i] << ".\n";
			return 1;
		}
	}

	PromptRebuildTheoryCurves(
		detectors, brill::kT0LayerCount, config);

	const std::string trigger_infix =
		brill::TriggerInfix(config.trigger);

	const std::string estimate_dir = brill::JoinPath(
		config.workspace, config.paths.estimate);
	std::string precal_path = brill::FindPreCalibrationFile(
		estimate_dir, trigger_infix, run);
	if (precal_path.empty()) {
		std::cerr << "Error: No pre_calibration file "
			"found for run "
			<< run << " in " << estimate_dir << ".\n";
		return 1;
	}
	std::cout << "Reading pre_calibration: "
		<< precal_path << "\n";

	TFile precal_file(precal_path.c_str(), "read");
	if (precal_file.IsZombie()) {
		std::cerr << "Error: Open pre_calibration file "
			<< precal_path << " failed.\n";
		return 1;
	}

	TGraph *graphs[4] = {nullptr};
	const char *graph_names[4] = {
		"g_d1d2", "g_d2d3", "g_d3d4", "g_d4s"};
	for (int i = 0; i < 4; ++i) {
		graphs[i] = dynamic_cast<TGraph*>(
			precal_file.Get(graph_names[i]));
		if (!graphs[i]) {
			std::cerr << "Error: TGraph "
				<< graph_names[i]
				<< " not found in "
				<< precal_path << ".\n";
			return 1;
		}
		std::cout << "  " << graph_names[i] << ": "
			<< graphs[i]->GetN() << " points\n";
	}

	const std::string cut_dir = "src/brill/Cut";
	if (!std::filesystem::exists(cut_dir)) {
		std::cerr << "Error: Cut directory "
			<< cut_dir << " not found.\n";
		return 1;
	}
	std::cout << "\nLoading cuts from: "
		<< cut_dir << "\n";

	std::vector<brill::LoadedCut> cuts;
	if (brill::LoadCuts(cut_dir, "_stop", cuts)) {
		std::cerr << "Error: No cuts loaded from "
			<< cut_dir << ".\n";
		return 1;
	}

	std::string calibration_dir = brill::JoinPath(
		config.workspace,
		config.paths.calibration);
	std::filesystem::create_directories(
		calibration_dir);

	TString output_path = TString::Format(
		"%s/t0_%s%04d.root",
		calibration_dir.c_str(),
		trigger_infix.c_str(),
		run
	);
	TFile opf(output_path, "recreate");
	std::cout << "\nOutput ROOT file: "
		<< output_path << "\n";

	std::cout
		<< "\n========================================\n"
		<< "Stage 1: Fit d2-d3, d3-d4, d4-s1 "
		<< "(8 parameters: d2,d3,d4,s1)\n"
		<< "========================================\n";

	std::set<int> stage1_layers = {1, 2, 3};
	std::set<std::pair<int, int>> stage1_projectile_set;
	TGraph gcali_stage1;
	FillGcaliForLayers(
		cuts, graphs, stage1_layers,
		gcali_stage1, stage1_projectile_set);

	if (gcali_stage1.GetN() == 0) {
		std::cerr << "Error: No calibration points "
			"for stage 1.\n";
		return 1;
	}
	std::cout << "Stage 1 calibration points: "
		<< gcali_stage1.GetN() << "\n";

	std::vector<std::pair<int, int>> stage1_projectiles(
		stage1_projectile_set.begin(),
		stage1_projectile_set.end());

	double stage1_init_params[8] = {
		0.0, 0.006,
		0.0, 0.006,
		0.0, 0.003,
		0.0, 0.003
	};

	PidFitFuncStage1 pid_fit_stage1(
		stage1_projectiles, config);
	TF1 fcali_stage1(
		"fcali_stage1", pid_fit_stage1,
		25000.0, 220000.0, 8);
	fcali_stage1.SetNpx(10000);
	for (int i = 0; i < 8; ++i) {
		fcali_stage1.SetParameter(
			i, stage1_init_params[i]);
	}
	for (int i = 0; i < 8; ++i) {
		fcali_stage1.SetParLimits(
			i,
			(i % 2 == 0) ? 0.0 : 0.0,
			(i % 2 == 0) ? 100.0 : 1.0);
	}

	std::cout << "Fitting stage 1 with \"R S\"...\n";
	gcali_stage1.Fit(&fcali_stage1, "R S");

	std::cout << "Stage 1 fit result:\n"
		<< "  chi2     = "
		<< fcali_stage1.GetChisquare() << "\n"
		<< "  ndf      = "
		<< fcali_stage1.GetNDF() << "\n"
		<< "  chi2/ndf = "
		<< fcali_stage1.GetChisquare()
			/ fcali_stage1.GetNDF() << "\n"
		<< "  prob     = "
		<< fcali_stage1.GetProb() << "\n";

	double *stage1_pars = fcali_stage1.GetParameters();
	std::cout << "Stage 1 parameters:\n";
	const char *stage1_layer_names[4] = {
		"t0d2", "t0d3", "t0d4", "t0s"};
	for (int i = 0; i < 4; ++i) {
		std::cout << "  " << stage1_layer_names[i]
			<< ": p0 = " << stage1_pars[i * 2]
			<< ", p1 = " << stage1_pars[i * 2 + 1]
			<< "\n";
	}

	double d2_p0 = stage1_pars[0];
	double d2_p1 = stage1_pars[1];

	std::cout
		<< "\n========================================\n"
		<< "Stage 2: Fit d1-d2 "
		<< "(2 parameters: d1 only, d2 fixed)\n"
		<< "  d2 p0 = " << d2_p0
		<< ", d2 p1 = " << d2_p1 << "\n"
		<< "========================================\n";

	std::set<int> stage2_layers = {0};
	std::set<std::pair<int, int>> stage2_projectile_set;
	TGraph gcali_stage2;
	FillGcaliForLayers(
		cuts, graphs, stage2_layers,
		gcali_stage2, stage2_projectile_set);

	if (gcali_stage2.GetN() == 0) {
		std::cerr << "Error: No calibration points "
			"for stage 2.\n";
		return 1;
	}
	std::cout << "Stage 2 calibration points: "
		<< gcali_stage2.GetN() << "\n";

	std::vector<std::pair<int, int>> stage2_projectiles(
		stage2_projectile_set.begin(),
		stage2_projectile_set.end());

	double stage2_init_params[2] = {0.0, 0.002};

	PidFitFuncStage2 pid_fit_stage2(
		stage2_projectiles, config,
		d2_p0, d2_p1);
	TF1 fcali_stage2(
		"fcali_stage2", pid_fit_stage2,
		0.0, 25000.0, 2);
	fcali_stage2.SetNpx(10000);
	for (int i = 0; i < 2; ++i) {
		fcali_stage2.SetParameter(
			i, stage2_init_params[i]);
	}
	for (int i = 0; i < 2; ++i) {
		fcali_stage2.SetParLimits(
			i,
			(i % 2 == 0) ? 0.0 : 0.0,
			(i % 2 == 0) ? 100.0 : 1.0);
	}

	std::cout << "Fitting stage 2 with \"R S\"...\n";
	gcali_stage2.Fit(&fcali_stage2, "R S");

	std::cout << "Stage 2 fit result:\n"
		<< "  chi2     = "
		<< fcali_stage2.GetChisquare() << "\n"
		<< "  ndf      = "
		<< fcali_stage2.GetNDF() << "\n"
		<< "  chi2/ndf = "
		<< fcali_stage2.GetChisquare()
			/ fcali_stage2.GetNDF() << "\n"
		<< "  prob     = "
		<< fcali_stage2.GetProb() << "\n";

	double *stage2_pars = fcali_stage2.GetParameters();
	std::cout << "Stage 2 parameters:\n"
		<< "  t0d1: p0 = " << stage2_pars[0]
		<< ", p1 = " << stage2_pars[1] << "\n";

	double final_parameters[10] = {
		stage2_pars[0], stage2_pars[1],
		d2_p0,          d2_p1,
		stage1_pars[2], stage1_pars[3],
		stage1_pars[4], stage1_pars[5],
		stage1_pars[6], stage1_pars[7]
	};

	std::cout
		<< "\n========================================\n"
		<< "Final combined parameters:\n"
		<< "========================================\n";
	for (int layer = 0;
		layer < brill::kT0LayerCount;
		++layer) {
		std::cout << "  "
			<< brill::kT0LayerNames[layer]
			<< ": p0 = "
			<< final_parameters[layer * 2]
			<< ", p1 = "
			<< final_parameters[layer * 2 + 1]
			<< "\n";
	}

	TString txt_path = TString::Format(
		"%s/t0_%04d.txt",
		calibration_dir.c_str(),
		run
	);
	brill::WriteCalibrationParameters(
		txt_path.Data(),
		final_parameters,
		brill::kT0LayerCount);
	std::cout << "\nCalibration written to: "
		<< txt_path << "\n";

	const char *th2_names[4] = {
		"d1d2", "d2d3", "d3d4", "d4s"};
	const char *th2_cal_names[4] = {
		"d1d2_cal", "d2d3_cal", "d3d4_cal",
		"d4s_cal"};

	std::cout << "\nCalibrating TH2Fs and "
		"drawing theory curves...\n";

	for (int i = 0; i < 4; ++i) {
		TH2F *h2_raw = dynamic_cast<TH2F*>(
			precal_file.Get(th2_names[i]));
		if (!h2_raw) {
			std::cerr << "Warning: TH2F "
				<< th2_names[i]
				<< " not found, skipping.\n";
			continue;
		}

		int deep_layer = i + 1;
		int shallow_layer = i;

		double p0_deep =
			final_parameters[deep_layer * 2];
		double p1_deep =
			final_parameters[deep_layer * 2 + 1];
		double p0_shallow =
			final_parameters[shallow_layer * 2];
		double p1_shallow =
			final_parameters[shallow_layer * 2 + 1];

		TString title = TString::Format(
			"%s Calibrated;%s [MeV];%s [MeV]",
			brill::kT0LayerPairs[i].name,
			brill::kT0LayerNames[deep_layer],
			brill::kT0LayerNames[shallow_layer]);

		CreateCalibratedTH2(
			h2_raw,
			th2_cal_names[i],
			title.Data(),
			p0_deep, p1_deep,
			p0_shallow, p1_shallow,
			opf);

		TH2F *h2_cal = dynamic_cast<TH2F*>(
			opf.Get(th2_cal_names[i]));

		if (h2_cal) {
			DrawTheoryCurves(
				h2_cal,
				detectors,
				shallow_layer,
				deep_layer,
				config,
				opf,
				TString::Format(
					"c_%s_v1",
					th2_names[i]),
				TString::Format(
					"%s Calibrated with "
					"Theory (v1)",
					brill::kT0LayerPairs[i]
						.name));
		}
	}

	opf.cd();
	gcali_stage1.Write("gcali_stage1");
	fcali_stage1.Write("fcali_stage1");
	gcali_stage2.Write("gcali_stage2");
	fcali_stage2.Write("fcali_stage2");
	opf.Close();

	std::cout << "\ncalibrate_t0_v1 done.\n";
	return 0;
}