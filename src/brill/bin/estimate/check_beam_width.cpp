#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TCutG.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TString.h>
#include <TTree.h>

#include "external/cxxopts.hpp"
#include "include/config.h"
#include "include/event/t0/dssd_match_event.h"
#include "include/utils.h"

namespace {

constexpr int kCalibrationLayers = 5;
constexpr int kMaxGroups = 10;

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

template<typename Func>
Long64_t ProcessEvents(
	TTree *d2_tree,
	TTree *d3_tree,
	brill::DssdMatchEvent &d2_event,
	brill::DssdMatchEvent &d3_event,
	const CalibrationParameters &calibration,
	TCutG *beam_cut,
	Long64_t max_passed,
	Func on_pass
) {
	Long64_t total = d2_tree->GetEntries();
	Long64_t passed = 0;
	Long64_t last_percentage = -1;

	for (Long64_t entry = 0; entry < total; ++entry) {
		Long64_t percentage = total > 0 ? entry * 100ll / total : 100ll;
		if (percentage > last_percentage) {
			last_percentage = percentage;
			std::printf("\b\b\b\b%3lld%%", percentage);
			std::fflush(stdout);
		}

		d2_tree->GetEntry(entry);
		d3_tree->GetEntry(entry);

		if (d2_event.num != 1 || d3_event.num != 1) continue;

		double dx = d3_event.x[0] - d2_event.x[0];
		double dy = d3_event.y[0] - d2_event.y[0];
		if (dx * dx + dy * dy >= 4.0) continue;

		double d2_cal_e = CalibrateEnergy(calibration, 1, d2_event.energy[0]);
		double d3_cal_e = CalibrateEnergy(calibration, 2, d3_event.energy[0]);

		if (!beam_cut->IsInside(d3_cal_e, d2_cal_e)) continue;

		double E = d2_cal_e + d3_cal_e;
		on_pass(E, passed);
		++passed;

		if (max_passed > 0 && passed >= max_passed) break;
	}
	std::printf("\b\b\b\b100%%\n");
	return passed;
}

int OpenRunFiles(
	const std::string &match_dir,
	const std::string &trigger_infix,
	int run,
	TFile *&d2_file,
	TFile *&d3_file,
	TTree *&d2_tree,
	TTree *&d3_tree,
	brill::DssdMatchEvent &d2_event,
	brill::DssdMatchEvent &d3_event
) {
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

	if (!std::filesystem::exists(d2_path.Data())) {
		std::cerr << "Error: d2 match file not found: " << d2_path << "\n";
		return 1;
	}
	if (!std::filesystem::exists(d3_path.Data())) {
		std::cerr << "Error: d3 match file not found: " << d3_path << "\n";
		return 1;
	}

	d2_file = new TFile(d2_path, "read");
	d3_file = new TFile(d3_path, "read");

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

struct GaussianFitResult {
	double peak = 0.0;
	double fwhm = 0.0;
	bool valid = false;
};

void DoGaussianFits(
	const std::vector<TH1D*> &histograms,
	const std::vector<double> &x_values,
	std::vector<GaussianFitResult> &results,
	TFile &output_file
) {
	results.clear();
	for (size_t i = 0; i < histograms.size(); ++i) {
		TH1D *h = histograms[i];
		GaussianFitResult r;
		if (!h || h->GetEntries() == 0) {
			results.push_back(r);
			continue;
		}

		int fit_status = h->Fit("gaus", "Q", "", 475.0, 485.0);
		if (fit_status == 0) {
			TF1 *fit_func = h->GetFunction("gaus");
			if (fit_func) {
				r.peak = fit_func->GetParameter(1);
				double sigma = std::abs(fit_func->GetParameter(2));
				r.fwhm = 2.355 * sigma;
				r.valid = true;
			}
		}
		results.push_back(r);
	}

	std::vector<double> peak_x, peak_y;
	std::vector<double> fwhm_x, fwhm_y;
	for (size_t i = 0; i < results.size(); ++i) {
		if (!results[i].valid) continue;
		peak_x.push_back(x_values[i]);
		peak_y.push_back(results[i].peak);
		fwhm_x.push_back(x_values[i]);
		fwhm_y.push_back(results[i].fwhm);
	}

	if (!peak_x.empty()) {
		TGraph g_peak(static_cast<int>(peak_x.size()), peak_x.data(), peak_y.data());
		g_peak.SetName("g_peak");
		g_peak.SetTitle("Gaussian Peak Position;Index;Peak (MeV)");
		g_peak.SetMarkerStyle(20);
		g_peak.SetMarkerSize(1.2);

		TCanvas c_peak("c_peak", "Gaussian Peak Position", 800, 600);
		g_peak.Draw("AP");
		output_file.cd();
		c_peak.Write();
		g_peak.Write();
	}

	if (!fwhm_x.empty()) {
		TGraph g_fwhm(static_cast<int>(fwhm_x.size()), fwhm_x.data(), fwhm_y.data());
		g_fwhm.SetName("g_fwhm");
		g_fwhm.SetTitle("Gaussian FWHM (2.355#sigma);Index;FWHM (MeV)");
		g_fwhm.SetMarkerStyle(20);
		g_fwhm.SetMarkerSize(1.2);

		TCanvas c_fwhm("c_fwhm", "Gaussian FWHM", 800, 600);
		g_fwhm.Draw("AP");
		output_file.cd();
		c_fwhm.Write();
		g_fwhm.Write();
	}
}

} // namespace

int main(int argc, char **argv) {
	cxxopts::Options options("check_beam_width", "Check beam width by energy spread.");
	options.add_options()
		("h,help", "Print help information.")
		(
			"c,config",
			"Config file path.",
			cxxopts::value<std::string>()->default_value("config.toml"),
			"file"
		)
		("r,run", "Run number.", cxxopts::value<int>(), "run")
		("e,end-run", "End run number for multi-run mode.", cxxopts::value<int>(), "end")
		("t,trigger", "Trigger type.", cxxopts::value<std::string>(), "trigger")
		(
			"n,group-power",
			"Group size as 10^n events per histogram (default: 6).",
			cxxopts::value<int>()->default_value("6"),
			"n"
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
	if (result.count("trigger")) {
		config.trigger = result["trigger"].as<std::string>();
	}

	const int group_power = result["group-power"].as<int>();
	const long long group_size = static_cast<long long>(std::pow(10.0, group_power));

	const bool multi_run = result.count("end-run") > 0;
	const int end_run = multi_run ? result["end-run"].as<int>() : run;

	const std::string trigger_infix = brill::TriggerInfix(config.trigger);
	const std::string match_dir = brill::JoinPath(config.workspace, config.paths.match);

	TCutG *beam_cut = LoadBeamCut();
	if (!beam_cut) return 1;

	const std::string output_dir = brill::JoinPath(config.workspace, config.paths.beam);
	std::filesystem::create_directories(output_dir);

	TString output_path;
	if (multi_run) {
		output_path = TString::Format(
			"%s/check_beam_width_%s%04d-%04d.root",
			output_dir.c_str(),
			trigger_infix.c_str(),
			run, end_run
		);
	} else {
		output_path = TString::Format(
			"%s/check_beam_width_%s%04d.root",
			output_dir.c_str(),
			trigger_infix.c_str(),
			run
		);
	}
	TFile output_file(output_path, "recreate");

	if (!multi_run) {
		if (brill::IsJumpRun(config, run)) {
			std::cout << "Skipping jump run " << run << ".\n";
			return 0;
		}

		TFile *d2_file = nullptr;
		TFile *d3_file = nullptr;
		TTree *d2_tree = nullptr;
		TTree *d3_tree = nullptr;
		brill::DssdMatchEvent d2_event;
		brill::DssdMatchEvent d3_event;

		if (OpenRunFiles(match_dir, trigger_infix, run,
			d2_file, d3_file, d2_tree, d3_tree, d2_event, d3_event)) {
			return 1;
		}

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
			"d2(p0=%+.6f,p1=%.6f) d3(p0=%+.6f,p1=%.6f)\n",
			calib_run,
			calibration.p0[1], calibration.p1[1],
			calibration.p0[2], calibration.p1[2]
		);

		std::printf("Run %d: group size = 10^%d = %lld, max %d groups\n",
			run, group_power, group_size, kMaxGroups);

		TH1D *h_groups[kMaxGroups] = {nullptr};
		for (int g = 0; g < kMaxGroups; ++g) {
			h_groups[g] = new TH1D(
				TString::Format("h_beam_E_group%d", g),
				TString::Format("Beam E group %d (10^{%d} evts);E (MeV);Counts", g, group_power),
				1000, 0.0, 1000.0
			);
			h_groups[g]->SetDirectory(0);
			h_groups[g]->SetLineColor(g + 1);
		}

		std::printf("Checking beam width   0%%");
		std::fflush(stdout);
		Long64_t passed_count = ProcessEvents(
			d2_tree, d3_tree, d2_event, d3_event,
			calibration, beam_cut,
			kMaxGroups * group_size,
			[&](double E, Long64_t passed_index) {
				int group = static_cast<int>(passed_index / group_size);
				if (group >= 0 && group < kMaxGroups) {
					h_groups[group]->Fill(E);
				}
			}
		);
		std::printf("Passed events: %lld (in %d groups of 10^%d)\n",
			passed_count,
			std::min(kMaxGroups, static_cast<int>(passed_count / group_size) + 1),
			group_power);

		d2_file->Close();
		d3_file->Close();
		delete d2_file;
		delete d3_file;

		TCanvas canvas("c_beam_width", "Beam Energy Width Check", 1400, 800);
		TLegend legend(0.75, 0.65, 0.90, 0.90);
		bool has_drawn = false;
		for (int g = 0; g < kMaxGroups; ++g) {
			if (h_groups[g]->GetEntries() == 0) continue;
			h_groups[g]->SetLineWidth(2);
			if (!has_drawn) {
				h_groups[g]->Draw("HIST");
				has_drawn = true;
			} else {
				h_groups[g]->Draw("HIST SAME");
			}
			legend.AddEntry(
				h_groups[g],
				TString::Format("Group %d (%lld evts)", g,
					static_cast<Long64_t>(h_groups[g]->GetEntries())),
				"l"
			);
		}
		legend.Draw();

		output_file.cd();
		canvas.Write();
		for (int g = 0; g < kMaxGroups; ++g) {
			h_groups[g]->Write();
		}

		std::vector<TH1D*> h_vec;
		std::vector<double> x_vals;
		for (int g = 0; g < kMaxGroups; ++g) {
			if (h_groups[g]->GetEntries() == 0) continue;
			h_vec.push_back(h_groups[g]);
			x_vals.push_back(static_cast<double>(g));
		}
		std::vector<GaussianFitResult> fit_results;
		DoGaussianFits(h_vec, x_vals, fit_results, output_file);
		std::printf("Gaussian fit (475-485 MeV):\n");
		for (size_t i = 0; i < fit_results.size(); ++i) {
			if (fit_results[i].valid) {
				std::printf("  Group %d: peak=%.4f MeV, FWHM=%.4f MeV\n",
					static_cast<int>(x_vals[i]),
					fit_results[i].peak,
					fit_results[i].fwhm);
			} else {
				std::printf("  Group %d: fit failed\n", static_cast<int>(x_vals[i]));
			}
		}

	} else {
		std::printf("Multi-run mode: run %d to %d, max 10^%d = %lld events per run\n",
			run, end_run, group_power, group_size);

		std::vector<TH1D*> h_runs;
		std::vector<int> processed_runs;
		TCanvas canvas("c_beam_width_multi", "Beam Energy Width Check (Multi-run)", 1400, 800);
		TLegend legend(0.75, 0.65, 0.90, 0.90);

		bool first_drawn = false;
		for (int r = run; r <= end_run; ++r) {
			if (brill::IsJumpRun(config, r)) {
				std::cout << "Skipping jump run " << r << ".\n";
				continue;
			}

			TFile *d2_file = nullptr;
			TFile *d3_file = nullptr;
			TTree *d2_tree = nullptr;
			TTree *d3_tree = nullptr;
			brill::DssdMatchEvent d2_event;
			brill::DssdMatchEvent d3_event;

			if (OpenRunFiles(match_dir, trigger_infix, r,
				d2_file, d3_file, d2_tree, d3_tree, d2_event, d3_event)) {
				continue;
			}

			int calib_run = brill::GetT0CalibrationRun(config, r);
			const std::string calibration_path = TString::Format(
				"%s/t0_%04d.txt",
				brill::JoinPath(config.workspace, config.paths.calibration).c_str(),
				calib_run
			).Data();

			CalibrationParameters calibration;
			if (ReadCalibrationParameters(calibration_path, calibration)) {
				d2_file->Close();
				d3_file->Close();
				delete d2_file;
				delete d3_file;
				continue;
			}
			std::printf(
				"Calibration (t0_%04d.txt): "
				"d2(p0=%+.6f,p1=%.6f) d3(p0=%+.6f,p1=%.6f)\n",
				calib_run,
				calibration.p0[1], calibration.p1[1],
				calibration.p0[2], calibration.p1[2]
			);

			TH1D *h = new TH1D(
				TString::Format("h_beam_E_run%04d", r),
				TString::Format("Beam E run %04d (10^{%d} evts);E (MeV);Counts", r, group_power),
				1000, 0.0, 1000.0
			);
			h->SetDirectory(0);
			int color = 1 + ((r - run) % 9);
			h->SetLineColor(color);

			std::printf("Run %04d   0%%", r);
			std::fflush(stdout);
			Long64_t passed_count = ProcessEvents(
				d2_tree, d3_tree, d2_event, d3_event,
				calibration, beam_cut,
				group_size,
				[&](double E, Long64_t passed_index) {
					(void)passed_index;
					h->Fill(E);
				}
			);
			std::printf("Run %04d: passed %lld events\n", r, passed_count);

			d2_file->Close();
			d3_file->Close();
			delete d2_file;
			delete d3_file;

			output_file.cd();
			h->Write();

			h->SetLineWidth(2);
			if (!first_drawn) {
				h->Draw("HIST");
				first_drawn = true;
			} else {
				h->Draw("HIST SAME");
			}
			legend.AddEntry(
				h,
				TString::Format("Run %04d (%lld evts)", r, passed_count),
				"l"
			);

			h_runs.push_back(h);
			processed_runs.push_back(r);
		}

		output_file.cd();
		legend.Draw();
		canvas.Write();

		std::vector<double> run_x_vals;
		for (size_t i = 0; i < processed_runs.size(); ++i) {
			run_x_vals.push_back(static_cast<double>(processed_runs[i]));
		}
		std::vector<GaussianFitResult> fit_results;
		DoGaussianFits(h_runs, run_x_vals, fit_results, output_file);
		std::printf("Gaussian fit (475-485 MeV):\n");
		for (size_t i = 0; i < fit_results.size(); ++i) {
			if (fit_results[i].valid) {
				std::printf("  Run %d: peak=%.4f MeV, FWHM=%.4f MeV\n",
					static_cast<int>(run_x_vals[i]),
					fit_results[i].peak,
					fit_results[i].fwhm);
			} else {
				std::printf("  Run %d: fit failed\n", static_cast<int>(run_x_vals[i]));
			}
		}
	}

	output_file.Close();
	std::printf("Output written to %s\n", output_path.Data());
	return 0;
}