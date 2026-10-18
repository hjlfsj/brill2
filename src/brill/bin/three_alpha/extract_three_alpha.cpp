#include "external/cxxopts.hpp"
#include "include/config.h"
#include "include/three_alpha/extract.h"
#include "include/three_alpha/three_alpha_event.h"
#include "include/event/beam/beam_sort.h"
#include "include/event/ppac/ppac_track_event.h"
#include "include/event/t0/dssd_match_event.h"
#include "include/physics/kinematics.h"
#include "include/utils.h"

#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TCutG.h>
#include <TROOT.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char **argv) {
	cxxopts::Options options("extract_three_alpha", "Extract three-alpha events from T0 detector");
	options.add_options()
		("h,help", "Print help information.")
		("r,run", "Start run number.", cxxopts::value<int>())
		("e,end-run", "End run number.", cxxopts::value<int>())
		("t,trigger", "Trigger type.", cxxopts::value<std::string>())
		("c,config", "Config file path.", cxxopts::value<std::string>()->default_value("config.toml"))
	;

	auto result = options.parse(argc, argv);

	if (result.count("help")) {
		std::cout << options.help() << std::endl;
		return 0;
	}

	if (!result.count("run") || !result.count("trigger")) {
		std::cerr << "Error: --run and --trigger are required.\n";
		std::cout << options.help() << std::endl;
		return 1;
	}

	int run_start = result["run"].as<int>();
	int run_end = result.count("end-run") ? result["end-run"].as<int>() : run_start;
	std::string trigger = result["trigger"].as<std::string>();
	std::string config_path = result["config"].as<std::string>();

	brill::AppConfig config;
	if (brill::LoadConfig(config_path, config)) {
		std::cerr << "Error: Load config failed.\n";
		return 1;
	}

	std::string cut_dir = "/home/ribll2026/ribll2026_www/github_code/brill2/src/brill/Cut";

	TString cut_d12_path = TString::Format("%s/cal_d1_d2_tran_He4_cut.C", cut_dir.c_str());
	gROOT->ProcessLine(TString::Format(".x %s", cut_d12_path.Data()));
	TCutG *d12_cut = (TCutG*)gROOT->FindObject("cal_d1_d2_tran_He4_cut");
	if (!d12_cut) {
		std::cerr << "Error: Load TCutG from " << cut_d12_path << " failed.\n";
		return 1;
	}
	printf("Cut loaded: %s\n", cut_d12_path.Data());

	TString cut_d23_path = TString::Format("%s/cal_d2_d3_tran_He4_cut.C", cut_dir.c_str());
	gROOT->ProcessLine(TString::Format(".x %s", cut_d23_path.Data()));
	TCutG *d23_cut = (TCutG*)gROOT->FindObject("cal_d2_d3_tran_He4_cut");
	if (!d23_cut) {
		std::cerr << "Error: Load TCutG from " << cut_d23_path << " failed.\n";
		return 1;
	}
	printf("Cut loaded: %s\n", cut_d23_path.Data());

	TString cut_d34_stop_path = TString::Format("%s/cal_d3_d4_stop_He4_cut.C", cut_dir.c_str());
	gROOT->ProcessLine(TString::Format(".x %s", cut_d34_stop_path.Data()));
	TCutG *d34_stop_cut = (TCutG*)gROOT->FindObject("cal_d3_d4_stop_He4_cut");
	if (!d34_stop_cut) {
		std::cerr << "Error: Load TCutG from " << cut_d34_stop_path << " failed.\n";
		return 1;
	}
	printf("Cut loaded: %s\n", cut_d34_stop_path.Data());

	TString cut_d34_tran_path = TString::Format("%s/cal_d3_d4_tran_He4_cut.C", cut_dir.c_str());
	gROOT->ProcessLine(TString::Format(".x %s", cut_d34_tran_path.Data()));
	TCutG *d34_tran_cut = (TCutG*)gROOT->FindObject("cal_d3_d4_tran_He4_cut");
	if (!d34_tran_cut) {
		std::cerr << "Error: Load TCutG from " << cut_d34_tran_path << " failed.\n";
		return 1;
	}
	printf("Cut loaded: %s\n", cut_d34_tran_path.Data());

	const brill::SquareDetectorConfig *d1_cfg = brill::FindDetectorConfig(config, "t0d1");
	const brill::SquareDetectorConfig *d2_cfg = brill::FindDetectorConfig(config, "t0d2");
	const brill::SquareDetectorConfig *d3_cfg = brill::FindDetectorConfig(config, "t0d3");
	const brill::SquareDetectorConfig *d4_cfg = brill::FindDetectorConfig(config, "t0d4");

	if (!d1_cfg || !d2_cfg || !d3_cfg || !d4_cfg) {
		std::cerr << "Error: t0d1/t0d2/t0d3/t0d4 config not found.\n";
		return 1;
	}

	double total_thickness_um = d1_cfg->thickness_um + d2_cfg->thickness_um
	                          + d3_cfg->thickness_um + d4_cfg->thickness_um;

	printf("Si total thickness: %.1f um (d1=%.1f, d2=%.1f, d3=%.1f, d4=%.1f)\n",
		total_thickness_um,
		d1_cfg->thickness_um, d2_cfg->thickness_um,
		d3_cfg->thickness_um, d4_cfg->thickness_um);

	std::string cache_dir = brill::JoinPath(config.workspace, config.paths.energy_calculator);
	std::string cache_path = TString::Format("%s/si_%.0fum_z2_a4.root",
		cache_dir.c_str(), total_thickness_um).Data();
	double e_test = brill::FindIncidentEnergy(42.122, total_thickness_um, cache_path);
	printf("=== Validation: alpha with 42.122 MeV deposited in %d-layer Si "
		"(total %.1f um) has incident energy = %.2f MeV ===\n",
		4, total_thickness_um, e_test);

	if (e_test <= 50.0 || std::isnan(e_test)) {
		std::cerr << "Error: Invalid punch-through correction. Check config.\n";
		return 1;
	}

	std::string match_dir = brill::JoinPath(config.workspace, config.paths.match);
	std::string track_dir = brill::JoinPath(config.workspace, config.paths.track);
	std::string beam_dir = brill::JoinPath(config.workspace, config.paths.beam);
	std::string output_dir = brill::JoinPath(config.workspace, config.paths.three_alpha);

	std::string trigger_infix = brill::TriggerInfix(trigger);

	TString output_path = TString::Format("%s/extract_three_alpha_%s%04d_%04d.root",
		output_dir.c_str(), trigger_infix.c_str(), run_start, run_end);

	TFile *output_file = new TFile(output_path, "recreate");
	if (!output_file || output_file->IsZombie()) {
		std::cerr << "Error: Create output file failed: " << output_path << "\n";
		return 1;
	}

	TTree *output_tree = new TTree("tree", "Three-alpha extracted events");
	brill::ThreeAlphaEvent out_event;
	brill::SetupOutputThreeAlpha(output_tree, out_event);

	TH2D *h_e1_He1_e2_He1 = new TH2D("h_e1_He1_e2_He1", "E1_He1 vs E2_He1;E2_He1 (MeV);E1_He1 (MeV)", 1000, 0, 400, 1000, 0, 25);
	TH2D *h_e2_He1_e3_He1 = new TH2D("h_e2_He1_e3_He1", "E2_He1 vs E3_He1;E3_He1 (MeV);E2_He1 (MeV)", 1000, 0, 350, 1000, 0, 400);
	TH2D *h_e3_He1_e4_He1 = new TH2D("h_e3_He1_e4_He1", "E3_He1 vs E4_He1;E4_He1 (MeV);E3_He1 (MeV)", 1000, 0, 250, 1000, 0, 300);

	TH2D *h_e1_He2_e2_He2 = new TH2D("h_e1_He2_e2_He2", "E1_He2 vs E2_He2;E2_He2 (MeV);E1_He2 (MeV)", 1000, 0, 400, 1000, 0, 25);
	TH2D *h_e2_He2_e3_He2 = new TH2D("h_e2_He2_e3_He2", "E2_He2 vs E3_He2;E3_He2 (MeV);E2_He2 (MeV)", 1000, 0, 350, 1000, 0, 400);
	TH2D *h_e3_He2_e4_He2 = new TH2D("h_e3_He2_e4_He2", "E3_He2 vs E4_He2;E4_He2 (MeV);E3_He2 (MeV)", 1000, 0, 250, 1000, 0, 300);

	TH2D *h_e1_He3_e2_He3 = new TH2D("h_e1_He3_e2_He3", "E1_He3 vs E2_He3;E2_He3 (MeV);E1_He3 (MeV)", 1000, 0, 400, 1000, 0, 25);
	TH2D *h_e2_He3_e3_He3 = new TH2D("h_e2_He3_e3_He3", "E2_He3 vs E3_He3;E3_He3 (MeV);E2_He3 (MeV)", 1000, 0, 350, 1000, 0, 400);
	TH2D *h_e3_He3_e4_He3 = new TH2D("h_e3_He3_e4_He3", "E3_He3 vs E4_He3;E4_He3 (MeV);E3_He3 (MeV)", 1000, 0, 250, 1000, 0, 300);

	for (int run = run_start; run <= run_end; ++run) {
		if (brill::IsJumpRun(config, run)) {
			std::cout << "Skipping run " << run << " (jump run)\n";
			continue;
		}

		TString d1_path = TString::Format("%s/t0d1_%s%04d.root",
			match_dir.c_str(), trigger_infix.c_str(), run);
		TString d2_path = TString::Format("%s/t0d2_%s%04d.root",
			match_dir.c_str(), trigger_infix.c_str(), run);
		TString d3_path = TString::Format("%s/t0d3_%s%04d.root",
			match_dir.c_str(), trigger_infix.c_str(), run);
		TString d4_path = TString::Format("%s/t0d4_%s%04d.root",
			match_dir.c_str(), trigger_infix.c_str(), run);
		TString beam_sort_path = TString::Format("%s/beam_%s%04d.root",
			beam_dir.c_str(), trigger_infix.c_str(), run);
		TString ppac_track_path = TString::Format("%s/ppac_%s%04d.root",
			track_dir.c_str(), trigger_infix.c_str(), run);

		if (!std::filesystem::exists(d1_path.Data())) {
			std::cerr << "Warning: t0d1 match file not found: " << d1_path << ", skipping run " << run << "\n";
			continue;
		}

		TFile *d1_file = new TFile(d1_path, "read");
		TTree *d1_tree = (TTree*)d1_file->Get("tree");
		if (!d1_tree) {
			std::cerr << "Warning: t0d1 tree not found in " << d1_path << ", skipping run " << run << "\n";
			d1_file->Close();
			continue;
		}
		brill::DssdMatchEvent d1_event;
		brill::SetupInput(d1_tree, d1_event, "");

		TFile *d2_file = new TFile(d2_path, "read");
		TTree *d2_tree = (TTree*)d2_file->Get("tree");
		brill::DssdMatchEvent d2_event;
		if (d2_tree) brill::SetupInput(d2_tree, d2_event, "");

		TFile *d3_file = new TFile(d3_path, "read");
		TTree *d3_tree = (TTree*)d3_file->Get("tree");
		brill::DssdMatchEvent d3_event;
		if (d3_tree) brill::SetupInput(d3_tree, d3_event, "");

		TFile *d4_file = new TFile(d4_path, "read");
		TTree *d4_tree = (TTree*)d4_file->Get("tree");
		brill::DssdMatchEvent d4_event;
		if (d4_tree) brill::SetupInput(d4_tree, d4_event, "");

		TFile *beam_file = nullptr;
		TTree *beam_tree = nullptr;
		bool beam_is_14O = false;
		bool beam_is_13N = false;
		bool beam_is_12C = false;
		if (std::filesystem::exists(beam_sort_path.Data())) {
			beam_file = new TFile(beam_sort_path, "read");
			beam_tree = (TTree*)beam_file->Get("tree");
			if (beam_tree) {
				brill::SetupInputSortBeamTree(beam_tree, beam_is_14O, beam_is_13N, beam_is_12C);
			}
		}

		TFile *ppac_file = nullptr;
		TTree *ppac_tree = nullptr;
		brill::PpacTrackEvent ppac_event;
		if (std::filesystem::exists(ppac_track_path.Data())) {
			ppac_file = new TFile(ppac_track_path, "read");
			ppac_tree = (TTree*)ppac_file->Get("tree");
			if (ppac_tree) {
				brill::SetupInput(ppac_tree, ppac_event, "");
			}
		}

		int calib_run = brill::GetT0CalibrationRun(config, run);
		std::string calibration_path = TString::Format(
			"%s/t0_%04d.txt",
			brill::JoinPath(config.workspace, config.paths.calibration).c_str(),
			calib_run
		).Data();
		brill::ThreeAlphaCalibration calib;
		if (brill::ReadThreeAlphaCalibration(calibration_path, calib)) {
			std::cerr << "Error: Read calibration from " << calibration_path << " failed, skipping run " << run << ".\n";
			continue;
		}
		printf("  calib (t0_%04d.txt): "
			"d1(p0=%+.6f,p1=%.6f) d2(p0=%+.6f,p1=%.6f) "
			"d3(p0=%+.6f,p1=%.6f) d4(p0=%+.6f,p1=%.6f)\n",
			calib_run,
			calib.p0[0], calib.p1[0],
			calib.p0[1], calib.p1[1],
			calib.p0[2], calib.p1[2],
			calib.p0[3], calib.p1[3]);

		Long64_t n_entries = d1_tree->GetEntries();
		printf("Run %d: %lld events", run, n_entries);
		fflush(stdout);

		Long64_t last_percentage = 0;
		for (Long64_t i = 0; i < n_entries; ++i) {
			if (i * 100ll / n_entries > last_percentage) {
				last_percentage = i * 100ll / n_entries;
				printf("\b\b\b\b%3lld%%", last_percentage);
				fflush(stdout);
			}
			d1_tree->GetEntry(i);
			if (d2_tree) d2_tree->GetEntry(i);
			if (d3_tree) d3_tree->GetEntry(i);
			if (d4_tree) d4_tree->GetEntry(i);
			if (beam_tree) beam_tree->GetEntry(i);
			if (ppac_tree) ppac_tree->GetEntry(i);

			if (!ppac_event.valid) continue;
			if (!brill::PassThreeAlphaCut(d1_event, d2_event, d3_event, d4_event)) continue;

			brill::ThreeAlphaResult cls = brill::ClassifyThreeAlpha(
				d1_event, d2_event, d3_event, d4_event,
				calib, d12_cut, d23_cut, d34_stop_cut, d34_tran_cut);
			if (!cls.passed) continue;

			double e1[3], e2[3], e3[3], e4[3];
			double et[3];
			bool stop[3];

			const int *he_indices[3] = {cls.He1_idx, cls.He2_idx, cls.He3_idx};

			for (int p = 0; p < 3; ++p) {
				e1[p] = brill::CalibrateThreeAlphaEnergy(calib, 0, d1_event.energy[he_indices[p][0]]);
				e2[p] = brill::CalibrateThreeAlphaEnergy(calib, 1, d2_event.energy[he_indices[p][1]]);
				e3[p] = brill::CalibrateThreeAlphaEnergy(calib, 2, d3_event.energy[he_indices[p][2]]);
				e4[p] = brill::CalibrateThreeAlphaEnergy(calib, 3, d4_event.energy[he_indices[p][3]]);

				double e_sum = e1[p] + e2[p] + e3[p] + e4[p];
				bool in_stop = d34_stop_cut->IsInside(e4[p], e3[p]);

				if (in_stop) {
					et[p] = e_sum;
					stop[p] = true;
				} else {
					et[p] = brill::FindIncidentEnergy(e_sum, total_thickness_um, cache_path);
					stop[p] = false;
				}
			}

			std::vector<int> sort_idx = {0, 1, 2};
			std::sort(sort_idx.begin(), sort_idx.end(),
				[&](int a, int b) { return et[a] > et[b]; });

			const int *he1 = he_indices[sort_idx[0]];
			const int *he2 = he_indices[sort_idx[1]];
			const int *he3 = he_indices[sort_idx[2]];

			out_event.run_number = run;
			out_event.entry = i;
			out_event.e1_He1 = e1[sort_idx[0]];
			out_event.e2_He1 = e2[sort_idx[0]];
			out_event.e3_He1 = e3[sort_idx[0]];
			out_event.e4_He1 = e4[sort_idx[0]];
			out_event.e1_He2 = e1[sort_idx[1]];
			out_event.e2_He2 = e2[sort_idx[1]];
			out_event.e3_He2 = e3[sort_idx[1]];
			out_event.e4_He2 = e4[sort_idx[1]];
			out_event.e1_He3 = e1[sort_idx[2]];
			out_event.e2_He3 = e2[sort_idx[2]];
			out_event.e3_He3 = e3[sort_idx[2]];
			out_event.e4_He3 = e4[sort_idx[2]];

			out_event.et_He1 = et[sort_idx[0]];
			out_event.et_He2 = et[sort_idx[1]];
			out_event.et_He3 = et[sort_idx[2]];

			out_event.stop_He1 = stop[sort_idx[0]];
			out_event.stop_He2 = stop[sort_idx[1]];
			out_event.stop_He3 = stop[sort_idx[2]];

			out_event.is_14O = beam_is_14O;
			out_event.is_13N = beam_is_13N;
			out_event.is_12C = beam_is_12C;

			out_event.ppac_valid = (ppac_event.valid != 0);
			out_event.target_x = ppac_event.target_x;
			out_event.target_y = ppac_event.target_y;
			out_event.dir_x = ppac_event.dir_x;
			out_event.dir_y = ppac_event.dir_y;

			out_event.t0d2_He1_x = d2_event.x[he1[1]];
			out_event.t0d2_He1_y = d2_event.y[he1[1]];
			out_event.t0d2_He1_z = d2_event.z[he1[1]];
			out_event.t0d2_He2_x = d2_event.x[he2[1]];
			out_event.t0d2_He2_y = d2_event.y[he2[1]];
			out_event.t0d2_He2_z = d2_event.z[he2[1]];
			out_event.t0d2_He3_x = d2_event.x[he3[1]];
			out_event.t0d2_He3_y = d2_event.y[he3[1]];
			out_event.t0d2_He3_z = d2_event.z[he3[1]];

			out_event.e_sum_He1 = e1[sort_idx[0]] + e2[sort_idx[0]] + e3[sort_idx[0]] + e4[sort_idx[0]];
			out_event.e_sum_He2 = e1[sort_idx[1]] + e2[sort_idx[1]] + e3[sort_idx[1]] + e4[sort_idx[1]];
			out_event.e_sum_He3 = e1[sort_idx[2]] + e2[sort_idx[2]] + e3[sort_idx[2]] + e4[sort_idx[2]];

			const double target_z = 0.0;
			double beam_dir[3] = {ppac_event.dir_x, ppac_event.dir_y, 1.0};
			out_event.theta_beam = brill::AngleWithZ(beam_dir);
			out_event.phi_beam = std::atan2(beam_dir[1], beam_dir[0]) * 180.0 / M_PI;

			double dir_He1[3] = {
				out_event.t0d2_He1_x - ppac_event.target_x,
				out_event.t0d2_He1_y - ppac_event.target_y,
				out_event.t0d2_He1_z - target_z
			};
			double dir_He2[3] = {
				out_event.t0d2_He2_x - ppac_event.target_x,
				out_event.t0d2_He2_y - ppac_event.target_y,
				out_event.t0d2_He2_z - target_z
			};
			double dir_He3[3] = {
				out_event.t0d2_He3_x - ppac_event.target_x,
				out_event.t0d2_He3_y - ppac_event.target_y,
				out_event.t0d2_He3_z - target_z
			};

			out_event.theta_He1 = brill::AngleBetween(dir_He1, beam_dir);
			out_event.theta_He2 = brill::AngleBetween(dir_He2, beam_dir);
			out_event.theta_He3 = brill::AngleBetween(dir_He3, beam_dir);

			out_event.opening_He1_He2 = brill::AngleBetween(dir_He1, dir_He2);
			out_event.opening_He1_He3 = brill::AngleBetween(dir_He1, dir_He3);
			out_event.opening_He2_He3 = brill::AngleBetween(dir_He2, dir_He3);

			output_tree->Fill();

			h_e1_He1_e2_He1->Fill(out_event.e2_He1, out_event.e1_He1);
			h_e2_He1_e3_He1->Fill(out_event.e3_He1, out_event.e2_He1);
			h_e3_He1_e4_He1->Fill(out_event.e4_He1, out_event.e3_He1);
			h_e1_He2_e2_He2->Fill(out_event.e2_He2, out_event.e1_He2);
			h_e2_He2_e3_He2->Fill(out_event.e3_He2, out_event.e2_He2);
			h_e3_He2_e4_He2->Fill(out_event.e4_He2, out_event.e3_He2);
			h_e1_He3_e2_He3->Fill(out_event.e2_He3, out_event.e1_He3);
			h_e2_He3_e3_He3->Fill(out_event.e3_He3, out_event.e2_He3);
			h_e3_He3_e4_He3->Fill(out_event.e4_He3, out_event.e3_He3);
		}
		printf("\b\b\b\b100%%\n");

		d1_file->Close();
		d2_file->Close();
		d3_file->Close();
		d4_file->Close();
		if (beam_file) beam_file->Close();
		if (ppac_file) ppac_file->Close();
	}

	output_file->cd();
	output_tree->Write();
	h_e1_He1_e2_He1->Write();
	h_e2_He1_e3_He1->Write();
	h_e3_He1_e4_He1->Write();
	h_e1_He2_e2_He2->Write();
	h_e2_He2_e3_He2->Write();
	h_e3_He2_e4_He2->Write();
	h_e1_He3_e2_He3->Write();
	h_e2_He3_e3_He3->Write();
	h_e3_He3_e4_He3->Write();

	output_file->Close();

	printf("Output written to %s\n", output_path.Data());
	return 0;
}