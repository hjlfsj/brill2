#include "external/cxxopts.hpp"
#include "include/config.h"
#include "include/C12_d_6Li/extract.h"
#include "include/C12_d_6Li/C12_d_6Li_event.h"
#include "include/event/beam/beam_sort.h"
#include "include/event/ingot/silicon_event.h"
#include "include/event/ppac/ppac_track_event.h"
#include "include/event/t0/dssd_match_event.h"
#include "include/physics/kinematics.h"
#include "include/utils.h"

#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TCutG.h>
#include <TROOT.h>

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <string>

int main(int argc, char **argv) {
	cxxopts::Options options("extract_6Li_two4He", "Extract C12+d->6Li+2alpha physics data");
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

	std::string match_dir = brill::JoinPath(config.workspace, config.paths.match);
	std::string track_dir = brill::JoinPath(config.workspace, config.paths.track);
	std::string beam_dir = brill::JoinPath(config.workspace, config.paths.beam);
	std::string ingot_dir = brill::JoinPath(config.workspace, config.paths.ingot);
	std::string output_dir = brill::JoinPath(config.workspace, config.paths.c12_d_6Li);
	std::string calibration_path = TString::Format(
		"%s/t0.txt",
		brill::JoinPath(config.workspace, config.paths.calibration).c_str()
	).Data();

	brill::C12D6LiCalibration calib;
	if (brill::ReadC12D6LiCalibration(calibration_path, calib)) {
		std::cerr << "Error: Read calibration from " << calibration_path << " failed.\n";
		return 1;
	}
	printf("Calibration loaded: %s\n", calibration_path.c_str());

	TString cut_li6_path = "/home/ribll2026/ribll2026_www/github_code/brill2/src/brill/Cut/cal_d1_d2_6Li_cut.C";
	gROOT->ProcessLine(TString::Format(".x %s", cut_li6_path.Data()));
	TCutG *li6_cut = (TCutG*)gROOT->FindObject("cal_d1_d2_6Li_cut");
	if (!li6_cut) {
		std::cerr << "Error: Load TCutG from " << cut_li6_path << " failed.\n";
		return 1;
	}
	printf("Cut loaded: %s\n", cut_li6_path.Data());

	std::string trigger_infix = brill::TriggerInfix(trigger);

	TString output_path = TString::Format("%s/extract_6Li_two4He_%s%04d_%04d.root",
		output_dir.c_str(), trigger_infix.c_str(), run_start, run_end);

	TFile *output_file = new TFile(output_path, "recreate");
	if (!output_file || output_file->IsZombie()) {
		std::cerr << "Error: Create output file failed: " << output_path << "\n";
		return 1;
	}

	TTree *output_tree = new TTree("tree", "C12+d->6Li+2alpha extracted events");
	brill::C12D6LiEvent out_event;
	brill::SetupOutputC12D6Li(output_tree, out_event);

	TH2D *h_e1_6Li_e2_6Li = new TH2D("h_e1_6Li_e2_6Li", "E1_6Li vs E2_6Li;E2_6Li (MeV);E1_6Li (MeV)", 1000, 0, 400, 1000, 0, 25);

	TH2D *h_e1_4He1_e2_4He1 = new TH2D("h_e1_4He1_e2_4He1", "E1_4He1 vs E2_4He1;E2_4He1 (MeV);E1_4He1 (MeV)", 1000, 0, 400, 1000, 0, 25);
	TH2D *h_e1_4He2_e2_4He2 = new TH2D("h_e1_4He2_e2_4He2", "E1_4He2 vs E2_4He2;E2_4He2 (MeV);E1_4He2 (MeV)", 1000, 0, 400, 1000, 0, 25);

	TH2D *h_e2_4He1_e3_4He1 = new TH2D("h_e2_4He1_e3_4He1", "E2_4He1 vs E3_4He1;E3_4He1 (MeV);E2_4He1 (MeV)", 1000, 0, 350, 1000, 0, 400);
	TH2D *h_e2_4He2_e3_4He2 = new TH2D("h_e2_4He2_e3_4He2", "E2_4He2 vs E3_4He2;E3_4He2 (MeV);E2_4He2 (MeV)", 1000, 0, 350, 1000, 0, 400);

	TH2D *h_e3_4He1_e4_4He1 = new TH2D("h_e3_4He1_e4_4He1", "E3_4He1 vs E4_4He1;E4_4He1 (MeV);E3_4He1 (MeV)", 1000, 0, 250, 1000, 0, 300);
	TH2D *h_e3_4He2_e4_4He2 = new TH2D("h_e3_4He2_e4_4He2", "E3_4He2 vs E4_4He2;E4_4He2 (MeV);E3_4He2 (MeV)", 1000, 0, 250, 1000, 0, 300);

	TH2D *h_e4sum_e5 = new TH2D("h_e4sum_e5", "E4_4He1+E4_4He2 vs E5;E5 (MeV);E4_4He1+E4_4He2 (MeV)", 1000, 0, 500, 1000, 0, 500);

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
		TString t0s_path = TString::Format("%s/t0s_%s%04d.root",
			ingot_dir.c_str(), trigger_infix.c_str(), run);
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

		TFile *d2_file = nullptr;
		TTree *d2_tree = nullptr;
		brill::DssdMatchEvent d2_event;
		if (std::filesystem::exists(d2_path.Data())) {
			d2_file = new TFile(d2_path, "read");
			d2_tree = (TTree*)d2_file->Get("tree");
			if (d2_tree) brill::SetupInput(d2_tree, d2_event, "");
		}

		TFile *d3_file = nullptr;
		TTree *d3_tree = nullptr;
		brill::DssdMatchEvent d3_event;
		if (std::filesystem::exists(d3_path.Data())) {
			d3_file = new TFile(d3_path, "read");
			d3_tree = (TTree*)d3_file->Get("tree");
			if (d3_tree) brill::SetupInput(d3_tree, d3_event, "");
		}

		TFile *d4_file = nullptr;
		TTree *d4_tree = nullptr;
		brill::DssdMatchEvent d4_event;
		if (std::filesystem::exists(d4_path.Data())) {
			d4_file = new TFile(d4_path, "read");
			d4_tree = (TTree*)d4_file->Get("tree");
			if (d4_tree) brill::SetupInput(d4_tree, d4_event, "");
		}

		TFile *t0s_file = nullptr;
		TTree *t0s_tree = nullptr;
		brill::SiliconEvent t0s_event;
		if (std::filesystem::exists(t0s_path.Data())) {
			t0s_file = new TFile(t0s_path, "read");
			t0s_tree = (TTree*)t0s_file->Get("tree");
			if (t0s_tree) brill::SetupInput(t0s_tree, t0s_event, "");
		}

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
			if (t0s_tree) t0s_tree->GetEntry(i);
			if (beam_tree) beam_tree->GetEntry(i);
			if (ppac_tree) ppac_tree->GetEntry(i);

			if (!brill::PassC12D6LiCut(d1_event, d2_event, d3_event, d4_event, calib, li6_cut)) continue;

			brill::Two4HeResult cls = brill::ClassifyTwo4He(d1_event, d2_event, d3_event, d4_event, calib, li6_cut);
			if (!cls.passed) continue;

			double e1_6Li = brill::CalibrateC12D6LiEnergy(calib, 0, d1_event.energy[cls.idx_6Li_d1]);
			double e2_6Li = brill::CalibrateC12D6LiEnergy(calib, 1, d2_event.energy[cls.idx_6Li_d2]);

			double e1_4He1 = brill::CalibrateC12D6LiEnergy(calib, 0, d1_event.energy[cls.e1_4He1_idx]);
			double e1_4He2 = brill::CalibrateC12D6LiEnergy(calib, 0, d1_event.energy[cls.e1_4He2_idx]);

			double e2_4He1 = brill::CalibrateC12D6LiEnergy(calib, 1, d2_event.energy[cls.e2_4He1_idx]);
			double e2_4He2 = brill::CalibrateC12D6LiEnergy(calib, 1, d2_event.energy[cls.e2_4He2_idx]);

			double e3_4He1 = brill::CalibrateC12D6LiEnergy(calib, 2, d3_event.energy[cls.e3_4He1_idx]);
			double e3_4He2 = brill::CalibrateC12D6LiEnergy(calib, 2, d3_event.energy[cls.e3_4He2_idx]);

			double e4_4He1 = brill::CalibrateC12D6LiEnergy(calib, 3, d4_event.energy[cls.e4_4He1_idx]);
			double e4_4He2 = brill::CalibrateC12D6LiEnergy(calib, 3, d4_event.energy[cls.e4_4He2_idx]);

			double e5 = brill::CalibrateC12D6LiEnergy(calib, 4, t0s_event.energy);

			out_event.run_number = run;
			out_event.entry = i;
			out_event.e1_6Li = e1_6Li;
			out_event.e2_6Li = e2_6Li;
			out_event.e1_4He1 = e1_4He1;
			out_event.e2_4He1 = e2_4He1;
			out_event.e3_4He1 = e3_4He1;
			out_event.e4_4He1 = e4_4He1;
			out_event.e1_4He2 = e1_4He2;
			out_event.e2_4He2 = e2_4He2;
			out_event.e3_4He2 = e3_4He2;
			out_event.e4_4He2 = e4_4He2;
			out_event.e5 = e5;
			out_event.is_14O = beam_is_14O;
			out_event.is_13N = beam_is_13N;
			out_event.is_12C = beam_is_12C;
			out_event.ppac_valid = (ppac_event.valid != 0);
			out_event.target_x = ppac_event.target_x;
			out_event.target_y = ppac_event.target_y;
			out_event.dir_x = ppac_event.dir_x;
			out_event.dir_y = ppac_event.dir_y;
			out_event.t0d2_6Li_x = d2_event.x[cls.idx_6Li_d2];
			out_event.t0d2_6Li_y = d2_event.y[cls.idx_6Li_d2];
			out_event.t0d2_6Li_z = d2_event.z[cls.idx_6Li_d2];
			out_event.t0d2_4He1_x = d2_event.x[cls.e2_4He1_idx];
			out_event.t0d2_4He1_y = d2_event.y[cls.e2_4He1_idx];
			out_event.t0d2_4He1_z = d2_event.z[cls.e2_4He1_idx];
			out_event.t0d2_4He2_x = d2_event.x[cls.e2_4He2_idx];
			out_event.t0d2_4He2_y = d2_event.y[cls.e2_4He2_idx];
			out_event.t0d2_4He2_z = d2_event.z[cls.e2_4He2_idx];

			const double target_z = 0.0;
			if (ppac_event.valid) {
				double beam_dir[3] = {ppac_event.dir_x, ppac_event.dir_y, 1.0};
				out_event.theta_beam = brill::AngleWithZ(beam_dir);
				out_event.phi_beam = std::atan2(beam_dir[1], beam_dir[0]) * 180.0 / M_PI;

				double dir_6Li[3] = {
					out_event.t0d2_6Li_x - ppac_event.target_x,
					out_event.t0d2_6Li_y - ppac_event.target_y,
					out_event.t0d2_6Li_z - target_z
				};
				double dir_4He1[3] = {
					out_event.t0d2_4He1_x - ppac_event.target_x,
					out_event.t0d2_4He1_y - ppac_event.target_y,
					out_event.t0d2_4He1_z - target_z
				};
				double dir_4He2[3] = {
					out_event.t0d2_4He2_x - ppac_event.target_x,
					out_event.t0d2_4He2_y - ppac_event.target_y,
					out_event.t0d2_4He2_z - target_z
				};

				out_event.theta_6Li = brill::AngleBetween(dir_6Li, beam_dir);
				out_event.theta_4He1 = brill::AngleBetween(dir_4He1, beam_dir);
				out_event.theta_4He2 = brill::AngleBetween(dir_4He2, beam_dir);
				out_event.opening_6Li_4He1 = brill::AngleBetween(dir_6Li, dir_4He1);
				out_event.opening_6Li_4He2 = brill::AngleBetween(dir_6Li, dir_4He2);
				out_event.opening_4He1_4He2 = brill::AngleBetween(dir_4He1, dir_4He2);
			} else {
				out_event.theta_beam = 0.0;
				out_event.phi_beam = 0.0;
				out_event.theta_6Li = 0.0;
				out_event.theta_4He1 = 0.0;
				out_event.theta_4He2 = 0.0;
				out_event.opening_6Li_4He1 = 0.0;
				out_event.opening_6Li_4He2 = 0.0;
				out_event.opening_4He1_4He2 = 0.0;
			}

			output_tree->Fill();

			h_e1_6Li_e2_6Li->Fill(e2_6Li, e1_6Li);

			h_e1_4He1_e2_4He1->Fill(e2_4He1, e1_4He1);
			h_e1_4He2_e2_4He2->Fill(e2_4He2, e1_4He2);

			h_e2_4He1_e3_4He1->Fill(e3_4He1, e2_4He1);
			h_e2_4He2_e3_4He2->Fill(e3_4He2, e2_4He2);

			h_e3_4He1_e4_4He1->Fill(e4_4He1, e3_4He1);
			h_e3_4He2_e4_4He2->Fill(e4_4He2, e3_4He2);

			h_e4sum_e5->Fill(e5, e4_4He1 + e4_4He2);
		}
		printf("\b\b\b\b100%%\n");

		d1_file->Close();
		if (d2_file) d2_file->Close();
		if (d3_file) d3_file->Close();
		if (d4_file) d4_file->Close();
		if (t0s_file) t0s_file->Close();
		if (beam_file) beam_file->Close();
		if (ppac_file) ppac_file->Close();
	}

	output_file->cd();
	output_tree->Write();
	h_e1_6Li_e2_6Li->Write();
	h_e1_4He1_e2_4He1->Write();
	h_e1_4He2_e2_4He2->Write();
	h_e2_4He1_e3_4He1->Write();
	h_e2_4He2_e3_4He2->Write();
	h_e3_4He1_e4_4He1->Write();
	h_e3_4He2_e4_4He2->Write();
	h_e4sum_e5->Write();

	output_file->Close();

	printf("Output written to %s\n", output_path.Data());
	return 0;
}