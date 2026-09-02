#include "external/cxxopts.hpp"
#include "include/config.h"
#include "include/10C+4He/extract.h"
#include "include/10C+4He/10C_4He_event.h"
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
	cxxopts::Options options("extract_10C_4He", "Extract 10C+4He physics data");
	options.add_options()
		("h,help", "Print help information.")
		("r,run", "Start run number.", cxxopts::value<int>())
		("e,end-run", "End run number.", cxxopts::value<int>())
		("t,trigger", "Trigger type.", cxxopts::value<std::string>())
		("c,config", "Config file path.", cxxopts::value<std::string>()->default_value("config.toml"))
		("d1_hit", "d1 hit count: 1=10C only (4He lost on d1), 2=both (default).", cxxopts::value<int>()->default_value("2"))
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
	int d1_hit = result["d1_hit"].as<int>();

	brill::AppConfig config;
	if (brill::LoadConfig(config_path, config)) {
		std::cerr << "Error: Load config failed.\n";
		return 1;
	}

	std::string match_dir = brill::JoinPath(config.workspace, config.paths.match);
	std::string track_dir = brill::JoinPath(config.workspace, config.paths.track);
	std::string beam_dir = brill::JoinPath(config.workspace, config.paths.beam);
	std::string ingot_dir = brill::JoinPath(config.workspace, config.paths.ingot);
	std::string output_dir = brill::JoinPath(config.workspace, config.paths.c10_he4);
	std::string calibration_path = TString::Format(
		"%s/t0.txt",
		brill::JoinPath(config.workspace, config.paths.calibration).c_str()
	).Data();

	brill::C10He4Calibration calib;
	if (brill::ReadC10He4Calibration(calibration_path, calib)) {
		std::cerr << "Error: Read calibration from " << calibration_path << " failed.\n";
		return 1;
	}
	printf("Calibration loaded: %s\n", calibration_path.c_str());

	TString cut_d2d3_path = "/home/ribll2026/ribll2026_www/github_code/brill2/src/brill/Cut/cal_d2_d3_stop_10C_cut.C";
	gROOT->ProcessLine(TString::Format(".x %s", cut_d2d3_path.Data()));
	TCutG *d2d3_cut = (TCutG*)gROOT->FindObject("cal_d2_d3_stop_10C_cut");
	if (!d2d3_cut) {
		std::cerr << "Error: Load TCutG from " << cut_d2d3_path << " failed.\n";
		return 1;
	}
	printf("Cut loaded: %s\n", cut_d2d3_path.Data());

	std::string trigger_infix = brill::TriggerInfix(trigger);

	TString d1_mark = (d1_hit == 1) ? "_d1hit1" : "";
	TString output_path = TString::Format("%s/extract_10C_4He_%s%04d_%04d%s.root",
		output_dir.c_str(), trigger_infix.c_str(), run_start, run_end, d1_mark.Data());

	TFile *output_file = new TFile(output_path, "recreate");
	if (!output_file || output_file->IsZombie()) {
		std::cerr << "Error: Create output file failed: " << output_path << "\n";
		return 1;
	}

	TTree *output_tree = new TTree("tree", "10C+4He extracted events");
	brill::C10He4Event out_event;
	brill::SetupOutputC10He4(output_tree, out_event);

	TH2D *h_e1_10C_e2_10C = new TH2D("h_e1_10C_e2_10C", "E1_10C vs E2_10C;E2_10C (MeV);E1_10C (MeV)", 1000, 0, 400, 1000, 0, 25);
	TH2D *h_e2_10C_e3_10C = new TH2D("h_e2_10C_e3_10C", "E2_10C vs E3_10C;E3_10C (MeV);E2_10C (MeV)", 1000, 0, 350, 1000, 0, 400);

	TH2D *h_e1_4He_e2_4He = new TH2D("h_e1_4He_e2_4He", "E1_4He vs E2_4He;E2_4He (MeV);E1_4He (MeV)", 1000, 0, 400, 1000, 0, 25);
	TH2D *h_e2_4He_e3_4He = new TH2D("h_e2_4He_e3_4He", "E2_4He vs E3_4He;E3_4He (MeV);E2_4He (MeV)", 1000, 0, 350, 1000, 0, 400);
	TH2D *h_e3_4He_e4_4He = new TH2D("h_e3_4He_e4_4He", "E3_4He vs E4_4He;E4_4He (MeV);E3_4He (MeV)", 1000, 0, 250, 1000, 0, 300);
	TH2D *h_e4_4He_e5_4He = new TH2D("h_e4_4He_e5_4He", "E4_4He vs E5_4He;E5_4He (MeV);E4_4He (MeV)", 1000, 0, 250, 1000, 0, 250);

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

			if (!brill::Pass10C_d3_4He_s1Cut(d1_event, d2_event, d3_event, d4_event, t0s_event, calib, d2d3_cut, d1_hit)) continue;

			int idx_10C_d1 = 0;
			int idx_4He_d1 = 1;
			int idx_10C_d2 = 0;
			int idx_4He_d2 = 1;
			int idx_10C_d3 = 0;
			int idx_4He_d3 = 1;

			if (d1_hit == 2) {
				if (d1_event.energy[1] > d1_event.energy[0]) {
					idx_10C_d1 = 1;
					idx_4He_d1 = 0;
				}
			}
			if (d2_event.energy[1] > d2_event.energy[0]) {
				idx_10C_d2 = 1;
				idx_4He_d2 = 0;
			}
			if (d3_event.energy[1] > d3_event.energy[0]) {
				idx_10C_d3 = 1;
				idx_4He_d3 = 0;
			}

			double e1_10C = brill::CalibrateC10He4Energy(calib, 0, d1_event.energy[idx_10C_d1]);
			double e2_10C = brill::CalibrateC10He4Energy(calib, 1, d2_event.energy[idx_10C_d2]);
			double e3_10C = brill::CalibrateC10He4Energy(calib, 2, d3_event.energy[idx_10C_d3]);

			double e1_4He = (d1_hit == 2) ? brill::CalibrateC10He4Energy(calib, 0, d1_event.energy[idx_4He_d1]) : 0.0;
			double e2_4He = brill::CalibrateC10He4Energy(calib, 1, d2_event.energy[idx_4He_d2]);
			double e3_4He = brill::CalibrateC10He4Energy(calib, 2, d3_event.energy[idx_4He_d3]);
			double e4_4He = brill::CalibrateC10He4Energy(calib, 3, d4_event.energy[0]);
			double e5_4He = brill::CalibrateC10He4Energy(calib, 4, t0s_event.energy);

			out_event.run_number = run;
			out_event.entry = i;
			out_event.e1_10C = e1_10C;
			out_event.e2_10C = e2_10C;
			out_event.e3_10C = e3_10C;
			out_event.e1_4He = e1_4He;
			out_event.e2_4He = e2_4He;
			out_event.e3_4He = e3_4He;
			out_event.e4_4He = e4_4He;
			out_event.e5_4He = e5_4He;
			out_event.is_14O = beam_is_14O;
			out_event.is_13N = beam_is_13N;
			out_event.is_12C = beam_is_12C;
			out_event.ppac_valid = (ppac_event.valid != 0);
			out_event.target_x = ppac_event.target_x;
			out_event.target_y = ppac_event.target_y;
			out_event.dir_x = ppac_event.dir_x;
			out_event.dir_y = ppac_event.dir_y;
			out_event.t0d2_10C_x = d2_event.x[idx_10C_d2];
			out_event.t0d2_10C_y = d2_event.y[idx_10C_d2];
			out_event.t0d2_10C_z = d2_event.z[idx_10C_d2];
			out_event.t0d2_4He_x = d2_event.x[idx_4He_d2];
			out_event.t0d2_4He_y = d2_event.y[idx_4He_d2];
			out_event.t0d2_4He_z = d2_event.z[idx_4He_d2];

			const double target_z = 0.0;
			if (ppac_event.valid) {
				double beam_dir[3] = {ppac_event.dir_x, ppac_event.dir_y, 1.0};
				out_event.theta_beam = brill::AngleWithZ(beam_dir);
				out_event.phi_beam = std::atan2(beam_dir[1], beam_dir[0]) * 180.0 / M_PI;

				double dir_10C[3] = {
					out_event.t0d2_10C_x - ppac_event.target_x,
					out_event.t0d2_10C_y - ppac_event.target_y,
					out_event.t0d2_10C_z - target_z
				};
				double dir_4He[3] = {
					out_event.t0d2_4He_x - ppac_event.target_x,
					out_event.t0d2_4He_y - ppac_event.target_y,
					out_event.t0d2_4He_z - target_z
				};

				out_event.theta_10C = brill::AngleBetween(dir_10C, beam_dir);
				out_event.theta_4He = brill::AngleBetween(dir_4He, beam_dir);
				out_event.opening_angle = brill::AngleBetween(dir_10C, dir_4He);
			} else {
				out_event.theta_beam = 0.0;
				out_event.phi_beam = 0.0;
				out_event.theta_10C = 0.0;
				out_event.theta_4He = 0.0;
				out_event.opening_angle = 0.0;
			}

			output_tree->Fill();

			h_e1_10C_e2_10C->Fill(e2_10C, e1_10C);
			h_e2_10C_e3_10C->Fill(e3_10C, e2_10C);

			h_e1_4He_e2_4He->Fill(e2_4He, e1_4He);
			h_e2_4He_e3_4He->Fill(e3_4He, e2_4He);
			h_e3_4He_e4_4He->Fill(e4_4He, e3_4He);
			h_e4_4He_e5_4He->Fill(e5_4He, e4_4He);
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
	h_e1_10C_e2_10C->Write();
	h_e2_10C_e3_10C->Write();
	h_e1_4He_e2_4He->Write();
	h_e2_4He_e3_4He->Write();
	h_e3_4He_e4_4He->Write();
	h_e4_4He_e5_4He->Write();

	output_file->Close();

	printf("Output written to %s\n", output_path.Data());
	return 0;
}