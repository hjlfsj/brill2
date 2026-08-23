#include "external/cxxopts.hpp"
#include "include/config.h"
#include "include/d_6Li/extract.h"
#include "include/d_6Li/d_6Li_event.h"
#include "include/event/t0/dssd_match_event.h"
#include "include/event/beam/beam_sort.h"
#include "include/event/ppac/ppac_track_event.h"
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
	cxxopts::Options options("extract_d_Li6", "Extract d+6Li physics data");
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
	std::string output_dir = brill::JoinPath(config.workspace, config.paths.d_Li6);
	std::string calibration_path = TString::Format(
		"%s/t0.txt",
		brill::JoinPath(config.workspace, config.paths.calibration).c_str()
	).Data();

	brill::D6LiCalibration calib;
	if (brill::ReadD6LiCalibration(calibration_path, calib)) {
		std::cerr << "Error: Read calibration from " << calibration_path << " failed.\n";
		return 1;
	}
	printf("Calibration loaded: %s\n", calibration_path.c_str());

	TString cut_path = "/home/ribll2026/ribll2026_www/github_code/brill2/src/brill/Cut/cal_d3_d4_10C_cut.C";
	gROOT->ProcessLine(TString::Format(".x %s", cut_path.Data()));
	TCutG *d3d4_cut = (TCutG*)gROOT->FindObject("cal_d3_d4_10C_cut");
	if (!d3d4_cut) {
		std::cerr << "Error: Load TCutG from " << cut_path << " failed.\n";
		return 1;
	}
	printf("Cut loaded: %s\n", cut_path.Data());

	TString cut_d2d3_path = "/home/ribll2026/ribll2026_www/github_code/brill2/src/brill/Cut/cal_d2_d3_10C_cut.C";
	gROOT->ProcessLine(TString::Format(".x %s", cut_d2d3_path.Data()));
	TCutG *d2d3_cut = (TCutG*)gROOT->FindObject("cal_d2_d3_10C_cut");
	if (!d2d3_cut) {
		std::cerr << "Error: Load TCutG from " << cut_d2d3_path << " failed.\n";
		return 1;
	}
	printf("Cut loaded: %s\n", cut_d2d3_path.Data());

	std::string trigger_infix = brill::TriggerInfix(trigger);

	TString output_path = TString::Format("%s/extract_d_Li6_%s%04d_%04d.root",
		output_dir.c_str(), trigger_infix.c_str(), run_start, run_end);

	TFile *output_file = new TFile(output_path, "recreate");
	if (!output_file || output_file->IsZombie()) {
		std::cerr << "Error: Create output file failed: " << output_path << "\n";
		return 1;
	}

	TTree *output_tree = new TTree("tree", "d+6Li extracted events");
	brill::D6LiEvent out_event;
	brill::SetupOutputD6Li(output_tree, out_event);

	TH2D *h_d1d2 = new TH2D("h_d1d2", "D1 vs D2 energy;D2 E_{1} (MeV);D1 E_{1} (MeV)", 1000, 0, 400, 1000, 0, 25);
	TH2D *h_d2d3 = new TH2D("h_d2d3", "D2 vs D3 energy;D3 E_{1} (MeV);D2 E_{1} (MeV)", 1000, 0, 350, 1000, 0, 400);
	TH2D *h_d3d4 = new TH2D("h_d3d4", "D3 vs D4 energy;D4 E_{1} (MeV);D3 E_{1} (MeV)", 1000, 0, 250, 1000, 0, 300);

	TH2D *h_e1_10C_e2_10C = new TH2D("h_e1_10C_e2_10C", "E1_10C vs E2_10C;E2_10C (MeV);E1_10C (MeV)", 1000, 0, 400, 1000, 0, 25);
	TH2D *h_e1_6Li_e2_6Li = new TH2D("h_e1_6Li_e2_6Li", "E1_6Li vs E2_6Li;E2_6Li (MeV);E1_6Li (MeV)", 1000, 0, 400, 1000, 0, 25);
	TH2D *h_e2_10C_e3_10C = new TH2D("h_e2_10C_e3_10C", "E2_10C vs E3_10C;E3_10C (MeV);E2_10C (MeV)", 1000, 0, 350, 1000, 0, 400);
	TH2D *h_e3_10C_e4_10C = new TH2D("h_e3_10C_e4_10C", "E3_10C vs E4_10C;E4_10C (MeV);E3_10C (MeV)", 1000, 0, 250, 1000, 0, 300);

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
			if (beam_tree) beam_tree->GetEntry(i);
			if (ppac_tree) ppac_tree->GetEntry(i);

			if (!brill::PassD6LiCut(d1_event, d2_event, d3_event, d4_event, calib, d3d4_cut)) continue;

			brill::D6LiAdvancedResult cls = brill::ClassifyD6Li(d1_event, d2_event, d3_event, d4_event, calib, d2d3_cut);
			if (!cls.passed) continue;

			double e1 = brill::CalibrateD6LiEnergy(calib, 0, d1_event.energy[0]);
			double e2 = brill::CalibrateD6LiEnergy(calib, 1, d2_event.energy[0]);
			double e3 = brill::CalibrateD6LiEnergy(calib, 2, d3_event.energy[0]);
			double e4 = brill::CalibrateD6LiEnergy(calib, 3, d4_event.energy[0]);

			double e1_10C = brill::CalibrateD6LiEnergy(calib, 0, d1_event.energy[cls.e1_10C_idx]);
			double e1_6Li = brill::CalibrateD6LiEnergy(calib, 0, d1_event.energy[cls.e1_6Li_idx]);
			double e2_10C = brill::CalibrateD6LiEnergy(calib, 1, d2_event.energy[cls.e2_10C_idx]);
			double e2_6Li = brill::CalibrateD6LiEnergy(calib, 1, d2_event.energy[cls.e2_6Li_idx]);
			double e3_10C = brill::CalibrateD6LiEnergy(calib, 2, d3_event.energy[0]);
			double e4_10C = brill::CalibrateD6LiEnergy(calib, 3, d4_event.energy[0]);

			out_event.run_number = run;
			out_event.entry = i;
			out_event.e1 = e1;
			out_event.e2 = e2;
			out_event.e3 = e3;
			out_event.e4 = e4;
			out_event.e1_10C = e1_10C;
			out_event.e2_10C = e2_10C;
			out_event.e3_10C = e3_10C;
			out_event.e4_10C = e4_10C;
			out_event.e1_6Li = e1_6Li;
			out_event.e2_6Li = e2_6Li;
			out_event.is_14O = beam_is_14O;
			out_event.is_13N = beam_is_13N;
			out_event.is_12C = beam_is_12C;
			out_event.ppac_valid = (ppac_event.valid != 0);
			out_event.target_x = ppac_event.target_x;
			out_event.target_y = ppac_event.target_y;
			out_event.dir_x = ppac_event.dir_x;
			out_event.dir_y = ppac_event.dir_y;
			out_event.t0d2_10C_x = d2_event.x[cls.e2_10C_idx];
			out_event.t0d2_10C_y = d2_event.y[cls.e2_10C_idx];
			out_event.t0d2_10C_z = d2_event.z[cls.e2_10C_idx];
			out_event.t0d2_6Li_x = d2_event.x[cls.e2_6Li_idx];
			out_event.t0d2_6Li_y = d2_event.y[cls.e2_6Li_idx];
			out_event.t0d2_6Li_z = d2_event.z[cls.e2_6Li_idx];

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
				double dir_6Li[3] = {
					out_event.t0d2_6Li_x - ppac_event.target_x,
					out_event.t0d2_6Li_y - ppac_event.target_y,
					out_event.t0d2_6Li_z - target_z
				};

				out_event.theta_10C = brill::AngleBetween(dir_10C, beam_dir);
				out_event.theta_6Li = brill::AngleBetween(dir_6Li, beam_dir);
				out_event.opening_angle = brill::AngleBetween(dir_10C, dir_6Li);
			} else {
				out_event.theta_beam = 0.0;
				out_event.phi_beam = 0.0;
				out_event.theta_10C = 0.0;
				out_event.theta_6Li = 0.0;
				out_event.opening_angle = 0.0;
			}

			output_tree->Fill();

			h_d1d2->Fill(e2, e1);
			h_d2d3->Fill(e3, e2);
			h_d3d4->Fill(e4, e3);

			h_e1_10C_e2_10C->Fill(e2_10C, e1_10C);
			h_e1_6Li_e2_6Li->Fill(e2_6Li, e1_6Li);
			h_e2_10C_e3_10C->Fill(e3_10C, e2_10C);
			h_e3_10C_e4_10C->Fill(e4_10C, e3_10C);
		}
		printf("\b\b\b\b100%%\n");

		d1_file->Close();
		if (d2_file) d2_file->Close();
		if (d3_file) d3_file->Close();
		if (d4_file) d4_file->Close();
		if (beam_file) beam_file->Close();
		if (ppac_file) ppac_file->Close();
	}

	output_file->cd();
	output_tree->Write();
	h_d1d2->Write();
	h_d2d3->Write();
	h_d3d4->Write();
	h_e1_10C_e2_10C->Write();
	h_e1_6Li_e2_6Li->Write();
	h_e2_10C_e3_10C->Write();
	h_e3_10C_e4_10C->Write();

	output_file->Close();

	printf("Output written to %s\n", output_path.Data());
	return 0;
}