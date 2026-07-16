#include <filesystem>
#include <iostream>
#include <string>

#include <TH2D.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TString.h>

#include "external/cxxopts.hpp"
#include "include/config.h"
#include "include/event/ingot/ppac_event.h"
#include "include/event/ppac/ppac_track.h"
#include "include/event/ppac/ppac_track_event.h"
#include "include/utils.h"

void PrintUsage(const cxxopts::Options &options) {
	std::cout << options.help() << "\n";
}

int main(int argc, char **argv) {
	cxxopts::Options options("track_ppac", "PPAC track reconstruction.");
	options.add_options()
		("h,help", "Print help information.")
		("r,run", "Run number.", cxxopts::value<int>(), "run")
		("t,trigger", "Trigger type.", cxxopts::value<std::string>()->default_value("main"), "trigger")
		("c,config", "Config file path.", cxxopts::value<std::string>()->default_value("config.toml"), "file")
		("draw", "Enable track line drawing histograms.");

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

	std::string trigger = result["trigger"].as<std::string>();
	config.trigger = trigger;
	const int run = result["run"].as<int>();
	bool draw_track = result.count("draw") > 0;

	std::string input_dir = brill::JoinPath(config.workspace, config.paths.ingot);
	std::string normalize_dir = brill::JoinPath(config.workspace, config.paths.normalize);
	std::string output_dir = brill::JoinPath(config.workspace, config.paths.track);
	std::filesystem::create_directories(output_dir);

	std::string prefix = (trigger == "main") ? "" : "t1_";

	TString input_path = TString::Format(
		"%s/ppac_%s%04d.root",
		input_dir.c_str(),
		brill::TriggerInfix(config.trigger).c_str(),
		run
	);
	TString offset_path = TString::Format(
		"%s/ppac_offset_%04d.txt",
		normalize_dir.c_str(),
		run
	);
	TString output_path = TString::Format(
		"%s/ppac_%s%04d.root",
		output_dir.c_str(),
		prefix.c_str(),
		run
	);

	std::cout << "Input:  " << input_path << "\n";
	std::cout << "Offset: " << offset_path << "\n";
	std::cout << "Output: " << output_path << "\n";

	brill::PpacOffsetParams offset[3];
	if (brill::ReadPpacOffsetParams(offset_path.Data(), offset)) {
		return 1;
	}
	for (int i = 0; i < 3; ++i) {
		std::cout << "PPAC" << (i + 1) << " offset: "
		          << "x_peak=" << offset[i].delay_x_peak
		          << " x_sigma=" << offset[i].delay_x_sigma
		          << " p0_x=" << offset[i].position_x_p0
		          << " p1_x=" << offset[i].position_x_p1
		          << " y_peak=" << offset[i].delay_y_peak
		          << " y_sigma=" << offset[i].delay_y_sigma
		          << " p0_y=" << offset[i].position_y_p0
		          << " p1_y=" << offset[i].position_y_p1 << "\n";
	}

	TChain chain("tree");
	if (std::filesystem::exists(input_path.Data())) {
		chain.Add(input_path);
	} else {
		std::cerr << "Error: Input file not found: " << input_path << "\n";
		return 1;
	}

	brill::PpacEvent event;
	brill::SetupInput(&chain, event, "");

	TFile opf(output_path, "recreate");
	TTree opt("tree", "ppac track");
	brill::PpacTrackEvent track;
	brill::SetupOutput(&opt, track);

	TH2D *h_track_zx = nullptr;
	TH2D *h_track_zy = nullptr;
	if (draw_track) {
		h_track_zx = new TH2D("h_track_zx", "Track Z-X;z (mm);x (mm)", 1500, -1000, 500, 100, -50, 50);
		h_track_zy = new TH2D("h_track_zy", "Track Z-Y;z (mm);y (mm)", 1500, -1000, 500, 100, -50, 50);
	}

	long long total = chain.GetEntries();
	long long last_percentage = 0;
	printf("Tracking  0%%");
	fflush(stdout);
	for (long long entry = 0; entry < total; ++entry) {
		if (entry * 100ll / total > last_percentage) {
			last_percentage = entry * 100ll / total;
			printf("\b\b\b\b%3lld%%", last_percentage);
			fflush(stdout);
		}
		chain.GetEntry(entry);
		brill::TrackPpac(event, config.ppac, offset, track);
		if (draw_track && track.valid) {
			for (int iz = 0; iz < 1500; ++iz) {
				double z = h_track_zx->GetXaxis()->GetBinCenter(iz + 1);
				double x = track.target_x + track.dir_x * z;
				double y = track.target_y + track.dir_y * z;
				h_track_zx->Fill(z, x);
				h_track_zy->Fill(z, y);
			}
		}
		opt.Fill();
	}
	printf("\b\b\b\b100%%\n");

	opf.cd();
	opt.Write();
	if (draw_track) {
		h_track_zx->Write();
		h_track_zy->Write();
		delete h_track_zx;
		delete h_track_zy;
	}
	opf.Close();

	return 0;
}