#include <filesystem>
#include <iostream>
#include <string>

#include <TString.h>

#include "external/cxxopts.hpp"
#include "include/config.h"
#include "include/event/t0/adjust_t0.h"
#include "include/utils.h"

void PrintUsage(const cxxopts::Options &options) {
	std::cout << options.help() << "\n";
}

int main(int argc, char **argv) {
	cxxopts::Options options("adjust_t0_step1", "Adjust T0D1/D2/D3 positions using beam events.");
	options.add_options()
		("h,help", "Print help information.")
		("r,run", "Run number.", cxxopts::value<int>(), "run")
		("t,trigger", "Trigger type.", cxxopts::value<std::string>()->default_value("t1"), "trigger")
		("c,config", "Config file path.",
			cxxopts::value<std::string>()->default_value("config.toml"), "file");

	auto result = options.parse(argc, argv);
	if (result.count("help")) { PrintUsage(options); return 0; }
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

	std::string track_dir = brill::JoinPath(config.workspace, config.paths.track);
	std::string match_dir = brill::JoinPath(config.workspace, config.paths.match);
	std::string normalize_dir = brill::JoinPath(config.workspace, config.paths.normalize);
	std::string track_path = TString::Format(
		"%s/ppac_%s%04d.root",
		track_dir.c_str(),
		brill::TriggerInfix(trigger).c_str(),
		run
	).Data();

	std::cout << "Track file: " << track_path << "\n";
	std::cout << "Match dir:  " << match_dir << "\n";

	if (!std::filesystem::exists(track_path)) {
		std::cerr << "Error: Track file not found: " << track_path << "\n";
		return 1;
	}

	std::vector<std::string> detector_names = {"t0d1", "t0d2", "t0d3"};
	std::vector<brill::SquareDetectorConfig> detectors;
	for (const auto &name : detector_names) {
		const auto *det = brill::FindDetectorConfig(config, name);
		if (!det) {
			std::cerr << "Error: Detector " << name << " not found in config.\n";
			return 1;
		}
		detectors.push_back(*det);
	}

	brill::T0AdjustConfig adjust_config;
	adjust_config.track_path = track_path;
	adjust_config.match_dir = match_dir;
	adjust_config.trigger = trigger;
	adjust_config.run = run;
	adjust_config.output_dir = normalize_dir;
	adjust_config.detector_names = detector_names;
	adjust_config.energy_cuts = {
		{0.0, 0.0},           // t0d1: no energy cut
		{43000.0, 47000.0},   // t0d2: 43000 < E < 47000
		{30000.0, 35000.0}    // t0d3: 30000 < E < 35000
	};

	std::vector<brill::T0AdjustResult> results;
	if (brill::AdjustT0Step1(adjust_config, detectors, results)) {
		return 1;
	}

	std::cout << "\n=== Results ===\n";
	for (const auto &r : results) {
		std::cout << "--- " << r.detector_name << " ---\n";
		std::cout << "  Selected events: " << r.num_events << "\n";
		std::cout << "  X residual mean: " << r.residual_x_mean
			<< "  sigma: " << r.residual_x_sigma << "\n";
		if (r.dsigma_x >= 0) {
			std::cout << "  X offset: " << r.dx << " +/- " << r.dsigma_x << " mm\n";
		} else {
			std::cout << "  X offset: " << r.dx << " mm (not well-determined)\n";
		}
		std::cout << "  Y residual mean: " << r.residual_y_mean
			<< "  sigma: " << r.residual_y_sigma << "\n";
		if (r.dsigma_y >= 0) {
			std::cout << "  Y offset: " << r.dy << " +/- " << r.dsigma_y << " mm\n";
		} else {
			std::cout << "  Y offset: " << r.dy << " mm (not well-determined)\n";
		}
		std::cout << "\n";
	}

	return 0;
}