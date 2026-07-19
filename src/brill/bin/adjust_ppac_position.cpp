#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

#include <TFile.h>
#include <TString.h>

#include "external/cxxopts.hpp"
#include "include/config.h"
#include "include/event/ppac/adjust_ppac.h"
#include "include/event/ppac/ppac_track_event.h"
#include "include/utils.h"

void PrintUsage(const cxxopts::Options &options) {
	std::cout << options.help() << "\n";
}

int main(int argc, char **argv) {
	cxxopts::Options options("adjust_ppac_position", "Adjust PPAC2/3 positions via grid search.");
	options.add_options()
		("h,help", "Print help information.")
		("r,run", "Run number.", cxxopts::value<int>(), "run")
		("t,trigger", "Trigger type.", cxxopts::value<std::string>()->default_value("main"), "trigger")
		("c,config", "Config file path.",
			cxxopts::value<std::string>()->default_value("config.toml"), "file")
		("n,grid", "Grid half-size (default 3).", cxxopts::value<int>()->default_value("3"), "n");

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
	const int n = result["n"].as<int>();

	std::string track_dir = brill::JoinPath(config.workspace, config.paths.track);
	std::string prefix = (trigger == "main") ? "" : "t1_";
	std::string track_path = TString::Format(
		"%s/ppac_%s%04d.root", track_dir.c_str(), prefix.c_str(), run
	).Data();

	std::cout << "Track file: " << track_path << "\n";
	std::cout << "Grid half-size n = " << n << " (coarse step=0.5mm, fine step=0.1mm)\n\n";

	brill::PpacPositionAdjust result_adj;
	if (brill::AdjustPpacPosition(track_path, config.ppac, n, result_adj)) {
		return 1;
	}

	std::cout << "=== X direction ===\n";
	std::cout << "  Baseline residual (dx2=0, dx3=0): " << result_adj.baseline_residual_x << "\n";
	if (result_adj.dsigma_x2 < 0) {
		std::cout << "  Best dx2: " << result_adj.dx2 << " mm (not well-determined)\n";
		std::cout << "  Best dx3: " << result_adj.dx3 << " mm (not well-determined)\n";
	} else {
		std::cout << "  Best dx2: " << result_adj.dx2 << " +/- " << result_adj.dsigma_x2 << " mm\n";
		std::cout << "  Best dx3: " << result_adj.dx3 << " +/- " << result_adj.dsigma_x3 << " mm\n";
	}
	std::cout << "  Best residual: " << result_adj.best_residual_x
		<< "  (improvement: " << result_adj.improvement_x << "%)\n\n";

	std::cout << "=== Y direction ===\n";
	std::cout << "  Baseline residual (dy2=0, dy3=0): " << result_adj.baseline_residual_y << "\n";
	if (result_adj.dsigma_y2 < 0) {
		std::cout << "  Best dy2: " << result_adj.dy2 << " mm (not well-determined)\n";
		std::cout << "  Best dy3: " << result_adj.dy3 << " mm (not well-determined)\n";
	} else {
		std::cout << "  Best dy2: " << result_adj.dy2 << " +/- " << result_adj.dsigma_y2 << " mm\n";
		std::cout << "  Best dy3: " << result_adj.dy3 << " +/- " << result_adj.dsigma_y3 << " mm\n";
	}
	std::cout << "  Best residual: " << result_adj.best_residual_y
		<< "  (improvement: " << result_adj.improvement_y << "%)\n\n";

	std::string output_path = TString::Format(
		"%s/ppac_offset_%s%04d.txt", track_dir.c_str(), prefix.c_str(), run
	).Data();
	std::ofstream ofs(output_path);
	if (ofs.is_open()) {
		ofs << "# PPAC position adjustment (PPAC1 as reference)\n";
		ofs << "# dx2 dx3 dy2 dy3\n";
		ofs << result_adj.dx2 << " " << result_adj.dx3 << " "
			<< result_adj.dy2 << " " << result_adj.dy3 << "\n";
		ofs.close();
		std::cout << "Saved to " << output_path << "\n";
	} else {
		std::cerr << "Error: Cannot write to " << output_path << "\n";
		return 1;
	}

	return 0;
}