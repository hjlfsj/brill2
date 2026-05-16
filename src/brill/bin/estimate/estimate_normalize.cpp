#include <filesystem>
#include <cstdio>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include <TFile.h>
#include <TH1F.h>
#include <TTree.h>
#include <TString.h>

#include "external/cxxopts.hpp"
#include "include/config.h"
#include "include/event/forge/dssd_event.h"
#include "include/t0/dssd.h"

std::string JoinPath(const std::string &left, const std::string &right) {
	if (left.empty()) return right;
	if (right.empty()) return left;
	if (left.back() == '/') return left + right;
	return left + "/" + right;
}

void PrintUsage(const cxxopts::Options &options) {
	std::cout << options.help() << "\n";
}

void RemoveAdjacentStrips(brill::DssdEvent &event) {
	std::vector<int> front_index, back_index;
	for (int i = 0; i < event.front_num; ++i) {
		bool found = false;
		for (int j = i+1; j < event.front_num; ++j) {
			if (abs(event.front_strip[i]-event.front_strip[j])==1) {
				found = true;
				break;
			}
		}
		if (!found) front_index.push_back(i);
	}
	for (int i = 0; i < event.back_num; ++i) {
		bool found = false;
		for (int j = i+1; j < event.back_num; ++j) {
			if (abs(event.back_strip[i]-event.back_strip[j])==1) {
				found = true;
				break;
			}
		}
		if (!found) back_index.push_back(i);
	}
	event.front_num = int(front_index.size());
	event.back_num = int(back_index.size());
	for (size_t i = 0; i < front_index.size(); ++i) {
		event.front_strip[i] = event.front_strip[front_index[i]];
		event.front_energy[i] = event.front_energy[front_index[i]];
		event.front_time[i] = event.front_time[front_index[i]];
	}
	for (size_t i = 0; i < back_index.size(); ++i) {
		event.back_strip[i] = event.back_strip[back_index[i]];
		event.back_energy[i] = event.back_energy[back_index[i]];
		event.back_time[i] = event.back_time[back_index[i]];
	}
}

int main(int argc, char **argv) {
	cxxopts::Options options("estimate_normalize", "Estimate normalized DSSD front-back result.");
	options.add_options()
		(
			"h,help",
			"Print help information."
		)
		(
			"c,config",
			"Config file path.",
			cxxopts::value<std::string>()->default_value("config.toml"),
			"file"
		)
		(
			"r,run",
			"Run number.",
			cxxopts::value<int>(),
			"run"
		)
		(
			"detector",
			"Detectors to estimate.",
			cxxopts::value<std::vector<std::string>>(),
			"detector"
		);
	options.parse_positional({"detector"});

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
	if (!result.count("detector")) {
		std::cerr << "Error: Missing detector positional arguments.\n";
		PrintUsage(options);
		return 1;
	}

	const std::set<std::string> allowed_detectors = {"t0d1", "t0d2", "t0d3", "t0d4"};
	std::vector<std::string> detectors = result["detector"].as<std::vector<std::string>>();
	for (const std::string &detector : detectors) {
		if (allowed_detectors.count(detector) == 0) {
			std::cerr << "Error: Unsupported detector " << detector << ".\n";
			return 1;
		}
	}

	brill::AppConfig config;
	if (brill::LoadConfig(result["config"].as<std::string>(), config)) {
		return 1;
	}
	const int run = result["run"].as<int>();

	for (const std::string &detector_name : detectors) {
		const brill::SquareDetectorConfig *detector =
			brill::FindDetectorConfig(config, detector_name);
		if (!detector) {
			std::cerr << "Error: Detector " << detector_name << " is not found in config.\n";
			return 1;
		}

		brill::DssdNormalizeParameters parameters;
		parameters.front_strips = detector->front_strips;
		parameters.back_strips = detector->back_strips;

		std::string normalize_dir = JoinPath(config.workspace, config.paths.normalize);
		std::string front_path = JoinPath(normalize_dir, detector_name + "_front.txt");
		std::string back_path = JoinPath(normalize_dir, detector_name + "_back.txt");
		if (brill::ReadDssdNormalizeParameters(front_path, back_path, parameters)) {
			return 1;
		}

		TString input_path = TString::Format(
			"%s/%s_%s%04d.root",
			JoinPath(config.workspace, config.paths.forge).c_str(),
			detector_name.c_str(),
			config.trigger == "main" ? "" : (config.trigger+"_").c_str(),
			run
		);
		TFile input_file(input_path, "read");
		TTree *input_tree = (TTree*)input_file.Get("tree");
		if (!input_tree) {
			std::cerr << "Error: Get tree from " << input_path << " failed.\n";
			return 1;
		}
		brill::DssdEvent raw_event;
		brill::SetupInput(input_tree, raw_event);

		TString output_path = TString::Format(
			"%s/%s_%s%04d.root",
			JoinPath(config.workspace, config.paths.estimate).c_str(),
			detector_name.c_str(),
			config.trigger == "main" ? "" : (config.trigger+"_").c_str(),
			run
		);
		TFile output_file(output_path, "recreate");
		TTree opt("tree", "normalized dssd");
		TH1F hist("hde", "normalized front-back energy difference", 1000, -5000, 5000);

		brill::DssdEvent normalized_event;
		for (long long entry = 0; entry < input_tree->GetEntriesFast(); ++entry) {
			input_tree->GetEntry(entry);
			brill::ApplyDssdNormalize(raw_event, parameters, normalized_event);
			RemoveAdjacentStrips(normalized_event);
			for (
				int i = 0;
				i < normalized_event.front_num && i < normalized_event.back_num;
				++i
			) {
				hist.Fill(normalized_event.front_energy[i] - normalized_event.back_energy[i]);
			}
		}

		hist.Write();
		output_file.Close();
		input_file.Close();
	}

	return 0;
}
