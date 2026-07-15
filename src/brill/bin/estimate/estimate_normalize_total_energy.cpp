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
#include "include/event/ingot/dssd_event.h"
#include "include/t0/dssd.h"
#include "include/utils.h"

void PrintUsage(const cxxopts::Options &options) {
	std::cout << options.help() << "\n";
}

int main(int argc, char **argv) {
	cxxopts::Options options("estimate_normalize_total_energy", "Estimate normalized DSSD total energy.");
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
			"t,trigger",
			"Trigger type.",
			cxxopts::value<std::string>(),
			"trigger"
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
	if (brill::IsJumpRun(config, run)) {
		std::cout << "Skipping jump run " << run << ".\n";
		return 0;
	}
	if (result.count("trigger")) {
		config.trigger = result["trigger"].as<std::string>();
	}
	int normalize_file_run = 0;
	for (const auto &[start, use] : config.normalize.runs) {
		if (run >= start) normalize_file_run = use;
	}

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

		std::string normalize_dir = brill::JoinPath(config.workspace, config.paths.normalize);
		TString front_path = TString::Format(
			"%s/%s_front_t1_%04d.txt",
			normalize_dir.c_str(),
			detector_name.c_str(),
			normalize_file_run
		);
		TString back_path = TString::Format(
			"%s/%s_back_t1_%04d.txt",
			normalize_dir.c_str(),
			detector_name.c_str(),
			normalize_file_run
		);
		if (brill::ReadDssdNormalizeParameters(
			front_path.Data(), back_path.Data(), parameters
		)) {
			return 1;
		}

		TString input_path = TString::Format(
			"%s/%s_%s%04d.root",
			brill::JoinPath(config.workspace, config.paths.ingot).c_str(),
			detector_name.c_str(),
			brill::TriggerInfix(config.trigger).c_str(),
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
			"%s/%s_normalize_total_energy_%s%04d.root",
			brill::JoinPath(config.workspace, config.paths.estimate).c_str(),
			detector_name.c_str(),
			brill::TriggerInfix(config.trigger).c_str(),
			run
		);

		std::cout << "Input file: " << input_path << std::endl;
		std::cout << "Output file: " << output_path << std::endl;

		TFile output_file(output_path, "recreate");
		TTree opt("tree", "normalized total energy");
		double fe_t = 0;
		double be_t = 0;
		int num_f = 0;
		int num_b = 0;
		opt.Branch("fe_t", &fe_t, "fe_t/D");
		opt.Branch("be_t", &be_t, "be_t/D");
		opt.Branch("num_f", &num_f, "num_f/I");
		opt.Branch("num_b", &num_b, "num_b/I");
		TH1F hist("hde", "normalized front-back energy difference", 1000, -5000, 5000);

		printf("Estimating %s   0%%", detector_name.c_str());
		fflush(stdout);
		long long total = input_tree->GetEntries();
		long long last_percentage = 0;
		for (long long entry = 0; entry < total; ++entry) {
			if (entry * 100ll / total > last_percentage) {
				last_percentage = entry * 100ll / total;
				printf("\b\b\b\b%3lld%%", last_percentage);
				fflush(stdout);
			}
			input_tree->GetEntry(entry);
			brill::DssdEvent normalized;
			brill::ApplyDssdNormalize(raw_event, parameters, normalized);

			fe_t = 0.0;
			be_t = 0.0;
			num_f = normalized.front_num;
			num_b = normalized.back_num;
			for (int i = 0; i < normalized.front_num; ++i) {
				fe_t += normalized.front_energy[i];
			}
			for (int i = 0; i < normalized.back_num; ++i) {
				be_t += normalized.back_energy[i];
			}
            hist.Fill(fe_t - be_t);
			opt.Fill();
		}
		printf("\b\b\b\b100%%\n");

		output_file.cd();
		hist.Write();
		opt.Write();
		output_file.Close();
		input_file.Close();
	}

	return 0;
}