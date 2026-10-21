#include <cmath>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <string>

#include <TChain.h>
#include <TFile.h>
#include <TH2F.h>
#include <TString.h>

#include "external/cxxopts.hpp"
#include "include/config.h"
#include "include/event/ingot/tafd_event.h"
#include "include/event/ingot/csi_event.h"
#include "include/10C+4He/tafd.h"
#include "include/10C+4He/tafcsi.h"
#include "include/utils.h"

void PrintUsage(const cxxopts::Options &options) {
	std::cout << options.help() << "\n";
}

int main(int argc, char **argv) {
	cxxopts::Options options("check_taf_pid",
		"Check TAFD PID: CsI vs TAFD per channel (t1 trigger only).");
	options.add_options()
		("h,help", "Print help information.")
		("r,run", "Start run number.", cxxopts::value<int>(), "run")
		("e,end-run", "End run number (optional, for multi-run mode).",
			cxxopts::value<int>(), "run")
		("c,config", "Config file path.",
			cxxopts::value<std::string>()->default_value("config.toml"), "file");

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

	config.trigger = "t1";
	const int run = result["run"].as<int>();
	const int end_run = result.count("end-run") ? result["end-run"].as<int>() : run;
	if (end_run < run) {
		std::cerr << "Error: end run " << end_run
			<< " is smaller than run " << run << ".\n";
		return 1;
	}

	const std::string ingot_dir = brill::JoinPath(config.workspace, config.paths.ingot);
	const std::string trigger_infix = brill::TriggerInfix(config.trigger);

	TChain chain_tafd("tree");
	TChain chain_tafcsi("tree");
	int added_runs = 0;
	for (int current_run = run; current_run <= end_run; ++current_run) {
		if (brill::IsJumpRun(config, current_run)) continue;

		TString tafd_path = TString::Format(
			"%s/tafd_%s%04d.root",
			ingot_dir.c_str(), trigger_infix.c_str(), current_run);
		TString tafcsi_path = TString::Format(
			"%s/tafcsi_%s%04d.root",
			ingot_dir.c_str(), trigger_infix.c_str(), current_run);

		if (!std::filesystem::exists(tafd_path.Data())) {
			std::cerr << "Warning: TAFD file not found: "
				<< tafd_path << ", skipping run " << current_run << "\n";
			continue;
		}
		if (!std::filesystem::exists(tafcsi_path.Data())) {
			std::cerr << "Warning: TAFCSI file not found: "
				<< tafcsi_path << ", skipping run " << current_run << "\n";
			continue;
		}

		++added_runs;
		chain_tafd.Add(tafd_path);
		chain_tafcsi.Add(tafcsi_path);
	}
	if (added_runs == 0) {
		std::cerr << "No runs to process.\n";
		return 1;
	}
	chain_tafd.AddFriend(&chain_tafcsi, "csi");

	brill::TafdEvent tafd;
	brill::CsiEvent<12> csi;
	brill::SetupInput(&chain_tafd, tafd, "");
	brill::SetupInput(&chain_tafd, csi, "csi.");

	brill::TafdCalibration calib;
	std::string calib_path = "calib_vme_alpha_0255.txt";
	if (brill::ReadTafdCalibration(calib_path, calib) != 0) {
		std::cerr << "Error: Failed to load TAFD calibration from "
			<< calib_path << "\n";
		return 1;
	}
	std::cout << "TAFD calibration loaded from " << calib_path << "\n";

	TString output_path = TString::Format(
		"%s/check_taf_pid_%s%04d",
		brill::JoinPath(config.workspace, config.paths.estimate).c_str(),
		trigger_infix.c_str(), run);
	if (end_run != run) {
		output_path += TString::Format("_%04d", end_run);
	}
	output_path += ".root";
	TFile opf(output_path, "recreate");

	constexpr int kNumChannels = 12;
	TH2F* h_ch[kNumChannels];
	for (int i = 0; i < kNumChannels; ++i) {
		h_ch[i] = new TH2F(
			TString::Format("h_ch%d", i),
			TString::Format("CsI ch%d vs TAFD;CsI E (ch);TAFD calibrated E (MeV)", i),
			400, 0, 4000, 400, 0, 40);
	}
	TH2F* h_all = new TH2F("h_all",
		"CsI vs TAFD (all channels);CsI E (ch);TAFD calibrated E (MeV)",
		400, 0, 4000, 400, 0, 40);

	const long long total = chain_tafd.GetEntries();
	long long last_percentage = 0;
	int matched = 0;
	std::printf("Filling TAFD PID      0%%");
	std::fflush(stdout);

	for (long long entry = 0; entry < total; ++entry) {
		long long percentage = total > 0 ? entry * 100ll / total : 100ll;
		if (percentage > last_percentage) {
			last_percentage = percentage;
			std::printf("\b\b\b\b%3lld%%", percentage);
			std::fflush(stdout);
		}
		chain_tafd.GetEntry(entry);

		// Custom TAFD selection: no front/back consistency cut,
		// no calib validity check. Use ring energy for TAFD 0-4,
		// pie energy for TAFD5 (front damaged).
		double best_e = -1e9;
		int best_tafd = -1, best_ring = -1, best_pie = -1;
		for (int i = 0; i < 6; ++i) {
			if (!tafd.valid[i]) continue;
			int fs = tafd.front_strip[i];
			int bs = tafd.back_strip[i];
			if (fs < 0 || fs >= 16) continue;
			if (bs < 0 || bs >= 8) continue;

			double e;
			if (i == 5) {
				e = calib.back_k[i][bs] * tafd.back_energy[i]
					+ calib.back_b[i][bs];
			} else {
				e = calib.front_k[i][fs] * tafd.front_energy[i]
					+ calib.front_b[i][fs];
			}
			if (e > best_e) {
				best_e = e;
				best_tafd = i;
				best_ring = fs;
				best_pie = bs;
			}
		}
		if (best_tafd < 0) continue;

		auto csi_target = brill::FindTafcsiTarget(csi, best_tafd);
		if (!csi_target.found) continue;
		if (csi_target.channel < 0 || csi_target.channel >= kNumChannels) continue;

		++matched;
		h_ch[csi_target.channel]->Fill(csi_target.energy, best_e);
		h_all->Fill(csi_target.energy, best_e);
	}

	std::printf("\b\b\b\b100%%\n");
	std::printf("Total entries: %lld, matched: %d (%.1f%%)\n",
		total, matched, total > 0 ? matched * 100.0 / total : 0.0);

	opf.cd();
	for (int i = 0; i < kNumChannels; ++i) h_ch[i]->Write();
	h_all->Write();
	opf.Close();

	std::cout << "Output saved to " << output_path << "\n";
	return 0;
}