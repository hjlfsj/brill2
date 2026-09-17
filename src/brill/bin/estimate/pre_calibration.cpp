#include <cmath>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <string>

#include <TChain.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH2F.h>
#include <TString.h>

#include "external/cxxopts.hpp"
#include "include/config.h"
#include "include/event/ingot/silicon_event.h"
#include "include/event/t0/dssd_match_event.h"
#include "include/utils.h"

namespace {

void PrintUsage(const cxxopts::Options &options) {
	std::cout << options.help() << "\n";
}

bool InTrackWindow(
	const brill::DssdMatchEvent &left,
	int left_index,
	const brill::DssdMatchEvent &right,
	int right_index,
	double max_dist_sq
) {
	double dx = right.x[right_index] - left.x[left_index];
	double dy = right.y[right_index] - left.y[right_index];
	return dx * dx + dy * dy <= max_dist_sq;
}

void FillPair(
	const brill::DssdMatchEvent &left,
	const brill::DssdMatchEvent &right,
	double max_dist_sq,
	TH2F &histogram,
	TGraph &graph
) {
	if (left.num < 1 || right.num < 1) return;
	if (!InTrackWindow(left, 0, right, 0, max_dist_sq)) return;
	histogram.Fill(right.energy[0], left.energy[0]);
	graph.SetPoint(graph.GetN(), right.energy[0], left.energy[0]);
}

void FillSilicon(
	const brill::DssdMatchEvent &d3,
	const brill::DssdMatchEvent &d4,
	const brill::SiliconEvent &silicon,
	double max_dist_sq,
	TH2F &histogram,
	TGraph &graph
) {
	if (!silicon.valid) return;
	if (d3.num != 1 || d4.num != 1) return;
	if (!InTrackWindow(d3, 0, d4, 0, max_dist_sq)) return;
	histogram.Fill(silicon.energy, d4.energy[0]);
	graph.SetPoint(graph.GetN(), silicon.energy, d4.energy[0]);
}

} // namespace

int main(int argc, char **argv) {
	cxxopts::Options options("pre_calibration", "Pre-calibration PID and TGraph generation.");
	options.add_options()
		("h,help", "Print help information.")
		("r,run", "Start run number.", cxxopts::value<int>(), "run")
		("e,end-run", "End run number.", cxxopts::value<int>(), "run")
		("t,trigger", "Trigger type.", cxxopts::value<std::string>(), "trigger")
		(
			"c,config",
			"Config file path.",
			cxxopts::value<std::string>()->default_value("config.toml"),
			"file"
		);

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
	if (result.count("trigger")) {
		config.trigger = result["trigger"].as<std::string>();
	}

	const int run = result["run"].as<int>();
	const int end_run = result.count("end-run") ? result["end-run"].as<int>() : run;
	if (end_run < run) {
		std::cerr << "Error: end run " << end_run << " is smaller than run " << run << ".\n";
		return 1;
	}

	const std::string match_dir = brill::JoinPath(config.workspace, config.paths.match);
	const std::string ingot_dir = brill::JoinPath(config.workspace, config.paths.ingot);
	const std::string trigger_infix = brill::TriggerInfix(config.trigger);

	TChain chain1("tree");
	TChain chain2("tree");
	TChain chain3("tree");
	TChain chain4("tree");
	TChain chain_s("tree");
	int added_runs = 0;
	for (int current_run = run; current_run <= end_run; ++current_run) {
		if (brill::IsJumpRun(config, current_run)) continue;
		++added_runs;
		chain1.Add(TString::Format(
			"%s/t0d1_%s%04d.root",
			match_dir.c_str(),
			trigger_infix.c_str(),
			current_run
		));
		chain2.Add(TString::Format(
			"%s/t0d2_%s%04d.root",
			match_dir.c_str(),
			trigger_infix.c_str(),
			current_run
		));
		chain3.Add(TString::Format(
			"%s/t0d3_%s%04d.root",
			match_dir.c_str(),
			trigger_infix.c_str(),
			current_run
		));
		chain4.Add(TString::Format(
			"%s/t0d4_%s%04d.root",
			match_dir.c_str(),
			trigger_infix.c_str(),
			current_run
		));
		chain_s.Add(TString::Format(
			"%s/t0s_%s%04d.root",
			ingot_dir.c_str(),
			trigger_infix.c_str(),
			current_run
		));
	}
	if (added_runs == 0) {
		std::cout << "No runs to process after applying jump_run.\n";
		return 0;
	}
	chain1.AddFriend(&chain2, "d2");
	chain1.AddFriend(&chain3, "d3");
	chain1.AddFriend(&chain4, "d4");
	chain1.AddFriend(&chain_s, "s");

	brill::DssdMatchEvent event1;
	brill::DssdMatchEvent event2;
	brill::DssdMatchEvent event3;
	brill::DssdMatchEvent event4;
	brill::SiliconEvent event_s;
	brill::SetupInput(&chain1, event1);
	brill::SetupInput(&chain1, event2, "d2.");
	brill::SetupInput(&chain1, event3, "d3.");
	brill::SetupInput(&chain1, event4, "d4.");
	brill::SetupInput(&chain1, event_s, "s.");

	TString output_path = TString::Format(
		"%s/pre_calibration_%s%04d_%04d.root",
		brill::JoinPath(config.workspace, config.paths.estimate).c_str(),
		trigger_infix.c_str(),
		run,
		end_run
	);
	TFile opf(output_path, "recreate");

	TH2F d1d2_pid("d1d2", "D1-D2 PID", 5000, 0.0, 60000.0, 5000, 0.0, 40000.0);
	TH2F d2d3_pid("d2d3", "D2-D3 PID", 5000, 0.0, 50000.0, 5000, 0.0, 70000.0);
	TH2F d3d4_pid("d3d4", "D3-D4 PID", 5000, 0.0, 50000.0, 5000, 0.0, 50000.0);
	TH2F d4s_pid("d4s", "D4-S PID", 5000, 0.0, 70000.0, 5000, 0.0, 50000.0);

	TGraph g_d1d2;
	TGraph g_d2d3;
	TGraph g_d3d4;
	TGraph g_d4s;
	g_d1d2.SetName("g_d1d2");
	g_d2d3.SetName("g_d2d3");
	g_d3d4.SetName("g_d3d4");
	g_d4s.SetName("g_d4s");

	const long long total = chain1.GetEntries();
	long long last_percentage = -1;
	std::printf("Pre-calibration   0%%");
	std::fflush(stdout);
	const double max_dist_sq = config.pre_calibration.max_distance_sq;

	for (long long entry = 0; entry < total; ++entry) {
		long long percentage = total > 0 ? entry * 100ll / total : 100ll;
		if (percentage > last_percentage) {
			last_percentage = percentage;
			std::printf("\b\b\b\b%3lld%%", percentage);
			std::fflush(stdout);
		}
		chain1.GetEntry(entry);
		FillPair(event1, event2, max_dist_sq, d1d2_pid, g_d1d2);
		FillPair(event2, event3, max_dist_sq, d2d3_pid, g_d2d3);
		FillPair(event3, event4, max_dist_sq, d3d4_pid, g_d3d4);
		FillSilicon(event3, event4, event_s, max_dist_sq, d4s_pid, g_d4s);
	}
	std::printf("\b\b\b\b100%%\n");

	opf.cd();
	d1d2_pid.Write();
	d2d3_pid.Write();
	d3d4_pid.Write();
	d4s_pid.Write();
	g_d1d2.Write();
	g_d2d3.Write();
	g_d3d4.Write();
	g_d4s.Write();
	opf.Close();

	return 0;
}