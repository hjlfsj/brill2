#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include <TChain.h>
#include <TCutG.h>
#include <TFile.h>
#include <TString.h>

#include "external/cxxopts.hpp"
#include "include/config.h"
#include "include/event/ingot/silicon_event.h"
#include "include/event/t0/dssd_match_event.h"
#include "include/event/t0/t0_event.h"
#include "include/utils.h"

struct CutInfo {
	std::string particle;
	int charge = 0;
	int mass = 0;
	bool stop = false;
	std::unique_ptr<TCutG> cut;
};

int CutFileExists(
	const std::string &workspace,
	const std::string &slice,
	const std::string &particle,
	bool tail
) {
	std::ifstream fin(brill::CutFilePath(workspace, slice, particle, tail));
	return fin.good();
}

std::string ParticleElement(const std::string &particle) {
	size_t index = 0;
	while (
		index < particle.size()
		&& std::isdigit(static_cast<unsigned char>(particle[index]))
	) {
		++index;
	}
	return particle.substr(index);
}

int Loadcut(
	const std::string &workspace,
	const std::string &slice,
	const std::string &particle,
	const bool tail,
	std::unique_ptr<TCutG> &cut
) {
	if (!CutFileExists(workspace, slice, particle, true)) {
		if (!tail) {
			std::cerr << "Error: Cut file for " << slice << "_" << particle
				<< " not exist.\n";
			return -2;
		}
		std::string element = ParticleElement(particle);
		if (!CutFileExists(workspace, slice, element, true)) {
			std::cerr << "Error: Cut file for " << slice << "_" << particle
				<< "_tail not exist.\n";
			return -2;
		}
		return brill::ParseCutFile(workspace, slice, element, tail, cut);
	}
	return brill::ParseCutFile(workspace, slice, particle, tail, cut);
}

bool ParseParticleName(const std::string &particle, int &charge, int &mass) {
	size_t index = 0;
	while (index < particle.size() && std::isdigit(static_cast<unsigned char>(particle[index]))) {
		++index;
	}
	if (index == 0 || index >= particle.size()) return false;
	mass = std::stoi(particle.substr(0, index));
	charge = ElementCharge(particle.substr(index));
	return charge > 0;
}

int BuildCut(
	const std::string &workspace,
	const std::string &slice,
	const std::string &particle,
	const bool tail,
	CutInfo &cut
) {
	cut.particle = particle;
	if (!ParseParticleName(particle, cut.charge, cut.mass)) {
		std::cerr << "Error: Parse particle name failed.\n";
		return -1;
	}
	cut.stop = !tail;
	return Loadcut(workspace, slice, particle, tail, cut.cut);
}

int ElementCharge(const std::string &element) {
	if (element == "H") return 1;
	if (element == "He") return 2;
	if (element == "Li") return 3;
	if (element == "Be") return 4;
	if (element == "B") return 5;
	if (element == "C") return 6;
	if (element == "N") return 7;
	if (element == "O") return 8;
	return 0;
}

void PrintUsage(const cxxopts::Options &options) {
	std::cout << options.help() << "\n";
}

int main(int argc, char **argv) {
	cxxopts::Options options("track_t0", "Track and identify T0 particles.");
	options.add_options()
		("h,help", "Print help information.")
		("r,run", "Run number.", cxxopts::value<int>(), "run")
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
	if (brill::IsJumpRun(config, run)) {
		std::cout << "Skipping jump run " << run << ".\n";
		return 0;
	}

	CutInfo d2d3_7Be_stopped;
	CutInfo d2d3_7Be_passed;
	CutInfo d2d3_3He_stopped;
	CutInfo d2d3_3He_passed;
	CutInfo d3d4_7Be_stopped;
	CutInfo d3d4_7Be_passed;
	CutInfo d3d4_3He_stopped;
	CutInfo d3d4_3He_passed;
	if (BuildCut(config.workspace, "t0d2d3", "7Be", false, d2d3_7Be_stopped)) return -1;
	if (BuildCut(config.workspace, "t0d2d3", "7Be", true, d2d3_7Be_passed)) return -1;
	if (BuildCut(config.workspace, "t0d2d3", "3He", false, d2d3_3He_stopped)) return -1;
	if (BuildCut(config.workspace, "t0d2d3", "3He", true, d2d3_3He_passed)) return -1;
	if (BuildCut(config.workspace, "t0d3d4", "7Be", false, d3d4_7Be_stopped)) return -1;
	if (BuildCut(config.workspace, "t0d3d4", "7Be", true, d3d4_7Be_passed)) return -1;
	if (BuildCut(config.workspace, "t0d3d4", "3He", false, d3d4_3He_stopped)) return -1;
	if (BuildCut(config.workspace, "t0d3d4", "3He", true, d3d4_3He_passed)) return -1;

	const std::string match_dir = brill::JoinPath(config.workspace, config.paths.match);
	const std::string ingot_dir = brill::JoinPath(config.workspace, config.paths.ingot);
	const std::string trigger_infix = brill::TriggerInfix(config.trigger);

	TChain chain1("tree");
	chain1.Add(TString::Format(
		"%s/t0d1_%s%04d.root",
		match_dir.c_str(),
		trigger_infix.c_str(),
		run
	));
	TChain chain2("tree");
	chain2.Add(TString::Format(
		"%s/t0d2_%s%04d.root",
		match_dir.c_str(),
		trigger_infix.c_str(),
		run
	));
	TChain chain3("tree");
	chain3.Add(TString::Format(
		"%s/t0d3_%s%04d.root",
		match_dir.c_str(),
		trigger_infix.c_str(),
		run
	));
	TChain chain4("tree");
	chain4.Add(TString::Format(
		"%s/t0d4_%s%04d.root",
		match_dir.c_str(),
		trigger_infix.c_str(),
		run
	));
	TChain chain_s("tree");
	chain_s.Add(TString::Format(
		"%s/t0s_%s%04d.root",
		ingot_dir.c_str(),
		trigger_infix.c_str(),
		run
	));
	chain1.AddFriend(&chain2, "d2");
	chain1.AddFriend(&chain3, "d3");
	chain1.AddFriend(&chain4, "d4");
	chain1.AddFriend(&chain_s, "s");

	brill::DssdMatchEvent d1_event;
	brill::DssdMatchEvent d2_event;
	brill::DssdMatchEvent d3_event;
	brill::DssdMatchEvent d4_event;
	brill::SiliconEvent s_event;
	brill::SetupInput(&chain1, d1_event);
	brill::SetupInput(&chain1, d2_event, "d2.");
	brill::SetupInput(&chain1, d3_event, "d3.");
	brill::SetupInput(&chain1, d4_event, "d4.");
	brill::SetupInput(&chain1, s_event, "s.");

	TString output_path = TString::Format(
		"%s/t0_%s%04d.root",
		brill::JoinPath(config.workspace, config.paths.track).c_str(),
		trigger_infix.c_str(),
		run
	);
	TFile opf(output_path, "recreate");
	TTree opt("tree", "tracked t0");
	int layer[2];
	double layer_energy[2];
	double kinetic[2], energy[2];
	double px[2], py[2], pz[2];
	int run, entry;
	opt.Branch("layer", layer, "layer[2]/I");
	opt.Branch("layer")

	const long long total = chain1.GetEntries();
	long long last_percentage = -1;
	std::printf("Tracking T0   0%%");
	std::fflush(stdout);
	for (long long entry = 0; entry < total; ++entry) {
		long long percentage = total > 0 ? entry * 100ll / total : 100ll;
		if (percentage > last_percentage) {
			last_percentage = percentage;
			std::printf("\b\b\b\b%3lld%%", percentage);
			std::fflush(stdout);
		}
		chain1.GetEntry(entry);

		if (d2_event.num < 2 || d3_event.num < 2) continue;
		int be7_layer = 0;
		int he3_layer = 0;
		if (d2d3_7Be_passed.cut->IsInside(d3_event.energy[0], d2_event.energy[0])) {
			if (d2d3_3He_stopped.cut->IsInside(d3_event.energy[1], d2_event.energy[1])) {
				for (int i = 0; i < d4_event.num; ++i) {
					if (d3d4_7Be_stopped.cut->IsInside(d4_event.energy[i], d3_event.energy[0])) {
						be7_layer = 4;
						he3_layer = 3;
						break;
					}
				}
				if (!be7_layer) {
					for (int i = 0; i < d4_event.num; ++i) {
						if (d3d4_7Be_passed.cut->IsInside(d4_event.energy[i], d3_event.energy[0])) {
							be7_layer = 5;
							he3_layer = 3;
							break;
						}
					}
				}
			}
		}
		if (d2d3)

		opt.Fill();
	}
	std::printf("\b\b\b\b100%%\n");

	opf.cd();
	opt.Write();
	opf.Close();
	return 0;
	return 0;
}