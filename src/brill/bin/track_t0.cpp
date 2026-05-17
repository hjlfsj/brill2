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
#include "include/event/forge/silicon_event.h"
#include "include/event/t0/dssd_match_event.h"
#include "include/event/t0/t0_event.h"
#include "include/utils.h"

namespace {

constexpr int kMaxHit = 8;
constexpr unsigned int kD2Offset = 0;
constexpr unsigned int kD3Offset = 8;
constexpr unsigned int kD4Offset = 16;
constexpr unsigned int kSOffset = 24;

struct TailCutInfo {
	std::string particle;
	int charge = 0;
	int mass = 0;
	std::unique_ptr<TCutG> cut;
};

struct Slice {
	int layer = -1;
	int charge = 0;
	int mass = 0;
	bool tail = false;
	int first = -1;
	int second = -1;
};

struct Chain {
	unsigned int mask = 0;
	int used_slices = 0;
	int charge = 0;
	int mass = 0;
	int layer = 0;
	int d1 = -1;
	int d2 = -1;
	int d3 = -1;
	int d4 = -1;
	bool has_s = false;
};

struct SelectionItem {
	unsigned int mask = 0;
	int used_slices = 0;
	double metric = 0.0;
	int payload = -1;
};

struct SelectionState {
	int prev = -1;
	int item = -1;
	int count = 0;
	int used_slices = 0;
	double metric = 0.0;
	unsigned int mask = 0;
};

struct D1D2Candidate {
	int d1 = -1;
	int d2 = -1;
	double metric = 0.0;
};

struct D1ExtensionCandidate {
	int particle = -1;
	int d1 = -1;
	double metric = 0.0;
};

void PrintUsage(const cxxopts::Options &options) {
	std::cout << options.help() << "\n";
}

bool InTrackWindow(
	double left_x,
	double left_y,
	double right_x,
	double right_y,
	const brill::TrackWindowConfig &window
) {
	double dx = right_x - left_x;
	double dy = right_y - left_y;
	return
		dx >= window.min
		&& dx <= window.max
		&& dy >= window.min
		&& dy <= window.max;
}

double TrackMetric(
	double left_x,
	double left_y,
	double right_x,
	double right_y
) {
	return std::fabs(right_x - left_x) + std::fabs(right_y - left_y);
}

double StraightEnergy(double de, double e, double a, double b) {
	return std::sqrt(de * e + a * de * de) + b * e;
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

std::string ParticleElement(const std::string &particle) {
	size_t index = 0;
	while (index < particle.size() && std::isdigit(static_cast<unsigned char>(particle[index]))) {
		++index;
	}
	return particle.substr(index);
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

bool CutFileExists(
	const std::string &workspace,
	const std::string &slice,
	const std::string &particle,
	bool tail
) {
	std::ifstream fin(brill::CutFilePath(workspace, slice, particle, tail));
	return fin.good();
}

int LoadTailCut(
	const std::string &workspace,
	const std::string &slice,
	const std::string &particle,
	std::unique_ptr<TCutG> &cut
) {
	if (CutFileExists(workspace, slice, particle, true)) {
		return brill::ParseCutFile(workspace, slice, particle, true, cut);
	}
	std::string element = ParticleElement(particle);
	if (
		!element.empty()
		&& element != particle
		&& CutFileExists(workspace, slice, element, true)
	) {
		return brill::ParseCutFile(workspace, slice, element, true, cut);
	}
	std::cerr
		<< "Error: Missing tail cut file for " << slice
		<< " particle " << particle << ".\n";
	return -1;
}

const brill::StraightParticleConfig *FindStraightParticle(
	const brill::StraightSliceConfig &slice,
	const std::string &particle
) {
	for (const auto &item : slice.particles) {
		if (item.particle == particle) return &item;
	}
	return nullptr;
}

void BuildStoppedDssdDssdSlices(
	const brill::DssdMatchEvent &left,
	const brill::DssdMatchEvent &right,
	const brill::TrackWindowConfig &window,
	const brill::StraightSliceConfig &straight,
	int layer,
	std::vector<Slice> &slices
) {
	for (int i = 0; i < left.num; ++i) {
		for (int j = 0; j < right.num; ++j) {
			if (!InTrackWindow(left.x[i], left.y[i], right.x[j], right.y[j], window)) continue;
			double ef = StraightEnergy(left.energy[i], right.energy[j], straight.A, straight.B);
			for (const auto &particle : straight.particles) {
				double low = particle.mean - particle.range * particle.sigma;
				double high = particle.mean + particle.range * particle.sigma;
				if (ef < low || ef > high) continue;
				int charge = 0;
				int mass = 0;
				if (!ParseParticleName(particle.particle, charge, mass)) continue;
				slices.push_back(Slice{layer, charge, mass, false, i, j});
			}
		}
	}
}

void BuildTailDssdDssdSlices(
	const brill::DssdMatchEvent &left,
	const brill::DssdMatchEvent &right,
	const brill::TrackWindowConfig &window,
	const std::vector<TailCutInfo> &cuts,
	int layer,
	std::vector<Slice> &slices
) {
	for (int i = 0; i < left.num; ++i) {
		for (int j = 0; j < right.num; ++j) {
			if (!InTrackWindow(left.x[i], left.y[i], right.x[j], right.y[j], window)) continue;
			for (const auto &cut : cuts) {
				if (!cut.cut || !cut.cut->IsInside(right.energy[j], left.energy[i])) continue;
				slices.push_back(Slice{layer, cut.charge, cut.mass, true, i, j});
			}
		}
	}
}

void BuildStoppedDssdSiliconSlices(
	const brill::DssdMatchEvent &dssd,
	const glimmer::SiliconEvent &silicon,
	const brill::StraightSliceConfig &straight,
	int layer,
	std::vector<Slice> &slices
) {
	if (!silicon.valid) return;
	for (int i = 0; i < dssd.num; ++i) {
		double ef = StraightEnergy(dssd.energy[i], silicon.energy, straight.A, straight.B);
		for (const auto &particle : straight.particles) {
			double low = particle.mean - particle.range * particle.sigma;
			double high = particle.mean + particle.range * particle.sigma;
			if (ef < low || ef > high) continue;
			int charge = 0;
			int mass = 0;
			if (!ParseParticleName(particle.particle, charge, mass)) continue;
			slices.push_back(Slice{layer, charge, mass, false, i, 0});
		}
	}
}

void BuildTailDssdSiliconSlices(
	const brill::DssdMatchEvent &dssd,
	const glimmer::SiliconEvent &silicon,
	const std::vector<TailCutInfo> &cuts,
	int layer,
	std::vector<Slice> &slices
) {
	if (!silicon.valid) return;
	for (int i = 0; i < dssd.num; ++i) {
		for (const auto &cut : cuts) {
			if (!cut.cut || !cut.cut->IsInside(silicon.energy, dssd.energy[i])) continue;
			slices.push_back(Slice{layer, cut.charge, cut.mass, true, i, 0});
		}
	}
}

void BuildChainsFromSlice(
	const std::vector<Slice> *slices,
	int layer,
	int index,
	Chain current,
	std::vector<Chain> &chains
) {
	const Slice &slice = slices[layer][index];
	current.used_slices += 1;
	current.charge = slice.charge;
	current.mass = slice.mass;
	if (layer == 0) {
		current.mask |= 1u << (kD2Offset + slice.first);
		current.mask |= 1u << (kD3Offset + slice.second);
		current.d2 = slice.first;
		current.d3 = slice.second;
	} else if (layer == 1) {
		current.mask |= 1u << (kD3Offset + slice.first);
		current.mask |= 1u << (kD4Offset + slice.second);
		current.d3 = slice.first;
		current.d4 = slice.second;
	} else if (layer == 2) {
		current.mask |= 1u << (kD4Offset + slice.first);
		current.mask |= 1u << kSOffset;
		current.d4 = slice.first;
		current.has_s = true;
	}

	bool extended = false;
	if (slice.tail && layer == 0) {
		for (size_t i = 0; i < slices[1].size(); ++i) {
			const Slice &next = slices[1][i];
			if (
				next.charge == slice.charge
				&& next.mass == slice.mass
				&& next.first == slice.second
			) {
				extended = true;
				BuildChainsFromSlice(slices, 1, int(i), current, chains);
			}
		}
	} else if (slice.tail && layer == 1) {
		for (size_t i = 0; i < slices[2].size(); ++i) {
			const Slice &next = slices[2][i];
			if (
				next.charge == slice.charge
				&& next.mass == slice.mass
				&& next.first == slice.second
			) {
				extended = true;
				BuildChainsFromSlice(slices, 2, int(i), current, chains);
			}
		}
	}

	if (extended) return;
	if (slice.tail && layer < 2) return;

	if (layer == 0) {
		current.layer = 3;
	} else if (layer == 1) {
		current.layer = 4;
	} else {
		current.layer = slice.tail ? 6 : 5;
	}
	chains.push_back(current);
}

void BuildChains(const std::vector<Slice> *slices, std::vector<Chain> &chains) {
	chains.clear();
	Chain empty;
	for (size_t i = 0; i < slices[0].size(); ++i) {
		BuildChainsFromSlice(slices, 0, int(i), empty, chains);
	}
	for (size_t i = 0; i < slices[1].size(); ++i) {
		BuildChainsFromSlice(slices, 1, int(i), empty, chains);
	}
	for (size_t i = 0; i < slices[2].size(); ++i) {
		BuildChainsFromSlice(slices, 2, int(i), empty, chains);
	}
}

void SelectBestSubset(
	const std::vector<SelectionItem> &items,
	bool prefer_slices,
	bool prefer_small_metric,
	std::vector<int> &selected
) {
	std::vector<SelectionState> groups;
	groups.push_back(SelectionState{0, -1, 0, 0, 0.0, 0});
	for (size_t i = 0; i < items.size(); ++i) {
		groups.push_back(SelectionState{
			0,
			int(i),
			1,
			items[i].used_slices,
			items[i].metric,
			items[i].mask
		});
	}

	size_t start = 1;
	size_t tail = groups.size();
	while (start < tail) {
		for (size_t i = start; i < tail; ++i) {
			for (size_t j = size_t(groups[i].item + 1); j < items.size(); ++j) {
				if ((groups[i].mask & items[j].mask) != 0) continue;
				groups.push_back(SelectionState{
					int(i),
					int(j),
					groups[i].count + 1,
					groups[i].used_slices + items[j].used_slices,
					groups[i].metric + items[j].metric,
					groups[i].mask | items[j].mask
				});
			}
		}
		start = tail;
		tail = groups.size();
	}

	int best = 0;
	for (size_t i = 1; i < groups.size(); ++i) {
		const SelectionState &candidate = groups[i];
		const SelectionState &current = groups[best];
		bool better = false;
		if (candidate.count > current.count) {
			better = true;
		} else if (candidate.count == current.count && prefer_slices && candidate.used_slices > current.used_slices) {
			better = true;
		} else if (
			candidate.count == current.count
			&& (!prefer_slices || candidate.used_slices == current.used_slices)
			&& prefer_small_metric
			&& candidate.metric < current.metric
		) {
			better = true;
		}
		if (better) best = int(i);
	}

	selected.clear();
	for (int index = best; index > 0; index = groups[index].prev) {
		selected.push_back(items[groups[index].item].payload);
	}
}

void BuildTailCuts(
	const std::string &workspace,
	const std::string &slice,
	const brill::StraightSliceConfig &straight,
	std::vector<TailCutInfo> &cuts
) {
	cuts.clear();
	std::set<std::string> loaded;
	for (const auto &particle : straight.particles) {
		if (loaded.count(particle.particle) != 0) continue;
		loaded.insert(particle.particle);
		TailCutInfo cut_info;
		cut_info.particle = particle.particle;
		if (!ParseParticleName(cut_info.particle, cut_info.charge, cut_info.mass)) continue;
		if (LoadTailCut(workspace, slice, cut_info.particle, cut_info.cut)) {
			throw std::runtime_error("load tail cut failed");
		}
		cuts.push_back(std::move(cut_info));
	}
}

void MarkUsed(const std::vector<Chain> &particles, bool *used_d2, bool *used_d3, bool *used_d4) {
	for (int i = 0; i < kMaxHit; ++i) {
		used_d2[i] = false;
		used_d3[i] = false;
		used_d4[i] = false;
	}
	for (const auto &particle : particles) {
		if (particle.d2 >= 0) used_d2[particle.d2] = true;
		if (particle.d3 >= 0) used_d3[particle.d3] = true;
		if (particle.d4 >= 0) used_d4[particle.d4] = true;
	}
}

int ParticleFlag(const Chain &particle) {
	int flag = 0;
	if (particle.d1 >= 0) flag |= 0x1;
	if (particle.d2 >= 0) flag |= 0x2;
	if (particle.d3 >= 0) flag |= 0x4;
	if (particle.d4 >= 0) flag |= 0x8;
	if (particle.has_s) flag |= 0x10;
	return flag;
}

void FillParticle(
	const Chain &particle,
	const brill::DssdMatchEvent &d1,
	const brill::DssdMatchEvent &d2,
	const brill::DssdMatchEvent &d3,
	const brill::DssdMatchEvent &d4,
	const glimmer::SiliconEvent &s,
	brill::T0Event &event
) {
	if (event.num >= 8) return;
	int index = event.num++;
	event.layer[index] = particle.layer;
	event.flag[index] = ParticleFlag(particle);
	event.charge[index] = particle.charge;
	event.mass[index] = particle.mass;

	if (particle.d1 >= 0) {
		event.energy[index][0] = d1.energy[particle.d1];
		event.time[index][0] = d1.time[particle.d1];
		event.x[index][0] = d1.x[particle.d1];
		event.y[index][0] = d1.y[particle.d1];
		event.z[index][0] = d1.z[particle.d1];
		event.last[index][0] = particle.d1;
	}
	if (particle.d2 >= 0) {
		event.energy[index][1] = d2.energy[particle.d2];
		event.time[index][1] = d2.time[particle.d2];
		event.x[index][1] = d2.x[particle.d2];
		event.y[index][1] = d2.y[particle.d2];
		event.z[index][1] = d2.z[particle.d2];
		event.last[index][1] = particle.d2;
	}
	if (particle.d3 >= 0) {
		event.energy[index][2] = d3.energy[particle.d3];
		event.time[index][2] = d3.time[particle.d3];
		event.x[index][2] = d3.x[particle.d3];
		event.y[index][2] = d3.y[particle.d3];
		event.z[index][2] = d3.z[particle.d3];
		event.last[index][2] = particle.d3;
	}
	if (particle.d4 >= 0) {
		event.energy[index][3] = d4.energy[particle.d4];
		event.time[index][3] = d4.time[particle.d4];
		event.x[index][3] = d4.x[particle.d4];
		event.y[index][3] = d4.y[particle.d4];
		event.z[index][3] = d4.z[particle.d4];
		event.last[index][3] = particle.d4;
	}
	if (particle.has_s && s.valid) {
		event.energy[index][4] = s.energy;
		event.time[index][4] = s.time;
	}
}

} // namespace

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

	const brill::StraightSliceConfig *d1d2_straight =
		brill::FindStraightSliceConfig(config, "t0d1d2");
	const brill::StraightSliceConfig *d2d3_straight =
		brill::FindStraightSliceConfig(config, "t0d2d3");
	const brill::StraightSliceConfig *d3d4_straight =
		brill::FindStraightSliceConfig(config, "t0d3d4");
	const brill::StraightSliceConfig *d4s_straight =
		brill::FindStraightSliceConfig(config, "t0d4s");
	if (!d1d2_straight || !d2d3_straight || !d3d4_straight || !d4s_straight) {
		std::cerr << "Error: Missing straight identify config.\n";
		return 1;
	}
	const brill::StraightParticleConfig *d1d2_he4 =
		FindStraightParticle(*d1d2_straight, "4He");
	if (!d1d2_he4) {
		std::cerr << "Error: Missing 4He straight identify window for t0d1d2.\n";
		return 1;
	}

	std::vector<TailCutInfo> d1d2_tail_cuts;
	std::vector<TailCutInfo> d2d3_tail_cuts;
	std::vector<TailCutInfo> d3d4_tail_cuts;
	std::vector<TailCutInfo> d4s_tail_cuts;
	try {
		BuildTailCuts(config.workspace, "t0d1d2", *d1d2_straight, d1d2_tail_cuts);
		BuildTailCuts(config.workspace, "t0d2d3", *d2d3_straight, d2d3_tail_cuts);
		BuildTailCuts(config.workspace, "t0d3d4", *d3d4_straight, d3d4_tail_cuts);
		BuildTailCuts(config.workspace, "t0d4s", *d4s_straight, d4s_tail_cuts);
	} catch (const std::runtime_error &) {
		return 1;
	}

	const std::string match_dir = brill::JoinPath(config.workspace, config.paths.match);
	const std::string forge_dir = brill::JoinPath(config.workspace, config.paths.forge);
	const std::string trigger_infix = brill::TriggerInfix(config.trigger);

	TString d1_path = TString::Format(
		"%s/t0d1_%s%04d.root",
		match_dir.c_str(),
		trigger_infix.c_str(),
		run
	);
	TFile ipf(d1_path, "read");
	TTree *ipt = (TTree*)ipf.Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from " << d1_path << " failed.\n";
		return 1;
	}
	ipt->AddFriend(
		"d2=tree",
		TString::Format("%s/t0d2_%s%04d.root", match_dir.c_str(), trigger_infix.c_str(), run)
	);
	ipt->AddFriend(
		"d3=tree",
		TString::Format("%s/t0d3_%s%04d.root", match_dir.c_str(), trigger_infix.c_str(), run)
	);
	ipt->AddFriend(
		"d4=tree",
		TString::Format("%s/t0d4_%s%04d.root", match_dir.c_str(), trigger_infix.c_str(), run)
	);
	ipt->AddFriend(
		"s=tree",
		TString::Format("%s/t0s_%s%04d.root", forge_dir.c_str(), trigger_infix.c_str(), run)
	);

	brill::DssdMatchEvent d1_event;
	brill::DssdMatchEvent d2_event;
	brill::DssdMatchEvent d3_event;
	brill::DssdMatchEvent d4_event;
	glimmer::SiliconEvent s_event;
	brill::SetupInput(ipt, d1_event);
	brill::SetupInput(ipt, d2_event, "d2.");
	brill::SetupInput(ipt, d3_event, "d3.");
	brill::SetupInput(ipt, d4_event, "d4.");
	glimmer::SetupInput(ipt, s_event, "s.");

	TString output_path = TString::Format(
		"%s/t0_%s%04d.root",
		brill::JoinPath(config.workspace, config.paths.track).c_str(),
		trigger_infix.c_str(),
		run
	);
	TFile opf(output_path, "recreate");
	TTree opt("tree", "tracked t0");
	brill::T0Event output_event;
	brill::Reset(output_event);
	brill::SetupOutput(&opt, output_event);

	const long long total = ipt->GetEntries();
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
		ipt->GetEntry(entry);
		brill::Reset(output_event);

		std::vector<Slice> slices[3];
		BuildStoppedDssdDssdSlices(
			d2_event, d3_event, config.track.d3d2_window, *d2d3_straight, 0, slices[0]
		);
		BuildTailDssdDssdSlices(
			d2_event, d3_event, config.track.d3d2_window, d2d3_tail_cuts, 0, slices[0]
		);
		BuildStoppedDssdDssdSlices(
			d3_event, d4_event, config.track.d4d3_window, *d3d4_straight, 1, slices[1]
		);
		BuildTailDssdDssdSlices(
			d3_event, d4_event, config.track.d4d3_window, d3d4_tail_cuts, 1, slices[1]
		);
		BuildStoppedDssdSiliconSlices(d4_event, s_event, *d4s_straight, 2, slices[2]);
		BuildTailDssdSiliconSlices(d4_event, s_event, d4s_tail_cuts, 2, slices[2]);

		std::vector<Chain> chains;
		BuildChains(slices, chains);

		std::vector<SelectionItem> chain_items;
		for (size_t i = 0; i < chains.size(); ++i) {
			chain_items.push_back(SelectionItem{
				chains[i].mask,
				chains[i].used_slices,
				0.0,
				int(i)
			});
		}
		std::vector<int> selected_chain_indices;
		SelectBestSubset(chain_items, true, false, selected_chain_indices);

		std::vector<Chain> particles;
		for (int index : selected_chain_indices) {
			particles.push_back(chains[index]);
		}

		bool used_d1[kMaxHit] = {false};
		bool used_d2[kMaxHit] = {false};
		bool used_d3[kMaxHit] = {false};
		bool used_d4[kMaxHit] = {false};
		MarkUsed(particles, used_d2, used_d3, used_d4);

		std::vector<SelectionItem> d1d2_items;
		std::vector<D1D2Candidate> d1d2_candidates;
		for (int i = 0; i < d1_event.num; ++i) {
			if (used_d1[i]) continue;
			for (int j = 0; j < d2_event.num; ++j) {
				if (used_d2[j]) continue;
				if (!InTrackWindow(
					d1_event.x[i], d1_event.y[i], d2_event.x[j], d2_event.y[j], config.track.d2d1_window
				)) {
					continue;
				}
				double ef = StraightEnergy(d1_event.energy[i], d2_event.energy[j], d1d2_straight->A, d1d2_straight->B);
				double low = d1d2_he4->mean - d1d2_he4->range * d1d2_he4->sigma;
				double high = d1d2_he4->mean + d1d2_he4->range * d1d2_he4->sigma;
				if (ef < low || ef > high) continue;
				double metric = TrackMetric(d1_event.x[i], d1_event.y[i], d2_event.x[j], d2_event.y[j]);
				d1d2_candidates.push_back(D1D2Candidate{i, j, metric});
				d1d2_items.push_back(SelectionItem{
					(1u << i) | (1u << (8 + j)),
					1,
					metric,
					int(d1d2_candidates.size() - 1)
				});
			}
		}
		std::vector<int> selected_d1d2;
		SelectBestSubset(d1d2_items, false, true, selected_d1d2);
		for (int candidate_index : selected_d1d2) {
			const D1D2Candidate &candidate = d1d2_candidates[candidate_index];
			Chain particle;
			particle.charge = 2;
			particle.mass = 4;
			particle.layer = 2;
			particle.d1 = candidate.d1;
			particle.d2 = candidate.d2;
			particles.push_back(particle);
			used_d1[candidate.d1] = true;
			used_d2[candidate.d2] = true;
		}

		std::unique_ptr<TCutG> he4_tail_cut;
		for (const auto &cut : d1d2_tail_cuts) {
			if (cut.particle == "4He") {
				he4_tail_cut.reset((TCutG*)cut.cut->Clone("t0d1d2_4He_tail_use"));
				break;
			}
		}
		if (!he4_tail_cut) {
			std::cerr << "Error: Missing usable t0d1d2 4He tail cut.\n";
			return 1;
		}

		std::vector<SelectionItem> extend_he4_items;
		std::vector<D1ExtensionCandidate> extend_he4_candidates;
		for (size_t i = 0; i < particles.size(); ++i) {
			if (particles[i].charge != 2 || particles[i].mass != 4) continue;
			if (particles[i].d2 < 0 || particles[i].d1 >= 0) continue;
			for (int j = 0; j < d1_event.num; ++j) {
				if (used_d1[j]) continue;
				if (!InTrackWindow(
					d1_event.x[j], d1_event.y[j],
					d2_event.x[particles[i].d2], d2_event.y[particles[i].d2],
					config.track.d2d1_window
				)) {
					continue;
				}
				if (!he4_tail_cut->IsInside(d2_event.energy[particles[i].d2], d1_event.energy[j])) continue;
				double metric = TrackMetric(
					d1_event.x[j], d1_event.y[j],
					d2_event.x[particles[i].d2], d2_event.y[particles[i].d2]
				);
				extend_he4_candidates.push_back(D1ExtensionCandidate{int(i), j, metric});
				extend_he4_items.push_back(SelectionItem{
					(1u << j) | (1u << (8 + int(i))),
					1,
					metric,
					int(extend_he4_candidates.size() - 1)
				});
			}
		}
		std::vector<int> selected_extend_he4;
		SelectBestSubset(extend_he4_items, false, true, selected_extend_he4);
		for (int candidate_index : selected_extend_he4) {
			const D1ExtensionCandidate &candidate = extend_he4_candidates[candidate_index];
			particles[candidate.particle].d1 = candidate.d1;
			used_d1[candidate.d1] = true;
		}

		std::vector<SelectionItem> extend_items;
		std::vector<D1ExtensionCandidate> extend_candidates;
		for (size_t i = 0; i < particles.size(); ++i) {
			if (particles[i].d2 < 0 || particles[i].d1 >= 0) continue;
			for (int j = 0; j < d1_event.num; ++j) {
				if (used_d1[j]) continue;
				if (!InTrackWindow(
					d1_event.x[j], d1_event.y[j],
					d2_event.x[particles[i].d2], d2_event.y[particles[i].d2],
					config.track.d2d1_window
				)) {
					continue;
				}
				double metric = TrackMetric(
					d1_event.x[j], d1_event.y[j],
					d2_event.x[particles[i].d2], d2_event.y[particles[i].d2]
				);
				extend_candidates.push_back(D1ExtensionCandidate{int(i), j, metric});
				extend_items.push_back(SelectionItem{
					(1u << j) | (1u << (8 + int(i))),
					1,
					metric,
					int(extend_candidates.size() - 1)
				});
			}
		}
		std::vector<int> selected_extensions;
		SelectBestSubset(extend_items, false, true, selected_extensions);
		for (int candidate_index : selected_extensions) {
			const D1ExtensionCandidate &candidate = extend_candidates[candidate_index];
			particles[candidate.particle].d1 = candidate.d1;
			used_d1[candidate.d1] = true;
		}

		for (const auto &particle : particles) {
			FillParticle(particle, d1_event, d2_event, d3_event, d4_event, s_event, output_event);
		}
		opt.Fill();
	}
	std::printf("\b\b\b\b100%%\n");

	opf.cd();
	opt.Write();
	opf.Close();
	ipf.Close();
	return 0;
}
