#include "include/analysis/t0.h"

#include <cmath>
#include <iostream>

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>

#include "include/analysis/dssd.h"
#include "include/analysis/file_name.h"
#include "include/event/analysis/ppac_track_event.h"
#include "include/event/analysis/t0_track_event.h"

namespace glimmer {

namespace {

void ResetTrackEvent(T0TrackEvent &event) {
	event.num = 0;
	for (int i = 0; i < kMaxTrackHit; ++i) {
		event.flag[i] = 0;
		for (int j = 0; j < kT0LayerCount; ++j) {
			event.energy[i][j] = 0.0;
			event.x[i][j] = 0.0;
			event.y[i][j] = 0.0;
			event.z[i][j] = 0.0;
		}
	}
}

void SeedTrackFromMerge(const DssdMergeEvent &merge, int layer, T0TrackEvent &track) {
	for (int i = 0; i < merge.num && track.num < kMaxTrackHit; ++i) {
		int index = track.num++;
		track.flag[index] = 1u << layer;
		track.energy[index][layer] = merge.energy[i];
		track.x[index][layer] = merge.x[i];
		track.y[index][layer] = merge.y[i];
		track.z[index][layer] = merge.z[i];
	}
}

void AttachLayer(
	const DssdMergeEvent &merge,
	int layer,
	double window_x,
	double window_y,
	T0TrackEvent &track
) {
	bool used[kMaxDssdHit] = {false};
	for (int i = 0; i < track.num; ++i) {
		int last_layer = -1;
		for (int j = layer - 1; j >= 0; --j) {
			if ((track.flag[i] & (1u << j)) != 0) {
				last_layer = j;
				break;
			}
		}
		if (last_layer < 0) continue;
		double best_metric = 1e9;
		int best = -1;
		for (int j = 0; j < merge.num; ++j) {
			if (used[j]) continue;
			double dx = merge.x[j] - track.x[i][last_layer];
			double dy = merge.y[j] - track.y[i][last_layer];
			if (std::fabs(dx) > window_x || std::fabs(dy) > window_y) continue;
			double metric = std::fabs(dx) + std::fabs(dy);
			if (metric < best_metric) {
				best_metric = metric;
				best = j;
			}
		}
		if (best >= 0) {
			used[best] = true;
			track.flag[i] |= 1u << layer;
			track.energy[i][layer] = merge.energy[best];
			track.x[i][layer] = merge.x[best];
			track.y[i][layer] = merge.y[best];
			track.z[i][layer] = merge.z[best];
		}
	}
}

void FillLayerCenter(const DssdMergeEvent &merge, TH1F &hx, TH1F &hy) {
	for (int i = 0; i < merge.num; ++i) {
		hx.Fill(merge.x[i]);
		hy.Fill(merge.y[i]);
	}
}

} // namespace

int RunT0Track(const AppConfig &config, const std::string &trigger, int run) {
	const char *names[kT0LayerCount] = {"t0d1", "t0d2", "t0d3", "t0d4"};
	const SquareDetectorConfig *detectors[kT0LayerCount];
	for (int i = 0; i < kT0LayerCount; ++i) {
		detectors[i] = FindDetectorConfig(config, names[i]);
		if (!detectors[i]) {
			std::cerr << "Error: Missing detector config for " << names[i] << ".\n";
			return -1;
		}
	}

	TFile files[kT0LayerCount] = {
		TFile(MergeFileName(config, names[0], trigger, run).c_str(), "read"),
		TFile(MergeFileName(config, names[1], trigger, run).c_str(), "read"),
		TFile(MergeFileName(config, names[2], trigger, run).c_str(), "read"),
		TFile(MergeFileName(config, names[3], trigger, run).c_str(), "read")
	};
	TTree *trees[kT0LayerCount];
	DssdMergeEvent merges[kT0LayerCount];
	for (int i = 0; i < kT0LayerCount; ++i) {
		trees[i] = dynamic_cast<TTree*>(files[i].Get("tree"));
		if (!trees[i]) {
			std::cerr << "Error: Open merge file for " << names[i] << " failed.\n";
			return -1;
		}
		SetupInput(trees[i], merges[i], "");
	}

	long long total = trees[0]->GetEntriesFast();
	for (int i = 1; i < kT0LayerCount; ++i) {
		total = std::min(total, trees[i]->GetEntriesFast());
	}

	std::string output_path = TrackFileName(config, "t0-track", trigger, run);
	TFile output(output_path.c_str(), "recreate");
	TTree opt("tree", "t0 track");
	T0TrackEvent track;
	SetupOutput(&opt, track);
	TH2F d1d2("d1d2", "D1-D2 PID", 1000, 0, 60000, 1000, 0, 60000);
	TH2F d2d3("d2d3", "D2-D3 PID", 1000, 0, 60000, 1000, 0, 60000);
	TH2F d3d4("d3d4", "D3-D4 PID", 1000, 0, 60000, 1000, 0, 60000);

	for (long long entry = 0; entry < total; ++entry) {
		for (int i = 0; i < kT0LayerCount; ++i) {
			trees[i]->GetEntry(entry);
		}
		ResetTrackEvent(track);
		SeedTrackFromMerge(merges[0], 0, track);
		AttachLayer(merges[1], 1, detectors[1]->track_window_x, detectors[1]->track_window_y, track);
		AttachLayer(merges[2], 2, detectors[2]->track_window_x, detectors[2]->track_window_y, track);
		AttachLayer(merges[3], 3, detectors[3]->track_window_x, detectors[3]->track_window_y, track);
		opt.Fill();
		FillDssdPidHistogram(merges[0], merges[1], d1d2);
		FillDssdPidHistogram(merges[1], merges[2], d2d3);
		FillDssdPidHistogram(merges[2], merges[3], d3d4);
	}
	opt.Write();
	d1d2.Write();
	d2d3.Write();
	d3d4.Write();
	output.Close();
	return 0;
}

int RunEstimatePpacCenter(const AppConfig &config, const std::string &trigger, int run) {
	std::string input_path = TrackFileName(config, "ppac-track", trigger, run);
	std::string output_path = ShowFileName(config, "estimate-ppac-center", trigger, run);
	TFile input(input_path.c_str(), "read");
	TTree *tree = dynamic_cast<TTree*>(input.Get("tree"));
	if (!tree) return -1;
	PpacTrackEvent event;
	SetupInput(tree, event, "");
	TFile output(output_path.c_str(), "recreate");
	TH1F hx("hx", "target x", 400, -40, 40);
	TH1F hy("hy", "target y", 400, -40, 40);
	for (long long entry = 0; entry < tree->GetEntriesFast(); ++entry) {
		tree->GetEntry(entry);
		if (!event.valid) continue;
		hx.Fill(event.target_x);
		hy.Fill(event.target_y);
	}
	hx.Write();
	hy.Write();
	output.Close();
	std::cout << "Estimated PPAC center x = " << hx.GetMean() << "\n";
	std::cout << "Estimated PPAC center y = " << hy.GetMean() << "\n";
	return 0;
}

int RunEstimateT0Center(const AppConfig &config, const std::string &trigger, int run) {
	const char *names[kT0LayerCount] = {"t0d1", "t0d2", "t0d3", "t0d4"};
	TFile files[kT0LayerCount] = {
		TFile(MergeFileName(config, names[0], trigger, run).c_str(), "read"),
		TFile(MergeFileName(config, names[1], trigger, run).c_str(), "read"),
		TFile(MergeFileName(config, names[2], trigger, run).c_str(), "read"),
		TFile(MergeFileName(config, names[3], trigger, run).c_str(), "read")
	};
	TTree *trees[kT0LayerCount];
	DssdMergeEvent merges[kT0LayerCount];
	for (int i = 0; i < kT0LayerCount; ++i) {
		trees[i] = dynamic_cast<TTree*>(files[i].Get("tree"));
		if (!trees[i]) return -1;
		SetupInput(trees[i], merges[i], "");
	}
	long long total = trees[0]->GetEntriesFast();
	for (int i = 1; i < kT0LayerCount; ++i) {
		total = std::min(total, trees[i]->GetEntriesFast());
	}
	std::string output_path = ShowFileName(config, "estimate-t0-center", trigger, run);
	TFile output(output_path.c_str(), "recreate");
	TH1F hx[kT0LayerCount] = {
		TH1F("d1x", "D1 x", 400, -40, 40),
		TH1F("d2x", "D2 x", 400, -40, 40),
		TH1F("d3x", "D3 x", 400, -40, 40),
		TH1F("d4x", "D4 x", 400, -40, 40)
	};
	TH1F hy[kT0LayerCount] = {
		TH1F("d1y", "D1 y", 400, -40, 40),
		TH1F("d2y", "D2 y", 400, -40, 40),
		TH1F("d3y", "D3 y", 400, -40, 40),
		TH1F("d4y", "D4 y", 400, -40, 40)
	};
	for (long long entry = 0; entry < total; ++entry) {
		for (int i = 0; i < kT0LayerCount; ++i) {
			trees[i]->GetEntry(entry);
			FillLayerCenter(merges[i], hx[i], hy[i]);
		}
	}
	for (int i = 0; i < kT0LayerCount; ++i) {
		hx[i].Write();
		hy[i].Write();
		std::cout << names[i] << " center x = " << hx[i].GetMean() << "\n";
		std::cout << names[i] << " center y = " << hy[i].GetMean() << "\n";
	}
	output.Close();
	return 0;
}

} // namespace glimmer
