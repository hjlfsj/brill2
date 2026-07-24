#include "include/analysis/ppac.h"

#include <cmath>
#include <iostream>
#include <vector>

#include <TFile.h>
#include <TH1F.h>
#include <TTree.h>

#include "include/analysis/file_name.h"
#include "include/event/analysis/ppac_track_event.h"

namespace glimmer {

namespace {

bool HasBit(int flag, int bit) {
	return (flag & (1 << bit)) != 0;
}

bool ReadNormalizeTree(const std::string &path, PpacNormalizeParameters &parameters) {
	TFile file(path.c_str(), "read");
	TTree *tree = dynamic_cast<TTree*>(file.Get("params"));
	if (!tree) return false;
	tree->SetBranchAddress("x_offset", parameters.x_offset);
	tree->SetBranchAddress("y_offset", parameters.y_offset);
	tree->SetBranchAddress("x_scale", parameters.x_scale);
	tree->SetBranchAddress("y_scale", parameters.y_scale);
	tree->GetEntry(0);
	return true;
}

} // namespace

bool ComputePpacRawPosition(const PpacEvent &event, int index, double &x, double &y) {
	int base = index * 4;
	if (!HasBit(event.flag, base) || !HasBit(event.flag, base + 1)) return false;
	if (!HasBit(event.flag, base + 2) || !HasBit(event.flag, base + 3)) return false;
	x = (event.time[base] - event.time[base + 1]) / 4.0;
	y = (event.time[base + 2] - event.time[base + 3]) / 4.0;
	return true;
}

void BuildPpacHitEvent(const PpacEvent &event, const PpacNormalizeParameters &parameters, PpacHitEvent &hit) {
	for (int i = 0; i < kMaxPpac; ++i) {
		double raw_x = 0.0;
		double raw_y = 0.0;
		hit.valid[i] = ComputePpacRawPosition(event, i, raw_x, raw_y) ? 1 : 0;
		hit.x[i] = raw_x * parameters.x_scale[i] + parameters.x_offset[i];
		hit.y[i] = raw_y * parameters.y_scale[i] + parameters.y_offset[i];
	}
}

bool FitPpacTrace(const double *z, const double *value, const bool *valid, double &intercept, double &slope) {
	double sum_z = 0.0;
	double sum_v = 0.0;
	double sum_zv = 0.0;
	double sum_z2 = 0.0;
	int count = 0;
	for (int i = 0; i < kMaxPpac; ++i) {
		if (!valid[i]) continue;
		sum_z += z[i];
		sum_v += value[i];
		sum_zv += z[i] * value[i];
		sum_z2 += z[i] * z[i];
		++count;
	}
	if (count < 2) return false;
	double denominator = count * sum_z2 - sum_z * sum_z;
	if (std::fabs(denominator) < 1e-9) return false;
	slope = (count * sum_zv - sum_z * sum_v) / denominator;
	intercept = (sum_v - slope * sum_z) / count;
	return true;
}

int RunPpacNormalize(const AppConfig &config, const std::string &trigger, int run, int end_run) {
	PpacNormalizeParameters parameters;
	TH1F hx[kMaxPpac] = {
		TH1F("hx0", "PPAC0 x offset", 400, -100, 100),
		TH1F("hx1", "PPAC1 x offset", 400, -100, 100),
		TH1F("hx2", "PPAC2 x offset", 400, -100, 100)
	};
	TH1F hy[kMaxPpac] = {
		TH1F("hy0", "PPAC0 y offset", 400, -100, 100),
		TH1F("hy1", "PPAC1 y offset", 400, -100, 100),
		TH1F("hy2", "PPAC2 y offset", 400, -100, 100)
	};
	double sum_x[kMaxPpac] = {0.0};
	double sum_y[kMaxPpac] = {0.0};
	int count_x[kMaxPpac] = {0};
	int count_y[kMaxPpac] = {0};

	for (int current = run; current <= end_run; ++current) {
		std::string ppac_path = ForgeFileName(config, "ppac", trigger, current);
		std::string vme_path = ForgeFileName(config, "vme_ppac", trigger, current);
		TFile ppac_file(ppac_path.c_str(), "read");
		TFile vme_file(vme_path.c_str(), "read");
		TTree *ppac_tree = dynamic_cast<TTree*>(ppac_file.Get("tree"));
		TTree *vme_tree = dynamic_cast<TTree*>(vme_file.Get("tree"));
		if (!ppac_tree || !vme_tree) {
			std::cerr << "Error: Open PPAC normalization inputs failed for run " << current << ".\n";
			return -1;
		}
		PpacEvent ppac;
		PpacEvent vme;
		SetupInput(ppac_tree, ppac, "");
		SetupInput(vme_tree, vme, "");
		long long total = std::min(ppac_tree->GetEntriesFast(), vme_tree->GetEntriesFast());
		for (long long entry = 0; entry < total; ++entry) {
			ppac_tree->GetEntry(entry);
			vme_tree->GetEntry(entry);
			for (int i = 0; i < kMaxPpac; ++i) {
				double px = 0.0;
				double py = 0.0;
				double vx = 0.0;
				double vy = 0.0;
				if (ComputePpacRawPosition(ppac, i, px, py) && ComputePpacRawPosition(vme, i, vx, vy)) {
					double dx = vx - px;
					double dy = vy - py;
					hx[i].Fill(dx);
					hy[i].Fill(dy);
					sum_x[i] += dx;
					sum_y[i] += dy;
					++count_x[i];
					++count_y[i];
				}
			}
		}
	}

	for (int i = 0; i < kMaxPpac; ++i) {
		parameters.x_offset[i] = count_x[i] > 0 ? sum_x[i] / count_x[i] : config.ppac.x_offset_mm[i];
		parameters.y_offset[i] = count_y[i] > 0 ? sum_y[i] / count_y[i] : config.ppac.y_offset_mm[i];
		parameters.x_scale[i] = config.ppac.x_scale[i];
		parameters.y_scale[i] = config.ppac.y_scale[i];
	}

	std::string output_path = NormalizeParamFileName(config, "ppac", trigger, run, end_run);
	TFile output(output_path.c_str(), "recreate");
	TTree params_tree("params", "ppac normalization parameters");
	params_tree.Branch("x_offset", parameters.x_offset, "x_offset[3]/D");
	params_tree.Branch("y_offset", parameters.y_offset, "y_offset[3]/D");
	params_tree.Branch("x_scale", parameters.x_scale, "x_scale[3]/D");
	params_tree.Branch("y_scale", parameters.y_scale, "y_scale[3]/D");
	params_tree.Fill();
	params_tree.Write();
	for (int i = 0; i < kMaxPpac; ++i) {
		hx[i].Write();
		hy[i].Write();
	}
	output.Close();
	return 0;
}

int ReadPpacNormalizeParameters(const std::string &path, PpacNormalizeParameters &parameters) {
	if (!ReadNormalizeTree(path, parameters)) {
		std::cerr << "Error: Read PPAC normalization parameters from " << path << " failed.\n";
		return -1;
	}
	return 0;
}

int RunPpacTrack(const AppConfig &config, const std::string &trigger, int run) {
	std::string ppac_path = ForgeFileName(config, "ppac", trigger, run);
	std::string param_path = NormalizeParamFileName(config, "ppac", trigger, run, run);
	TFile input(ppac_path.c_str(), "read");
	TTree *tree = dynamic_cast<TTree*>(input.Get("tree"));
	if (!tree) {
		std::cerr << "Error: Open PPAC forge file failed.\n";
		return -1;
	}
	PpacNormalizeParameters parameters;
	if (ReadPpacNormalizeParameters(param_path, parameters)) {
		for (int i = 0; i < kMaxPpac; ++i) {
			parameters.x_offset[i] = config.ppac.x_offset_mm[i];
			parameters.y_offset[i] = config.ppac.y_offset_mm[i];
			parameters.x_scale[i] = config.ppac.x_scale[i];
			parameters.y_scale[i] = config.ppac.y_scale[i];
		}
	}

	PpacEvent input_event;
	SetupInput(tree, input_event, "");
	std::string output_path = TrackFileName(config, "ppac-track", trigger, run);
	TFile output(output_path.c_str(), "recreate");
	TTree opt("tree", "ppac track");
	PpacTrackEvent track;
	SetupOutput(&opt, track);
	PpacHitEvent hit;
	double hz[kMaxPpac] = {
		config.ppac.z_mm[0],
		config.ppac.z_mm[1],
		config.ppac.z_mm[2]
	};
	for (long long entry = 0; entry < tree->GetEntriesFast(); ++entry) {
		tree->GetEntry(entry);
		BuildPpacHitEvent(input_event, parameters, hit);
		track.valid = 0;
		track.used_ppac = 0;
		bool valid[kMaxPpac] = {false, false, false};
		for (int i = 0; i < kMaxPpac; ++i) {
			track.ppac_valid[i] = hit.valid[i];
			track.ppac_x[i] = hit.x[i];
			track.ppac_y[i] = hit.y[i];
			valid[i] = hit.valid[i] != 0;
			if (valid[i]) track.used_ppac |= 1u << i;
		}
		double x0 = 0.0;
		double y0 = 0.0;
		double kx = 0.0;
		double ky = 0.0;
		if (FitPpacTrace(hz, hit.x, valid, x0, kx) && FitPpacTrace(hz, hit.y, valid, y0, ky)) {
			track.valid = 1;
			track.target_x = x0;
			track.target_y = y0;
			track.dir_x = kx;
			track.dir_y = ky;
		}
		opt.Fill();
	}
	opt.Write();
	output.Close();
	return 0;
}

} // namespace glimmer