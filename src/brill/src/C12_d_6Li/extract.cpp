#include "include/C12_d_6Li/extract.h"

#include <TCutG.h>

#include <cmath>
#include <fstream>
#include <iostream>

namespace brill {

int ReadC12D6LiCalibration(const std::string &path, C12D6LiCalibration &calib) {
	std::ifstream fin(path);
	if (!fin.good()) {
		std::cerr << "Error: Open calibration file " << path << " failed.\n";
		return -1;
	}

	std::string header;
	std::getline(fin, header);
	int index = -1;
	double p0 = 0.0;
	double p1 = 1.0;
	while (fin >> index >> p0 >> p1) {
		if (index < 0 || index >= kC12D6LiCalibrationLayers) continue;
		calib.p0[index] = p0;
		calib.p1[index] = p1;
	}
	return 0;
}

bool PassC12D6LiCut(
	const DssdMatchEvent &d1,
	const DssdMatchEvent &d2,
	const DssdMatchEvent &d3,
	const DssdMatchEvent &d4,
	const C12D6LiCalibration &calib,
	TCutG *li6_cut
) {
	if (d1.num != 3) return false;
	if (d2.num != 3) return false;
	if (d3.num != 2) return false;
	if (d4.num != 2) return false;

	int idx_6Li_d1 = 0;
	int idx_6Li_d2 = 0;
	for (int i = 1; i < d1.num; ++i) {
		if (d1.energy[i] > d1.energy[idx_6Li_d1]) idx_6Li_d1 = i;
	}
	for (int i = 1; i < d2.num; ++i) {
		if (d2.energy[i] > d2.energy[idx_6Li_d2]) idx_6Li_d2 = i;
	}

	if (li6_cut) {
		double e1_cal = CalibrateC12D6LiEnergy(calib, 0, d1.energy[idx_6Li_d1]);
		double e2_cal = CalibrateC12D6LiEnergy(calib, 1, d2.energy[idx_6Li_d2]);
		if (!li6_cut->IsInside(e2_cal, e1_cal)) return false;
	}

	return true;
}

namespace {

int FindMaxEnergyIndex(const DssdMatchEvent &event) {
	int idx = 0;
	for (int i = 1; i < event.num; ++i) {
		if (event.energy[i] > event.energy[idx]) idx = i;
	}
	return idx;
}

void GetRemainingIndices(const DssdMatchEvent &event, int exclude, int out[2]) {
	int cnt = 0;
	for (int i = 0; i < event.num; ++i) {
		if (i != exclude) out[cnt++] = i;
	}
}

double SqDist(const DssdMatchEvent &layer, int idx, const DssdMatchEvent &ref_layer, int ref_idx) {
	double dx = layer.x[idx] - ref_layer.x[ref_idx];
	double dy = layer.y[idx] - ref_layer.y[ref_idx];
	return dx * dx + dy * dy;
}

void AssignRemaining(
	const DssdMatchEvent &layer,
	const int layer_rem[2],
	const DssdMatchEvent &d2,
	const int d2_rem[2],
	int &out_4He1_idx,
	int &out_4He2_idx
) {
	double dist_a = SqDist(layer, layer_rem[0], d2, d2_rem[0])
	              + SqDist(layer, layer_rem[1], d2, d2_rem[1]);
	double dist_b = SqDist(layer, layer_rem[0], d2, d2_rem[1])
	              + SqDist(layer, layer_rem[1], d2, d2_rem[0]);

	if (dist_a <= dist_b) {
		out_4He1_idx = layer_rem[0];
		out_4He2_idx = layer_rem[1];
	} else {
		out_4He1_idx = layer_rem[1];
		out_4He2_idx = layer_rem[0];
	}
}

bool CheckAdjacent(
	const DssdMatchEvent &layer_a, int idx_a,
	const DssdMatchEvent &layer_b, int idx_b,
	double max_dx = 2.0, double max_dy = 2.0
) {
	double dx = std::abs(layer_a.x[idx_a] - layer_b.x[idx_b]);
	double dy = std::abs(layer_a.y[idx_a] - layer_b.y[idx_b]);
	return dx < max_dx && dy < max_dy;
}

} // namespace

Two4HeResult ClassifyTwo4He(
	const DssdMatchEvent &d1,
	const DssdMatchEvent &d2,
	const DssdMatchEvent &d3,
	const DssdMatchEvent &d4,
	const C12D6LiCalibration &calib,
	TCutG *li6_cut
) {
	Two4HeResult result;

	if (d1.num != 3 || d2.num != 3 || d3.num != 2 || d4.num != 2) return result;

	result.idx_6Li_d1 = FindMaxEnergyIndex(d1);
	result.idx_6Li_d2 = FindMaxEnergyIndex(d2);

	if (li6_cut) {
		double e1_cal = CalibrateC12D6LiEnergy(calib, 0, d1.energy[result.idx_6Li_d1]);
		double e2_cal = CalibrateC12D6LiEnergy(calib, 1, d2.energy[result.idx_6Li_d2]);
		if (!li6_cut->IsInside(e2_cal, e1_cal)) return result;
	}

	int d2_rem[2];
	GetRemainingIndices(d2, result.idx_6Li_d2, d2_rem);

	int d1_rem[2];
	GetRemainingIndices(d1, result.idx_6Li_d1, d1_rem);

	AssignRemaining(d1, d1_rem, d2, d2_rem, result.e1_4He1_idx, result.e1_4He2_idx);

	result.e2_4He1_idx = d2_rem[0];
	result.e2_4He2_idx = d2_rem[1];

	int d3_rem[2] = {0, 1};
	AssignRemaining(d3, d3_rem, d2, d2_rem, result.e3_4He1_idx, result.e3_4He2_idx);

	int d4_rem[2] = {0, 1};
	AssignRemaining(d4, d4_rem, d2, d2_rem, result.e4_4He1_idx, result.e4_4He2_idx);

	if (!CheckAdjacent(d1, result.e1_4He1_idx, d2, result.e2_4He1_idx)) return result;
	if (!CheckAdjacent(d1, result.e1_4He2_idx, d2, result.e2_4He2_idx)) return result;
	if (!CheckAdjacent(d2, result.e2_4He1_idx, d3, result.e3_4He1_idx)) return result;
	if (!CheckAdjacent(d2, result.e2_4He2_idx, d3, result.e3_4He2_idx)) return result;
	if (!CheckAdjacent(d3, result.e3_4He1_idx, d4, result.e4_4He1_idx)) return result;
	if (!CheckAdjacent(d3, result.e3_4He2_idx, d4, result.e4_4He2_idx)) return result;

	result.passed = true;
	return result;
}

} // namespace brill