#include "include/three_alpha/extract.h"
#include "include/energy_calculator/lost_energy_calculator.h"
#include "include/energy_calculator/range_energy_calculator.h"

#include <TCutG.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>

namespace brill {

static catima::Material GetSiliconMaterial() {
	return SiliconMaterial();
}

int ReadThreeAlphaCalibration(const std::string &path, ThreeAlphaCalibration &calib) {
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
		if (index < 0 || index >= kThreeAlphaLayers) continue;
		calib.p0[index] = p0;
		calib.p1[index] = p1;
	}
	return 0;
}

bool PassThreeAlphaCut(
	const DssdMatchEvent &d1,
	const DssdMatchEvent &d2,
	const DssdMatchEvent &d3,
	const DssdMatchEvent &d4
) {
	if (d1.num != 3) return false;
	if (d2.num != 3) return false;
	if (d3.num != 3) return false;
	if (d4.num != 3) return false;
	return true;
}

namespace {

double SqDist(const DssdMatchEvent &a, int ia, const DssdMatchEvent &b, int ib) {
	double dx = a.x[ia] - b.x[ib];
	double dy = a.y[ia] - b.y[ib];
	return dx * dx + dy * dy;
}

void AssignThree(
	const DssdMatchEvent &layer,
	const int layer_indices[3],
	const DssdMatchEvent &d2,
	const int d2_indices[3],
	int out[3]
) {
	int perm[3] = {0, 1, 2};
	int best_perm[3] = {0, 1, 2};
	double best_dist = std::numeric_limits<double>::max();

	do {
		double dist = SqDist(layer, layer_indices[perm[0]], d2, d2_indices[0])
		            + SqDist(layer, layer_indices[perm[1]], d2, d2_indices[1])
		            + SqDist(layer, layer_indices[perm[2]], d2, d2_indices[2]);
		if (dist < best_dist) {
			best_dist = dist;
			best_perm[0] = perm[0];
			best_perm[1] = perm[1];
			best_perm[2] = perm[2];
		}
	} while (std::next_permutation(perm, perm + 3));

	for (int i = 0; i < 3; ++i) {
		out[i] = layer_indices[best_perm[i]];
	}
}

} // namespace

ThreeAlphaResult ClassifyThreeAlpha(
	const DssdMatchEvent &d1,
	const DssdMatchEvent &d2,
	const DssdMatchEvent &d3,
	const DssdMatchEvent &d4,
	const ThreeAlphaCalibration &calib,
	TCutG *d12_cut,
	TCutG *d23_cut,
	TCutG *d34_stop_cut,
	TCutG *d34_tran_cut
) {
	ThreeAlphaResult result;

	if (d1.num != 3 || d2.num != 3 || d3.num != 3 || d4.num != 3) return result;

	std::vector<int> d2_order = {0, 1, 2};
	std::sort(d2_order.begin(), d2_order.end(),
		[&](int a, int b) { return d2.energy[a] > d2.energy[b]; });

	int d2_He[3] = {d2_order[0], d2_order[1], d2_order[2]};

	int d1_all[3] = {0, 1, 2};
	int d1_out[3] = {-1, -1, -1};
	AssignThree(d1, d1_all, d2, d2_He, d1_out);

	int d3_all[3] = {0, 1, 2};
	int d3_out[3] = {-1, -1, -1};
	AssignThree(d3, d3_all, d2, d2_He, d3_out);

	int d4_all[3] = {0, 1, 2};
	int d4_out[3] = {-1, -1, -1};
	AssignThree(d4, d4_all, d2, d2_He, d4_out);

	for (int i = 0; i < 3; ++i) {
		double e1_cal = CalibrateThreeAlphaEnergy(calib, 0, d1.energy[d1_out[i]]);
		double e2_cal = CalibrateThreeAlphaEnergy(calib, 1, d2.energy[d2_He[i]]);
		double e3_cal = CalibrateThreeAlphaEnergy(calib, 2, d3.energy[d3_out[i]]);
		double e4_cal = CalibrateThreeAlphaEnergy(calib, 3, d4.energy[d4_out[i]]);

		if (!d12_cut->IsInside(e2_cal, e1_cal)) return result;
		if (!d23_cut->IsInside(e3_cal, e2_cal)) return result;

		bool in_stop = d34_stop_cut->IsInside(e4_cal, e3_cal);
		bool in_tran = d34_tran_cut->IsInside(e4_cal, e3_cal);
		if (!in_stop && !in_tran) return result;
	}

	result.He1_idx[0] = d1_out[0];
	result.He1_idx[1] = d2_He[0];
	result.He1_idx[2] = d3_out[0];
	result.He1_idx[3] = d4_out[0];

	result.He2_idx[0] = d1_out[1];
	result.He2_idx[1] = d2_He[1];
	result.He2_idx[2] = d3_out[1];
	result.He2_idx[3] = d4_out[1];

	result.He3_idx[0] = d1_out[2];
	result.He3_idx[1] = d2_He[2];
	result.He3_idx[2] = d3_out[2];
	result.He3_idx[3] = d4_out[2];

	result.passed = true;
	return result;
}

double FindIncidentEnergy(double e_total, double total_thickness_um,
	const std::string &cache_path) {
	catima::Material si_mat = GetSiliconMaterial();
	LostEnergyCalculator lec(2, 4, si_mat, total_thickness_um, cache_path);

	double lo = e_total;
	double hi = e_total * 5.0;

	while (lec.EnergyLost(hi) >= e_total) {
		hi *= 2.0;
	}

	for (int iter = 0; iter < 80; ++iter) {
		double mid = (lo + hi) * 0.5;
		double elost = lec.EnergyLost(mid);

		if (elost > e_total) {
			lo = mid;
		} else {
			hi = mid;
		}
	}

	return (lo + hi) * 0.5;
}

} // namespace brill