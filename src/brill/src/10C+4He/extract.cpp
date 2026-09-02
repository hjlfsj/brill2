#include "include/10C+4He/extract.h"

#include <TCutG.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <utility>

namespace brill {

int ReadC10He4Calibration(const std::string &path, C10He4Calibration &calib) {
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
		if (index < 0 || index >= kC10He4CalibrationLayers) continue;
		calib.p0[index] = p0;
		calib.p1[index] = p1;
	}
	return 0;
}

namespace {

void SortHitsByEnergy(DssdMatchEvent &event) {
	if (event.num <= 1) return;

	int indices[8];
	for (int i = 0; i < event.num; ++i) indices[i] = i;

	std::sort(indices, indices + event.num, [&](int a, int b) {
		return event.energy[a] > event.energy[b];
	});

	double tmp_front_strip[8];
	double tmp_back_strip[8];
	double tmp_energy[8];
	double tmp_time[8];
	double tmp_x[8];
	double tmp_y[8];
	double tmp_z[8];
	int tmp_merge_tag[8];
	double tmp_energy_diff[8];

	for (int i = 0; i < event.num; ++i) {
		tmp_front_strip[i] = event.front_strip[indices[i]];
		tmp_back_strip[i] = event.back_strip[indices[i]];
		tmp_energy[i] = event.energy[indices[i]];
		tmp_time[i] = event.time[indices[i]];
		tmp_x[i] = event.x[indices[i]];
		tmp_y[i] = event.y[indices[i]];
		tmp_z[i] = event.z[indices[i]];
		tmp_merge_tag[i] = event.merge_tag[indices[i]];
		tmp_energy_diff[i] = event.energy_diff[indices[i]];
	}

	for (int i = 0; i < event.num; ++i) {
		event.front_strip[i] = tmp_front_strip[i];
		event.back_strip[i] = tmp_back_strip[i];
		event.energy[i] = tmp_energy[i];
		event.time[i] = tmp_time[i];
		event.x[i] = tmp_x[i];
		event.y[i] = tmp_y[i];
		event.z[i] = tmp_z[i];
		event.merge_tag[i] = tmp_merge_tag[i];
		event.energy_diff[i] = tmp_energy_diff[i];
	}
}

} // namespace

bool Pass10C_d3_4He_s1Cut(
	const DssdMatchEvent &d1,
	const DssdMatchEvent &d2,
	const DssdMatchEvent &d3,
	const DssdMatchEvent &d4,
	const SiliconEvent &t0s,
	const C10He4Calibration &calib,
	TCutG *d2d3_cut,
	int d1_hit
) {
	if (d1.num != d1_hit) return false;
	if (d2.num != 2) return false;
	if (d3.num != 2) return false;
	if (d4.num != 1) return false;

	DssdMatchEvent d1_sorted = d1;
	DssdMatchEvent d2_sorted = d2;
	DssdMatchEvent d3_sorted = d3;
	SortHitsByEnergy(d1_sorted);
	SortHitsByEnergy(d2_sorted);
	SortHitsByEnergy(d3_sorted);

	if (d2d3_cut) {
		double e3 = CalibrateC10He4Energy(calib, 2, d3_sorted.energy[0]);
		double e2 = CalibrateC10He4Energy(calib, 1, d2_sorted.energy[0]);
		if (!d2d3_cut->IsInside(e3, e2)) return false;
	}

	if (!t0s.valid) return false;

	return true;
}

} // namespace brill