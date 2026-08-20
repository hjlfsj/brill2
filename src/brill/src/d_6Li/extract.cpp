#include "include/d_6Li/extract.h"

#include <TCutG.h>

#include <cmath>
#include <fstream>
#include <iostream>

namespace brill {

int ReadD6LiCalibration(const std::string &path, D6LiCalibration &calib) {
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
		if (index < 0 || index >= kD6LiCalibrationLayers) continue;
		calib.p0[index] = p0;
		calib.p1[index] = p1;
	}
	return 0;
}

bool PassD6LiCut(
	const DssdMatchEvent &d1,
	const DssdMatchEvent &d2,
	const DssdMatchEvent &d3,
	const DssdMatchEvent &d4,
	const D6LiCalibration &calib,
	TCutG *d3d4_cut
) {
	if (d1.num != 2) return false;
	if (d2.num != 2) return false;
	if (d3.num != 1) return false;
	if (d4.num != 1) return false;

	if (d3d4_cut) {
		double e3 = CalibrateD6LiEnergy(calib, 2, d3.energy[0]);
		double e4 = CalibrateD6LiEnergy(calib, 3, d4.energy[0]);
		if (!d3d4_cut->IsInside(e4, e3)) return false;
	}

	return true;
}

D6LiAdvancedResult ClassifyD6Li(
	const DssdMatchEvent &d1,
	const DssdMatchEvent &d2,
	const DssdMatchEvent &d3,
	const DssdMatchEvent &d4,
	const D6LiCalibration &calib,
	TCutG *d2d3_cut
) {
	D6LiAdvancedResult result;

	if (std::abs(d3.x[0] - d4.x[0]) >= 2.0 || std::abs(d3.y[0] - d4.y[0]) >= 2.0) return result;

	double e3_cal = CalibrateD6LiEnergy(calib, 2, d3.energy[0]);

	if (d2d3_cut) {
		int in_cut_count = 0;
		int in_cut_idx = -1;
		int out_cut_idx = -1;
		for (int i = 0; i < 2; ++i) {
			double e2_cal = CalibrateD6LiEnergy(calib, 1, d2.energy[i]);
			bool in_cut = d2d3_cut->IsInside(e3_cal, e2_cal);
			bool pos_match = std::abs(d2.x[i] - d3.x[0]) < 2.0 && std::abs(d2.y[i] - d3.y[0]) < 2.0;
			if (in_cut && pos_match) {
				in_cut_count++;
				in_cut_idx = i;
			} else {
				out_cut_idx = i;
			}
		}
		if (in_cut_count != 1) return result;
		result.e2_10C_idx = in_cut_idx;
		result.e2_6Li_idx = out_cut_idx;
	} else {
		return result;
	}

	double dx0 = d1.x[0] - d2.x[result.e2_10C_idx];
	double dy0 = d1.y[0] - d2.y[result.e2_10C_idx];
	double dist0 = dx0 * dx0 + dy0 * dy0;
	double dx1 = d1.x[1] - d2.x[result.e2_10C_idx];
	double dy1 = d1.y[1] - d2.y[result.e2_10C_idx];
	double dist1 = dx1 * dx1 + dy1 * dy1;

	if (dist0 < dist1) {
		result.e1_10C_idx = 0;
		result.e1_6Li_idx = 1;
	} else {
		result.e1_10C_idx = 1;
		result.e1_6Li_idx = 0;
	}

	result.passed = true;
	return result;
}

} // namespace brill