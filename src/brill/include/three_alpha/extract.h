#pragma once

#include "include/event/t0/dssd_match_event.h"

#include <string>

class TCutG;

namespace brill {

constexpr int kThreeAlphaLayers = 4;

struct ThreeAlphaCalibration {
	double p0[kThreeAlphaLayers] = {0.0};
	double p1[kThreeAlphaLayers] = {1.0};
};

int ReadThreeAlphaCalibration(const std::string &path, ThreeAlphaCalibration &calib);

inline double CalibrateThreeAlphaEnergy(
	const ThreeAlphaCalibration &calib,
	int layer,
	double raw_energy
) {
	return calib.p0[layer] + calib.p1[layer] * raw_energy;
}

bool PassThreeAlphaCut(
	const DssdMatchEvent &d1,
	const DssdMatchEvent &d2,
	const DssdMatchEvent &d3,
	const DssdMatchEvent &d4
);

struct ThreeAlphaResult {
	bool passed = false;
	int He1_idx[4] = {-1, -1, -1, -1};
	int He2_idx[4] = {-1, -1, -1, -1};
	int He3_idx[4] = {-1, -1, -1, -1};
};

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
);

double FindIncidentEnergy(double e_total, double total_thickness_um,
	const std::string &cache_path);

} // namespace brill