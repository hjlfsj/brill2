#pragma once

#include "include/event/t0/dssd_match_event.h"

#include <string>

class TCutG;

namespace brill {

constexpr int kD6LiCalibrationLayers = 5;

struct D6LiCalibration {
	double p0[kD6LiCalibrationLayers] = {0.0};
	double p1[kD6LiCalibrationLayers] = {1.0};
};

int ReadD6LiCalibration(const std::string &path, D6LiCalibration &calib);

inline double CalibrateD6LiEnergy(
	const D6LiCalibration &calib,
	int layer,
	double raw_energy
) {
	return calib.p0[layer] + calib.p1[layer] * raw_energy;
}

bool PassD6LiCut(
	const DssdMatchEvent &d1,
	const DssdMatchEvent &d2,
	const DssdMatchEvent &d3,
	const DssdMatchEvent &d4,
	const D6LiCalibration &calib,
	TCutG *d3d4_cut = nullptr
);

struct D6LiAdvancedResult {
	bool passed = false;
	int e1_10C_idx = -1;
	int e1_6Li_idx = -1;
	int e2_10C_idx = -1;
	int e2_6Li_idx = -1;
};

D6LiAdvancedResult ClassifyD6Li(
	const DssdMatchEvent &d1,
	const DssdMatchEvent &d2,
	const DssdMatchEvent &d3,
	const DssdMatchEvent &d4,
	const D6LiCalibration &calib,
	TCutG *d2d3_cut = nullptr
);

} // namespace brill