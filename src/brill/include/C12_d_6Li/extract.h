#pragma once

#include "include/event/ingot/silicon_event.h"
#include "include/event/t0/dssd_match_event.h"

#include <string>

class TCutG;

namespace brill {

constexpr int kC12D6LiCalibrationLayers = 5;

struct C12D6LiCalibration {
	double p0[kC12D6LiCalibrationLayers] = {0.0};
	double p1[kC12D6LiCalibrationLayers] = {1.0};
};

int ReadC12D6LiCalibration(const std::string &path, C12D6LiCalibration &calib);

inline double CalibrateC12D6LiEnergy(
	const C12D6LiCalibration &calib,
	int layer,
	double raw_energy
) {
	return calib.p0[layer] + calib.p1[layer] * raw_energy;
}

bool PassC12D6LiCut(
	const DssdMatchEvent &d1,
	const DssdMatchEvent &d2,
	const DssdMatchEvent &d3,
	const DssdMatchEvent &d4,
	const C12D6LiCalibration &calib,
	TCutG *li6_cut = nullptr
);

struct Two4HeResult {
	bool passed = false;
	int idx_6Li_d1 = -1;
	int idx_6Li_d2 = -1;
	int e1_4He1_idx = -1;
	int e1_4He2_idx = -1;
	int e2_4He1_idx = -1;
	int e2_4He2_idx = -1;
	int e3_4He1_idx = -1;
	int e3_4He2_idx = -1;
	int e4_4He1_idx = -1;
	int e4_4He2_idx = -1;
};

Two4HeResult ClassifyTwo4He(
	const DssdMatchEvent &d1,
	const DssdMatchEvent &d2,
	const DssdMatchEvent &d3,
	const DssdMatchEvent &d4,
	const C12D6LiCalibration &calib,
	TCutG *li6_cut = nullptr
);

} // namespace brill