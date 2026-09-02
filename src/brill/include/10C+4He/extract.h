#pragma once

#include "include/event/ingot/silicon_event.h"
#include "include/event/t0/dssd_match_event.h"

#include <string>

class TCutG;

namespace brill {

constexpr int kC10He4CalibrationLayers = 5;

struct C10He4Calibration {
	double p0[kC10He4CalibrationLayers] = {0.0};
	double p1[kC10He4CalibrationLayers] = {1.0};
};

int ReadC10He4Calibration(const std::string &path, C10He4Calibration &calib);

inline double CalibrateC10He4Energy(
	const C10He4Calibration &calib,
	int layer,
	double raw_energy
) {
	return calib.p0[layer] + calib.p1[layer] * raw_energy;
}

bool Pass10C_d3_4He_s1Cut(
	const DssdMatchEvent &d1,
	const DssdMatchEvent &d2,
	const DssdMatchEvent &d3,
	const DssdMatchEvent &d4,
	const SiliconEvent &t0s,
	const C10He4Calibration &calib,
	TCutG *d2d3_cut = nullptr,
	int d1_hit = 2
);

} // namespace brill