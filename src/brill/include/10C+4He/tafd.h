#pragma once

#include <string>

#include "include/event/ingot/tafd_event.h"

namespace brill {

struct TAFDConfig;

constexpr int kMaxTafd = 6;
constexpr int kMaxTafdFrontStrips = 16;
constexpr int kMaxTafdBackStrips = 8;

struct TafdCalibration {
	double front_k[kMaxTafd][kMaxTafdFrontStrips] = {};
	double front_b[kMaxTafd][kMaxTafdFrontStrips] = {};
	double back_k[kMaxTafd][kMaxTafdBackStrips] = {};
	double back_b[kMaxTafd][kMaxTafdBackStrips] = {};
	bool front_valid[kMaxTafd][kMaxTafdFrontStrips] = {};
	bool back_valid[kMaxTafd][kMaxTafdBackStrips] = {};
};

struct TafdTarget {
	bool found = false;
	double energy = 0.0;
	double x = 0.0;
	double y = 0.0;
	double z = 0.0;
	int tafd_index = -1;
	int ring = -1;
	int pie = -1;
};

int ReadTafdCalibration(const std::string &path, TafdCalibration &calib);

double TAFDGetRingRadius(const TAFDConfig &config, int ring);

double TAFDGetPieAngle(const TAFDConfig &config, int pie);

void TAFDGetPosition(
	const TAFDConfig &config,
	int tafd_index,
	int ring,
	int pie,
	double &x,
	double &y,
	double &z
);

TafdTarget FindTafdTarget(
	const TafdEvent &event,
	const TafdCalibration &calib,
	const TAFDConfig &config
);

}  // namespace brill