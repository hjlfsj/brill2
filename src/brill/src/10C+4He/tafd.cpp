#include "include/10C+4He/tafd.h"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include "include/config.h"

namespace brill {

namespace {

int ParseTafdIndex(const std::string &name) {
	if (name.size() < 5 || name.substr(0, 4) != "tafd") return -1;
	int index = std::atoi(name.substr(4).c_str());
	if (index < 0 || index >= kMaxTafd) return -1;
	return index;
}

}  // namespace

int ReadTafdCalibration(const std::string &path, TafdCalibration &calib) {
	std::ifstream fin(path);
	if (!fin.good()) {
		std::cerr << "Error: Open TAFD calibration file " << path << " failed.\n";
		return -1;
	}

	calib = TafdCalibration();

	std::string line;
	while (std::getline(fin, line)) {
		if (line.empty() || line[0] == '#') continue;

		std::istringstream iss(line);
		std::string detector, face, valid_str;
		int strip = -1;
		double k = 0.0, b = 0.0;

		iss >> detector >> face >> strip;
		for (int i = 0; i < 4; ++i) {
			std::string dummy;
			iss >> dummy;
		}
		iss >> k >> b >> valid_str;

		int idx = ParseTafdIndex(detector);
		if (idx < 0) continue;

		bool valid = (valid_str == "yes");
		if (strip < 0) continue;

		if (face == "front") {
			if (strip >= kMaxTafdFrontStrips) continue;
			calib.front_k[idx][strip] = k;
			calib.front_b[idx][strip] = b;
			calib.front_valid[idx][strip] = valid;
		} else if (face == "back") {
			if (strip >= kMaxTafdBackStrips) continue;
			calib.back_k[idx][strip] = k;
			calib.back_b[idx][strip] = b;
			calib.back_valid[idx][strip] = valid;
		}
	}

	return 0;
}

double TAFDGetRingRadius(const TAFDConfig &config, int ring) {
	if (ring < 0 || ring >= config.num_rings) {
		throw std::out_of_range("TAFD ring index out of range");
	}
	double pitch = (config.outer_radius_mm - config.inner_radius_mm) / config.num_rings;
	return config.inner_radius_mm + (ring + 0.5) * pitch;
}

double TAFDGetPieAngle(const TAFDConfig &config, int pie) {
	if (pie < 0 || pie >= config.num_pies) {
		throw std::out_of_range("TAFD pie index out of range");
	}
	double half_fan = config.fan_angle_deg * 0.5;
	double half_pie = config.fan_angle_deg / config.num_pies * 0.5;
	return half_fan - half_pie - pie * (config.fan_angle_deg / config.num_pies);
}

void TAFDGetPosition(
	const TAFDConfig &config,
	int tafd_index,
	int ring,
	int pie,
	double &x,
	double &y,
	double &z
) {
	if (tafd_index < 0 || tafd_index >= static_cast<int>(config.angles_deg.size())) {
		throw std::out_of_range("TAFD index out of range");
	}

	double deg2rad = M_PI / 180.0;

	double theta_center_deg = config.angles_deg[tafd_index];
	double theta_center = theta_center_deg * deg2rad;

	double r = TAFDGetRingRadius(config, ring);

	double phi_rel_deg = TAFDGetPieAngle(config, pie);
	double phi_global_deg = theta_center_deg + phi_rel_deg;
	double phi_global = phi_global_deg * deg2rad;

	double cx = config.center_distance_mm * std::cos(theta_center);
	double cy = config.center_distance_mm * std::sin(theta_center);

	x = cx + r * std::cos(phi_global);
	y = cy + r * std::sin(phi_global);
	z = config.z_mm;
}

TafdTarget FindTafdTarget(
	const TafdEvent &event,
	const TafdCalibration &calib,
	const TAFDConfig &config
) {
	TafdTarget target;

	for (int i = 0; i < kMaxTafd; ++i) {
		if (!event.valid[i]) continue;

		int fs = event.front_strip[i];
		int bs = event.back_strip[i];
		if (fs < 0 || fs >= kMaxTafdFrontStrips) continue;
		if (bs < 0 || bs >= kMaxTafdBackStrips) continue;

		if (!calib.front_valid[i][fs] || !calib.back_valid[i][bs]) continue;

		double fe_cal = calib.front_k[i][fs] * event.front_energy[i] + calib.front_b[i][fs];
		double be_cal = calib.back_k[i][bs] * event.back_energy[i] + calib.back_b[i][bs];

		if (std::fabs(fe_cal - be_cal) >= 0.2) continue;

		if (fe_cal > target.energy) {
			target.found = true;
			target.energy = fe_cal;
			target.tafd_index = i;
			target.ring = fs;
			target.pie = bs;
		}
	}

	if (target.found) {
		TAFDGetPosition(config, target.tafd_index, target.ring, target.pie,
			target.x, target.y, target.z);
	}

	return target;
}

}  // namespace brill