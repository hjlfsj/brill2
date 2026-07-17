#pragma once

#include <map>
#include <string>
#include <vector>

namespace brill {

constexpr int kMaxPpac = 3;
constexpr int kMaxStrips = 128;

struct SquareDetectorConfig {
	std::string name;
	std::string type = "dssd";
	int front_strips = 32;
	int back_strips = 32;
	double thickness_um = 0.0;
	double size_x_mm = 64.0;
	double size_y_mm = 64.0;
	double center_x_mm = 0.0;
	double center_y_mm = 0.0;
	double z_mm = 0.0;
	double match_tolerance = 1000.0;
};

struct NromalizeStripsConfig {
	int norm_side = 0;
	int ref[2] = {0, 0};
	int norm[2] = {0, 0};
	double ref_energy[2] = {0.0, 60000.0};
	double norm_energy[2] = {0.0, 60000.0};
};

struct NormalizeDetectorConfig {
	std::vector<NromalizeStripsConfig> strips;
};

struct NormalizeConfig {
	std::vector<std::pair<int, int>> runs;
	std::map<std::string, NormalizeDetectorConfig> detectors;
};

struct TrackWindowConfig {
	double min = -4.0;
	double max = 4.0;
};

struct TrackConfig {
	TrackWindowConfig d2d1_window;
	TrackWindowConfig d3d2_window;
	TrackWindowConfig d4d3_window;
};

struct StraightParticleConfig {
	std::string particle;
	double mean = 0.0;
	double sigma = 0.0;
	double range = 0.0;
};

struct StraightSliceConfig {
	double A = 0.0;
	double B = 0.0;
	std::vector<StraightParticleConfig> particles;
};

struct IdentifyConfig {
	std::map<std::string, StraightSliceConfig> straight;
};

struct PpacConfig {
	double z_mm[kMaxPpac] = {0.0, 0.0, 0.0};
	double z_x_mm[kMaxPpac] = {0.0, 0.0, 0.0};
	double z_y_mm[kMaxPpac] = {0.0, 0.0, 0.0};
	double x_offset_mm[kMaxPpac] = {0.0, 0.0, 0.0};
	double y_offset_mm[kMaxPpac] = {0.0, 0.0, 0.0};
	double x_scale[kMaxPpac] = {1.0, 1.0, 1.0};
	double y_scale[kMaxPpac] = {1.0, 1.0, 1.0};
};

struct T0Config {
	std::vector<std::string> silicon;
};

struct AppPaths {
	std::string ore = "ore";
	std::string grit = "grit";
	std::string ingot = "ingot";
	std::string normalize = "normalize";
	std::string match = "match";
	std::string track = "track";
	std::string particle = "particle";
	std::string estimate = "estimate";
	std::string spectrum = "spectrum";
	std::string calibration = "calibration";
	std::string energy_calculator = "energy_calculator";
};

struct AppConfig {
	std::string workspace = "/data/";
	std::string trigger = "";
	std::string assets = "assets";
	std::vector<int> jump_run;
	AppPaths paths;
	T0Config t0;
	NormalizeConfig normalize;
	TrackConfig track;
	IdentifyConfig identify;
	PpacConfig ppac;
	std::map<std::string, SquareDetectorConfig> detectors;
};

int LoadConfig(const std::string &path, AppConfig &config);

const SquareDetectorConfig *FindDetectorConfig(
	const AppConfig &config,
	const std::string &name
);

const StraightSliceConfig *FindStraightSliceConfig(
	const AppConfig &config,
	const std::string &name
);

} // namespace brill