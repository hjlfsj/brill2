#pragma once

#include <map>
#include <string>
#include <vector>

namespace brill {

constexpr int kMaxPpac = 3;
constexpr int kMaxStrips = 128;

struct SquareDetectorConfig {
	std::string name;
	int front_strips = 32;
	int back_strips = 32;
	double thickness_um = 0.0;
	double size_x_mm = 64.0;
	double size_y_mm = 64.0;
	double center_x_mm = 0.0;
	double center_y_mm = 0.0;
	double z_mm = 0.0;
	double merge_tolerance = 1000.0;
	double track_window_x = 4.0;
	double track_window_y = 4.0;
};

struct PpacConfig {
	double z_mm[kMaxPpac] = {0.0, 0.0, 0.0};
	double x_offset_mm[kMaxPpac] = {0.0, 0.0, 0.0};
	double y_offset_mm[kMaxPpac] = {0.0, 0.0, 0.0};
	double x_scale[kMaxPpac] = {1.0, 1.0, 1.0};
	double y_scale[kMaxPpac] = {1.0, 1.0, 1.0};
};

struct AppPaths {
	std::string decode = "decode";
	std::string forge = "forge";
	std::string fuse = "fuse";
	std::string normalize = "normalize";
	std::string match = "match";
	std::string track = "track";
	std::string estimate = "estimate";
};

struct AppConfig {
	std::string workspace = "/data/";
	std::string trigger = "";
	AppPaths paths;
	PpacConfig ppac;
	std::map<std::string, SquareDetectorConfig> detectors;
};

int LoadConfig(const std::string &path, AppConfig &config);

const SquareDetectorConfig *FindDetectorConfig(
	const AppConfig &config,
	const std::string &name
);

} // namespace brill
