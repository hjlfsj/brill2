#pragma once

#include <map>
#include <string>
#include <vector>

namespace glimmer {

constexpr int kMaxPpac = 3;
constexpr int kMaxStrips = 128;
constexpr int kMaxDssdHit = 8;
constexpr int kMaxTrackHit = 8;
constexpr int kT0LayerCount = 4;

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
	double merge_tolerance = 3000.0;
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

struct AnalysisPaths {
	std::string forge = "forge";
	std::string normalize = "analysis/normalize";
	std::string merge = "analysis/merge";
	std::string track = "analysis/track";
	std::string show = "analysis/show";
};

struct AppConfig {
	std::string workspace = ".";
	std::string trigger_default;
	AnalysisPaths paths;
	PpacConfig ppac;
	std::map<std::string, SquareDetectorConfig> detectors;
};

int LoadConfig(const std::string &path, AppConfig &config);

const SquareDetectorConfig *FindDetectorConfig(
	const AppConfig &config,
	const std::string &name
);

} // namespace glimmer
