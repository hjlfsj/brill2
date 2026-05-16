#include "include/config.h"

#include <iostream>

#include "external/toml.hpp"

namespace brill {

namespace {

void LoadPath(const toml::table &table, const char *key, std::string &value) {
	if (auto node = table[key].value<std::string>()) {
		value = *node;
	}
}

void LoadDetectorDouble(
	const toml::table &table,
	const char *key,
	double &value
) {
	if (auto node = table[key].value<double>()) {
		value = *node;
	}
}

void LoadDetectorInt(
	const toml::table &table,
	const char *key,
	int &value
) {
	if (auto node = table[key].value<int>()) {
		value = *node;
	}
}

void LoadPaths(const toml::table &table, AppPaths &paths) {
	LoadPath(table, "decode", paths.decode);
	LoadPath(table, "forge", paths.forge);
	LoadPath(table, "fuse", paths.fuse);
	LoadPath(table, "normalize", paths.normalize);
	LoadPath(table, "match", paths.match);
	LoadPath(table, "track", paths.track);
	LoadPath(table, "estimate", paths.estimate);
}

void LoadPpac(const toml::table &table, PpacConfig &ppac) {
	if (const auto *node = table["z_mm"].as_array()) {
		for (size_t i = 0; i < node->size() && i < kMaxPpac; ++i) {
			if (auto value = (*node)[i].value<double>()) {
				ppac.z_mm[i] = *value;
			}
		}
	}
	if (const auto *node = table["x_offset_mm"].as_array()) {
		for (size_t i = 0; i < node->size() && i < kMaxPpac; ++i) {
			if (auto value = (*node)[i].value<double>()) {
				ppac.x_offset_mm[i] = *value;
			}
		}
	}
	if (const auto *node = table["y_offset_mm"].as_array()) {
		for (size_t i = 0; i < node->size() && i < kMaxPpac; ++i) {
			if (auto value = (*node)[i].value<double>()) {
				ppac.y_offset_mm[i] = *value;
			}
		}
	}
	if (const auto *node = table["x_scale"].as_array()) {
		for (size_t i = 0; i < node->size() && i < kMaxPpac; ++i) {
			if (auto value = (*node)[i].value<double>()) {
				ppac.x_scale[i] = *value;
			}
		}
	}
	if (const auto *node = table["y_scale"].as_array()) {
		for (size_t i = 0; i < node->size() && i < kMaxPpac; ++i) {
			if (auto value = (*node)[i].value<double>()) {
				ppac.y_scale[i] = *value;
			}
		}
	}
}

void LoadDetector(const toml::table &table, const std::string &name, SquareDetectorConfig &detector) {
	detector.name = name;
	LoadDetectorInt(table, "front_strips", detector.front_strips);
	LoadDetectorInt(table, "back_strips", detector.back_strips);
	LoadDetectorDouble(table, "thickness_um", detector.thickness_um);
	LoadDetectorDouble(table, "size_x_mm", detector.size_x_mm);
	LoadDetectorDouble(table, "size_y_mm", detector.size_y_mm);
	LoadDetectorDouble(table, "center_x_mm", detector.center_x_mm);
	LoadDetectorDouble(table, "center_y_mm", detector.center_y_mm);
	LoadDetectorDouble(table, "z_mm", detector.z_mm);
	LoadDetectorDouble(table, "merge_tolerance", detector.merge_tolerance);
	LoadDetectorDouble(table, "match_tolerance", detector.merge_tolerance);
	LoadDetectorDouble(table, "track_window_x", detector.track_window_x);
	LoadDetectorDouble(table, "track_window_y", detector.track_window_y);
}

} // namespace

int LoadConfig(const std::string &path, AppConfig &config) {
	try {
		toml::table table = toml::parse_file(path);
		if (auto node = table["workspace"].value<std::string>()) {
			config.workspace = *node;
		} else {
			std::cerr << "Error: Missing required key workspace in " << path << ".\n";
			return -1;
		}
		if (auto node = table["trigger"].value<std::string>()) {
			config.trigger = *node;
		}

		if (const auto *paths = table["paths"].as_table()) {
			LoadPaths(*paths, config.paths);
		}

		if (const auto *detectors = table["detectors"].as_table()) {
			if (const auto *ppac = (*detectors)["ppac"].as_table()) {
				LoadPpac(*ppac, config.ppac);
			}
			for (const char *name : {"t0d1", "t0d2", "t0d3", "t0d4", "t0s"}) {
				if (const auto *detector = (*detectors)[name].as_table()) {
					SquareDetectorConfig config_detector;
					LoadDetector(*detector, name, config_detector);
					config.detectors[config_detector.name] = config_detector;
				}
			}
		}
	} catch (const toml::parse_error &err) {
		std::cerr << "Error: Parsing " << path << " failed:\n" << err << "\n";
		return -1;
	}
	return 0;
}

const SquareDetectorConfig *FindDetectorConfig(const AppConfig &config, const std::string &name) {
	auto iter = config.detectors.find(name);
	if (iter == config.detectors.end()) return nullptr;
	return &iter->second;
}

} // namespace brill
