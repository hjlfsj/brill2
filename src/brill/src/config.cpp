#include "include/config.h"

#include <iostream>
#include <algorithm>

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
	LoadPath(table, "ore", paths.ore);
	LoadPath(table, "grit", paths.grit);
	LoadPath(table, "ingot", paths.ingot);
	LoadPath(table, "normalize", paths.normalize);
	LoadPath(table, "match", paths.match);
	LoadPath(table, "track", paths.track);
	LoadPath(table, "particle", paths.particle);
	LoadPath(table, "estimate", paths.estimate);
	LoadPath(table, "spectrum", paths.spectrum);
	LoadPath(table, "calibration", paths.calibration);
	LoadPath(table, "energy_calculator", paths.energy_calculator);
}

template<typename T>
void LoadNormalizeStripWindow(
	const toml::table &table,
	const char *key,
	T *fill
) {
	const auto *range = table[key].as_array();
	if (!range || range->size() != 2) return;
	std::optional<T> range_start = (*range)[0].template value<T>();
	std::optional<T> range_end = (*range)[1].template value<T>();
	if (range_start) fill[0] = *range_start;
	if (range_end) fill[1] = *range_end;
}

void LoadNormalizeStrips(
	const toml::table &table,
	NromalizeStripsConfig &strip,
	int &index
) {
	LoadDetectorInt(table, "index", index);
	LoadDetectorInt(table, "norm_side", strip.norm_side);
	LoadNormalizeStripWindow(table, "ref", strip.ref);
	LoadNormalizeStripWindow(table, "norm", strip.norm);
	LoadNormalizeStripWindow(table, "ref_energy", strip.ref_energy);
	LoadNormalizeStripWindow(table, "norm_energy", strip.norm_energy);
}

void LoadNormalizeDetector(
	const toml::table &table,
	NormalizeDetectorConfig &detector
) {
	if (const auto *strips = table["strips"].as_array()) {
		detector.strips.resize(strips->size());
		for (size_t i = 0; i < strips->size(); ++i) {
			const auto *strips_table = (*strips)[i].as_table();
			if (!strips_table) continue;
			NromalizeStripsConfig strip_config;
			int index = -1;
			LoadNormalizeStrips(*strips_table, strip_config, index);
			detector.strips[index] = strip_config;
		}
	}
}

void LoadNormalize(const toml::table &table, NormalizeConfig &normalize) {
	if (const auto *run_array = table["run"].as_array()) {
		for (const auto &run_item : *run_array) {
			const auto &run_table = *run_item.as_table();
			int start_run = 0;
			int use_run = 0;
			LoadDetectorInt(run_table, "run", start_run);
			LoadDetectorInt(run_table, "use", use_run);
			normalize.runs.push_back(std::make_pair(start_run, use_run));
		}
	}
	for (const auto &item : table) {
		if (!item.second.is_table()) continue;
		NormalizeDetectorConfig detector;
		LoadNormalizeDetector(*item.second.as_table(), detector);
		normalize.detectors[std::string(item.first.str())] = detector;
	}

	std::sort(
		normalize.runs.begin(),
		normalize.runs.end(),
		[](const std::pair<int, int>& a, const std::pair<int, int> &b) {
			return a.first < b.first;
		}
	);
}
void LoadTrackWindow(
	const toml::table &table,
	const char *key,
	TrackWindowConfig &window
) {
	const auto *node = table[key].as_array();
	if (!node || node->size() < 2) return;
	auto min = (*node)[0].value<double>();
	auto max = (*node)[1].value<double>();
	if (min) window.min = *min;
	if (max) window.max = *max;
}

void LoadTrack(const toml::table &table, TrackConfig &track) {
	LoadTrackWindow(table, "d2d1_window", track.d2d1_window);
	LoadTrackWindow(table, "d3d2_window", track.d3d2_window);
	LoadTrackWindow(table, "d4d3_window", track.d4d3_window);
}

void LoadT0(const toml::table &table, T0Config &t0) {
	if (const auto *silicon = table["silicon"].as_array()) {
		t0.silicon.clear();
		for (size_t i = 0; i < silicon->size(); ++i) {
			if (auto value = (*silicon)[i].value<std::string>()) {
				t0.silicon.push_back(*value);
			}
		}
	}
}

void LoadStraightParticle(
	const toml::table &table,
	StraightParticleConfig &particle
) {
	if (auto node = table["particle"].value<std::string>()) {
		particle.particle = *node;
	}
	LoadDetectorDouble(table, "mean", particle.mean);
	LoadDetectorDouble(table, "sigma", particle.sigma);
	LoadDetectorDouble(table, "range", particle.range);
}

void LoadStraightSlice(
	const toml::table &table,
	StraightSliceConfig &slice
) {
	LoadDetectorDouble(table, "A", slice.A);
	LoadDetectorDouble(table, "B", slice.B);
	if (const auto *particles = table["particles"].as_array()) {
		for (size_t i = 0; i < particles->size(); ++i) {
			const auto *particle_table = (*particles)[i].as_table();
			if (!particle_table) continue;
			StraightParticleConfig particle;
			LoadStraightParticle(*particle_table, particle);
			if (!particle.particle.empty()) {
				slice.particles.push_back(particle);
			}
		}
	}
}

void LoadIdentify(const toml::table &table, IdentifyConfig &identify) {
	const auto *straight = table["straight"].as_table();
	if (!straight) return;
	for (const auto &item : *straight) {
		if (!item.second.is_table()) continue;
		StraightSliceConfig slice;
		LoadStraightSlice(*item.second.as_table(), slice);
		identify.straight[std::string(item.first.str())] = slice;
	}
}

void LoadPpac(const toml::table &table, PpacConfig &ppac) {
	if (const auto *node = table["z_mm"].as_array()) {
		for (size_t i = 0; i < node->size() && i < kMaxPpac; ++i) {
			if (auto value = (*node)[i].value<double>()) {
				ppac.z_mm[i] = *value;
			}
		}
	}
	if (const auto *node = table["z_x_mm"].as_array()) {
		for (size_t i = 0; i < node->size() && i < kMaxPpac; ++i) {
			if (auto value = (*node)[i].value<double>()) {
				ppac.z_x_mm[i] = *value;
			}
		}
	} else {
		for (size_t i = 0; i < kMaxPpac; ++i) ppac.z_x_mm[i] = ppac.z_mm[i];
	}
	if (const auto *node = table["z_y_mm"].as_array()) {
		for (size_t i = 0; i < node->size() && i < kMaxPpac; ++i) {
			if (auto value = (*node)[i].value<double>()) {
				ppac.z_y_mm[i] = *value;
			}
		}
	} else {
		for (size_t i = 0; i < kMaxPpac; ++i) ppac.z_y_mm[i] = ppac.z_mm[i];
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
	if (auto node = table["type"].value<std::string>()) {
		detector.type = *node;
	}
	LoadDetectorInt(table, "front_strips", detector.front_strips);
	LoadDetectorInt(table, "back_strips", detector.back_strips);
	LoadDetectorDouble(table, "thickness_um", detector.thickness_um);
	LoadDetectorDouble(table, "size_x_mm", detector.size_x_mm);
	LoadDetectorDouble(table, "size_y_mm", detector.size_y_mm);
	LoadDetectorDouble(table, "center_x_mm", detector.center_x_mm);
	LoadDetectorDouble(table, "center_y_mm", detector.center_y_mm);
	LoadDetectorDouble(table, "z_mm", detector.z_mm);
	LoadDetectorDouble(table, "match_tolerance", detector.match_tolerance);
}

} // namespace

int LoadConfig(const std::string &path, AppConfig &config) {
	try {
		config = AppConfig();
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
		if (auto node = table["assets"].value<std::string>()) {
			config.assets = *node;
		}
		if (const auto *jump_run = table["jump_run"].as_array()) {
			for (size_t i = 0; i < jump_run->size(); ++i) {
				if (auto value = (*jump_run)[i].value<int>()) {
					config.jump_run.push_back(*value);
				}
			}
		}

		if (const auto *paths = table["paths"].as_table()) {
			LoadPaths(*paths, config.paths);
		}
		if (const auto *t0 = table["t0"].as_table()) {
			LoadT0(*t0, config.t0);
		}
		if (const auto *normalize = table["normalize"].as_table()) {
			LoadNormalize(*normalize, config.normalize);
		}
		if (const auto *track = table["track"].as_table()) {
			LoadTrack(*track, config.track);
		}
		if (const auto *identify = table["identify"].as_table()) {
			LoadIdentify(*identify, config.identify);
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

const StraightSliceConfig *FindStraightSliceConfig(
	const AppConfig &config,
	const std::string &name
) {
	auto iter = config.identify.straight.find(name);
	if (iter == config.identify.straight.end()) return nullptr;
	return &iter->second;
}

} // namespace brill