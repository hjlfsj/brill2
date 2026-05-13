#include "include/analysis/config.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <sstream>

namespace glimmer {

namespace {

std::string Trim(const std::string &text) {
	size_t left = 0;
	while (left < text.size() && std::isspace(static_cast<unsigned char>(text[left]))) {
		++left;
	}
	size_t right = text.size();
	while (right > left && std::isspace(static_cast<unsigned char>(text[right - 1]))) {
		--right;
	}
	return text.substr(left, right - left);
}

std::string RemoveComment(const std::string &line) {
	bool in_quote = false;
	for (size_t i = 0; i < line.size(); ++i) {
		if (line[i] == '"') in_quote = !in_quote;
		if (!in_quote && line[i] == '#') return line.substr(0, i);
	}
	return line;
}

std::string Unquote(const std::string &value) {
	if (value.size() >= 2 && value.front() == '"' && value.back() == '"') {
		return value.substr(1, value.size() - 2);
	}
	return value;
}

bool GetString(
	const std::map<std::string, std::string> &values,
	const std::string &key,
	std::string &value
) {
	auto iter = values.find(key);
	if (iter == values.end()) return false;
	value = Unquote(iter->second);
	return true;
}

bool GetInt(
	const std::map<std::string, std::string> &values,
	const std::string &key,
	int &value
) {
	auto iter = values.find(key);
	if (iter == values.end()) return false;
	value = std::stoi(iter->second);
	return true;
}

bool GetDouble(
	const std::map<std::string, std::string> &values,
	const std::string &key,
	double &value
) {
	auto iter = values.find(key);
	if (iter == values.end()) return false;
	value = std::stod(iter->second);
	return true;
}

void LoadDetector(const std::map<std::string, std::string> &values, const std::string &name, SquareDetectorConfig &detector) {
	detector.name = name;
	GetInt(values, "detectors." + name + ".front_strips", detector.front_strips);
	GetInt(values, "detectors." + name + ".back_strips", detector.back_strips);
	GetDouble(values, "detectors." + name + ".thickness_um", detector.thickness_um);
	GetDouble(values, "detectors." + name + ".size_x_mm", detector.size_x_mm);
	GetDouble(values, "detectors." + name + ".size_y_mm", detector.size_y_mm);
	GetDouble(values, "detectors." + name + ".center_x_mm", detector.center_x_mm);
	GetDouble(values, "detectors." + name + ".center_y_mm", detector.center_y_mm);
	GetDouble(values, "detectors." + name + ".z_mm", detector.z_mm);
	GetDouble(values, "detectors." + name + ".merge_tolerance", detector.merge_tolerance);
	GetDouble(values, "detectors." + name + ".track_window_x", detector.track_window_x);
	GetDouble(values, "detectors." + name + ".track_window_y", detector.track_window_y);
}

} // namespace

int LoadConfig(const std::string &path, AppConfig &config) {
	std::ifstream fin(path);
	if (!fin.good()) {
		std::cerr << "Error: Open config file " << path << " failed.\n";
		return -1;
	}

	std::map<std::string, std::string> values;
	std::string line;
	std::string section;
	while (std::getline(fin, line)) {
		line = Trim(RemoveComment(line));
		if (line.empty()) continue;
		if (line.front() == '[' && line.back() == ']') {
			section = Trim(line.substr(1, line.size() - 2));
			continue;
		}
		size_t split = line.find('=');
		if (split == std::string::npos) continue;
		std::string key = Trim(line.substr(0, split));
		std::string value = Trim(line.substr(split + 1));
		if (!section.empty()) key = section + "." + key;
		values[key] = value;
	}

	GetString(values, "workspace", config.workspace);
	GetString(values, "analysis.trigger_default", config.trigger_default);
	GetString(values, "analysis.paths.forge", config.paths.forge);
	GetString(values, "analysis.paths.normalize", config.paths.normalize);
	GetString(values, "analysis.paths.merge", config.paths.merge);
	GetString(values, "analysis.paths.track", config.paths.track);
	GetString(values, "analysis.paths.show", config.paths.show);

	for (int i = 0; i < kMaxPpac; ++i) {
		GetDouble(values, "detectors.ppac.z" + std::to_string(i) + "_mm", config.ppac.z_mm[i]);
		GetDouble(values, "detectors.ppac.x_offset" + std::to_string(i) + "_mm", config.ppac.x_offset_mm[i]);
		GetDouble(values, "detectors.ppac.y_offset" + std::to_string(i) + "_mm", config.ppac.y_offset_mm[i]);
		GetDouble(values, "detectors.ppac.x_scale" + std::to_string(i), config.ppac.x_scale[i]);
		GetDouble(values, "detectors.ppac.y_scale" + std::to_string(i), config.ppac.y_scale[i]);
	}

	const char *names[] = {"t0d1", "t0d2", "t0d3", "t0d4", "t0s"};
	for (const char *name : names) {
		SquareDetectorConfig detector;
		LoadDetector(values, name, detector);
		config.detectors[detector.name] = detector;
	}
	return 0;
}

const SquareDetectorConfig *FindDetectorConfig(const AppConfig &config, const std::string &name) {
	auto iter = config.detectors.find(name);
	if (iter == config.detectors.end()) return nullptr;
	return &iter->second;
}

} // namespace glimmer
