#include "include/analysis/file_name.h"

namespace glimmer {

std::string JoinPath(const std::string &left, const std::string &right) {
	if (left.empty()) return right;
	if (right.empty()) return left;
	if (left.back() == '/') return left + right;
	return left + "/" + right;
}

std::string TriggerSuffix(const std::string &trigger) {
	return trigger.empty() ? "" : ("_" + trigger);
}

std::string ForgeFileName(const AppConfig &config, const std::string &detector, const std::string &trigger, int run) {
	return JoinPath(
		JoinPath(config.workspace, config.paths.forge),
		detector + TriggerSuffix(trigger) + "_" + (run < 1000 ? (run < 100 ? (run < 10 ? "000" : "00") : "0") : "") + std::to_string(run) + ".root"
	);
}

std::string NormalizeFileName(const AppConfig &config, const std::string &detector, const std::string &trigger, int run) {
	return JoinPath(
		JoinPath(config.workspace, config.paths.normalize),
		detector + "-normalize" + TriggerSuffix(trigger) + "_" + (run < 1000 ? (run < 100 ? (run < 10 ? "000" : "00") : "0") : "") + std::to_string(run) + ".root"
	);
}

std::string NormalizeParamFileName(const AppConfig &config, const std::string &detector, const std::string &trigger, int run, int end_run) {
	return JoinPath(
		JoinPath(config.workspace, config.paths.normalize),
		detector + "-params" + TriggerSuffix(trigger) + "_" +
		(run < 1000 ? (run < 100 ? (run < 10 ? "000" : "00") : "0") : "") + std::to_string(run) +
		"-" +
		(end_run < 1000 ? (end_run < 100 ? (end_run < 10 ? "000" : "00") : "0") : "") + std::to_string(end_run) +
		".root"
	);
}

std::string MergeFileName(const AppConfig &config, const std::string &detector, const std::string &trigger, int run) {
	return JoinPath(
		JoinPath(config.workspace, config.paths.merge),
		detector + "-merge" + TriggerSuffix(trigger) + "_" + (run < 1000 ? (run < 100 ? (run < 10 ? "000" : "00") : "0") : "") + std::to_string(run) + ".root"
	);
}

std::string TrackFileName(const AppConfig &config, const std::string &name, const std::string &trigger, int run) {
	return JoinPath(
		JoinPath(config.workspace, config.paths.track),
		name + TriggerSuffix(trigger) + "_" + (run < 1000 ? (run < 100 ? (run < 10 ? "000" : "00") : "0") : "") + std::to_string(run) + ".root"
	);
}

std::string ShowFileName(const AppConfig &config, const std::string &name, const std::string &trigger, int run) {
	return JoinPath(
		JoinPath(config.workspace, config.paths.show),
		name + TriggerSuffix(trigger) + "_" + (run < 1000 ? (run < 100 ? (run < 10 ? "000" : "00") : "0") : "") + std::to_string(run) + ".root"
	);
}

} // namespace glimmer
