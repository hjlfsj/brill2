#pragma once

#include <string>

#include "include/analysis/config.h"

namespace glimmer {

std::string JoinPath(const std::string &left, const std::string &right);
std::string TriggerSuffix(const std::string &trigger);
std::string ForgeFileName(
	const AppConfig &config,
	const std::string &detector,
	const std::string &trigger,
	int run
);
std::string NormalizeFileName(
	const AppConfig &config,
	const std::string &detector,
	const std::string &trigger,
	int run
);
std::string NormalizeParamFileName(
	const AppConfig &config,
	const std::string &detector,
	const std::string &trigger,
	int run,
	int end_run
);
std::string MergeFileName(
	const AppConfig &config,
	const std::string &detector,
	const std::string &trigger,
	int run
);
std::string TrackFileName(
	const AppConfig &config,
	const std::string &name,
	const std::string &trigger,
	int run
);
std::string ShowFileName(
	const AppConfig &config,
	const std::string &name,
	const std::string &trigger,
	int run
);

} // namespace glimmer
