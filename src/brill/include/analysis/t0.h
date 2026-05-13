#pragma once

#include <string>

#include "include/analysis/config.h"

namespace glimmer {

int RunT0Track(const AppConfig &config, const std::string &trigger, int run);
int RunEstimatePpacCenter(const AppConfig &config, const std::string &trigger, int run);
int RunEstimateT0Center(const AppConfig &config, const std::string &trigger, int run);

} // namespace glimmer
