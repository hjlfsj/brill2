#pragma once

#include <string>

#include <TTree.h>

#include "include/analysis/config.h"

namespace glimmer {

struct PpacHitEvent {
	int valid[kMaxPpac] = {0};
	double x[kMaxPpac] = {0.0};
	double y[kMaxPpac] = {0.0};
};

void SetupInput(TTree *tree, PpacHitEvent &event, const std::string &prefix = "");
void SetupOutput(TTree *tree, PpacHitEvent &event);

} // namespace glimmer
