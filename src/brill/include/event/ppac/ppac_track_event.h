#pragma once

#include <string>

#include <TTree.h>

#include "include/config.h"

namespace brill {

struct PpacTrackEvent {
	int valid = 0;
	unsigned int used_ppac = 0;
	double target_x = 0.0;
	double target_y = 0.0;
	double dir_x = 0.0;
	double dir_y = 0.0;
	double ppac_x[kMaxPpac] = {0.0};
	double ppac_y[kMaxPpac] = {0.0};
	int ppac_valid[kMaxPpac] = {0};
};

void SetupInput(TTree *tree, PpacTrackEvent &event, const std::string &prefix = "");
void SetupOutput(TTree *tree, PpacTrackEvent &event);

} // namespace glimmer
