#pragma once

#include <cmath>
#include <string>

#include <TTree.h>

#include "include/config.h"

namespace brill {

struct PpacTrackEvent {
	int valid = 0;
	unsigned int x_used_ppac = 0;
	unsigned int y_used_ppac = 0;
	double target_x = 0.0;
	double target_y = 0.0;
	double dir_x = 0.0;
	double dir_y = 0.0;
	double ppac_x[kMaxPpac] = {NAN, NAN, NAN};
	double ppac_y[kMaxPpac] = {NAN, NAN, NAN};
};

void SetupInput(TTree *tree, PpacTrackEvent &event, const std::string &prefix = "");
void SetupOutput(TTree *tree, PpacTrackEvent &event);

} // namespace brill