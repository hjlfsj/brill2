#pragma once

#include <string>

#include <TTree.h>

#include "include/analysis/config.h"

namespace glimmer {

struct T0TrackEvent {
	int num = 0;
	unsigned int flag[kMaxTrackHit] = {0};
	double energy[kMaxTrackHit][kT0LayerCount] = {{0.0}};
	double x[kMaxTrackHit][kT0LayerCount] = {{0.0}};
	double y[kMaxTrackHit][kT0LayerCount] = {{0.0}};
	double z[kMaxTrackHit][kT0LayerCount] = {{0.0}};
};

void SetupInput(TTree *tree, T0TrackEvent &event, const std::string &prefix = "");
void SetupOutput(TTree *tree, T0TrackEvent &event);

} // namespace glimmer
