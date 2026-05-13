#pragma once

#include <string>

#include <TTree.h>

namespace glimmer {

struct BeamEvent {
	int flag;
	bool valid;
	double tof;
};

void SetupInput(TTree *tree, BeamEvent &event, const std::string &prefix = "");

void SetupOutput(TTree *tree, BeamEvent &event);

}
