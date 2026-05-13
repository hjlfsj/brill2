#ifndef __PPAC_EVENT_H__
#define __PPAC_EVENT_H__

#include <string>
#include <TTree.h>

namespace glimmer {

struct PpacEvent {
	int flag;
	double time[15];
	int energy[15];
	// bool cv[15];
};

void SetupInput(TTree *tree, PpacEvent &event, const std::string &prefix);

void SetupOutput(TTree *tree, PpacEvent &event);

}

#endif