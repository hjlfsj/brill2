#ifndef __CSI_EVENT_H__
#define __CSI_EVENT_H__

#include <string>

#include <TTree.h>

namespace glimmer {

struct CsiEvent {
	unsigned long long flag;
	bool valid[36];
	double time[36];
	int energy[36];
};

void SetupInput(TTree *tree, CsiEvent &event, const std::string &prefix = "");

void SetupOutput(TTree *tree, CsiEvent &event);

void Reset(CsiEvent &event);

void Update(CsiEvent &event, const int index, const int energy, const double time);

}

#endif
