#ifndef __VME_EVENT_H__
#define __VME_EVNET_H__

#include <string>
#include <TTree.h>

namespace glimmer {

struct VmeEvent {
	long long vme_entry;
	long long xia_entry;
	long long vme_ts;
	long long xia_ts;
};

void SetupInput(TTree *tree, VmeEvent &event, const std::string &prefix);

void SetupOutput(TTree *tree, VmeEvent &event);


}

#endif