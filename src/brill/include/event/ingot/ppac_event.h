#ifndef __PPAC_EVENT_H__
#define __PPAC_EVENT_H__

#include <string>
#include <TTree.h>

namespace brill {

struct PpacEvent {
	int flag;
	bool valid[15];
	double time[15];
	int energy[15];
};

struct PpacChannelMap {
	static constexpr int ch_x1[3] = {0, 4, 8};
	static constexpr int ch_x2[3] = {1, 5, 9};
	static constexpr int ch_y1[3] = {2, 6, 10};
	static constexpr int ch_y2[3] = {3, 7, 11};
	static constexpr int ch_anode[3] = {12, 13, 14};
};

void SetupInput(TTree *tree, PpacEvent &event, const std::string &prefix);

void SetupOutput(TTree *tree, PpacEvent &event);

}

#endif