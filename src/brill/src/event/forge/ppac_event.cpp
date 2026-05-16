#include "include/event/forge/ppac_event.h"

namespace glimmer {

void SetupInput(TTree *tree, PpacEvent &event, const std::string &prefix) {
	tree->SetBranchAddress((prefix+"flag").c_str(), &event.flag);
	tree->SetBranchAddress((prefix+"valid").c_str(), event.valid);
	tree->SetBranchAddress((prefix+"time").c_str(), event.time);
	tree->SetBranchAddress((prefix+"energy").c_str(), event.energy);
	// tree->SetBranchAddress((prefix+"cv").c_str(), event.cv, "cv[15]/O");
}

void SetupOutput(TTree *tree, PpacEvent &event) {
	tree->Branch("flag", &event.flag, "flag/I");
	tree->Branch("valid", event.valid, "v[15]/O");
	tree->Branch("time", event.time, "t[15]/D");
	tree->Branch("energy", event.energy, "e[15]/I");
}

}