#include "include/event/forge/beam_event.h"

namespace glimmer {

void SetupInput(TTree *tree, BeamEvent &event, const std::string &prefix) {
	tree->SetBranchAddress((prefix+"flag").c_str(), &event.flag);
	tree->SetBranchAddress((prefix+"valid").c_str(), &event.valid);
	tree->SetBranchAddress((prefix+"tof").c_str(), &event.tof);
}

void SetupOutput(TTree *tree, BeamEvent &event) {
	tree->Branch("flag", &event.flag, "flag/I");
	tree->Branch("valid", &event.valid, "valid/O");
	tree->Branch("tof", &event.tof, "tof/D");
}

}
