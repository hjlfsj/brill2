#include "include/event/ppac/ppac_hit_event.h"

#include "include/config.h"

namespace brill {

void SetupInput(TTree *tree, PpacHitEvent &event, const std::string &prefix) {
	tree->SetBranchAddress((prefix + "valid").c_str(), event.valid);
	tree->SetBranchAddress((prefix + "x").c_str(), event.x);
	tree->SetBranchAddress((prefix + "y").c_str(), event.y);
}

void SetupOutput(TTree *tree, PpacHitEvent &event) {
	tree->Branch("valid", event.valid, "valid[3]/I");
	tree->Branch("x", event.x, "x[3]/D");
	tree->Branch("y", event.y, "y[3]/D");
}

} // namespace glimmer
