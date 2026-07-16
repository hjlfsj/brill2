#include "include/event/ppac/ppac_track_event.h"

namespace brill {

void SetupInput(TTree *tree, PpacTrackEvent &event, const std::string &prefix) {
	tree->SetBranchAddress((prefix + "valid").c_str(), &event.valid);
	tree->SetBranchAddress((prefix + "x_used_ppac").c_str(), &event.x_used_ppac);
	tree->SetBranchAddress((prefix + "y_used_ppac").c_str(), &event.y_used_ppac);
	tree->SetBranchAddress((prefix + "target_x").c_str(), &event.target_x);
	tree->SetBranchAddress((prefix + "target_y").c_str(), &event.target_y);
	tree->SetBranchAddress((prefix + "dir_x").c_str(), &event.dir_x);
	tree->SetBranchAddress((prefix + "dir_y").c_str(), &event.dir_y);
	tree->SetBranchAddress((prefix + "ppac_x").c_str(), event.ppac_x);
	tree->SetBranchAddress((prefix + "ppac_y").c_str(), event.ppac_y);
}

void SetupOutput(TTree *tree, PpacTrackEvent &event) {
	tree->Branch("valid", &event.valid, "valid/I");
	tree->Branch("x_used_ppac", &event.x_used_ppac, "x_used_ppac/i");
	tree->Branch("y_used_ppac", &event.y_used_ppac, "y_used_ppac/i");
	tree->Branch("target_x", &event.target_x, "target_x/D");
	tree->Branch("target_y", &event.target_y, "target_y/D");
	tree->Branch("dir_x", &event.dir_x, "dir_x/D");
	tree->Branch("dir_y", &event.dir_y, "dir_y/D");
	tree->Branch("ppac_x", event.ppac_x, "ppac_x[3]/D");
	tree->Branch("ppac_y", event.ppac_y, "ppac_y[3]/D");
}

} // namespace brill