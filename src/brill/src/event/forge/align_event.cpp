#include "include/event/forge/align_event.h"

namespace glimmer {

void SetupInput(TTree *tree, AlignEvent &event, const std::string &prefix) {
	tree->SetBranchAddress((prefix+"vme_entry").c_str(), &event.vme_entry);
	tree->SetBranchAddress((prefix+"xia_entry").c_str(), &event.xia_entry);
	tree->SetBranchAddress((prefix+"vme_time").c_str(), &event.vme_time);
	tree->SetBranchAddress((prefix+"xia_time").c_str(), &event.xia_time);
}

void SetupOutput(TTree *tree, AlignEvent &event) {
	tree->Branch("vme_entry", &event.vme_entry, "ve/L");
	tree->Branch("xia_entry", &event.xia_entry, "xe/L");
	tree->Branch("vme_time", &event.vme_time, "vt/L");
	tree->Branch("xia_time", &event.xia_time, "xt/L");
}



}