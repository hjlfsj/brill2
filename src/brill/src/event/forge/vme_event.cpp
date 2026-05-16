#include "include/event/forge/vme_event.h"

namespace glimmer {

void SetupInput(TTree *tree, VmeEvent &event, const std::string &prefix) {
	tree->SetBranchAddress((prefix+"vme_entry").c_str(), &event.vme_entry);
	tree->SetBranchAddress((prefix+"xia_entry").c_str(), &event.xia_entry);
	tree->SetBranchAddress((prefix+"vme_ts").c_str(), &event.vme_ts);
	tree->SetBranchAddress((prefix+"xia_ts").c_str(), &event.xia_ts);
}

void SetupOutput(TTree *tree, VmeEvent &event) {
	tree->Branch("vme_entry", &event.vme_entry, "ve/L");
	tree->Branch("xia_entry", &event.xia_entry, "xe/L");
	tree->Branch("vme_ts", &event.vme_ts, "vts/L");
	tree->Branch("xia_ts", &event.xia_ts, "xts/L");
}

}