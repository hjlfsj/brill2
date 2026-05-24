// #include "include/event/ingot/csi_event.h"

// namespace brill {

// void SetupInput(TTree *tree, CsiEvent &event, const std::string &prefix) {
// 	tree->SetBranchAddress((prefix+"flag").c_str(), &event.flag);
// 	tree->SetBranchAddress((prefix+"valid").c_str(), event.valid);
// 	tree->SetBranchAddress((prefix+"time").c_str(), event.time);
// 	tree->SetBranchAddress((prefix+"energy").c_str(), event.energy);
// }

// void SetupOutput(TTree *tree, CsiEvent &event) {
// 	tree->Branch("flag", &event.flag, "flag/l");
// 	tree->Branch("valid", event.valid, "v[36]/O");
// 	tree->Branch("time", event.time, "t[36]/D");
// 	tree->Branch("energy", event.energy, "e[36]/I");
// }

// void Reset(CsiEvent &event) {
// 	event.flag = 0;
// 	for (int i = 0; i < 36; ++i) event.valid[i] = false;
// }

// void Update(
// 	CsiEvent &event,
// 	const int index,
// 	const int energy,
// 	const double time
// ) {
// 	event.flag |= 1ll << index;
// 	event.valid[index] = true;
// 	event.energy[index] = energy;
// 	event.time[index] = time;
// }

// }
