#include "include/three_alpha/three_alpha_event.h"

namespace brill {

void SetupInputThreeAlpha(TTree *tree, ThreeAlphaEvent &event) {
	tree->SetBranchAddress("run_number", &event.run_number);
	tree->SetBranchAddress("entry", &event.entry);
	tree->SetBranchAddress("e1_He1", &event.e1_He1);
	tree->SetBranchAddress("e2_He1", &event.e2_He1);
	tree->SetBranchAddress("e3_He1", &event.e3_He1);
	tree->SetBranchAddress("e4_He1", &event.e4_He1);
	tree->SetBranchAddress("e1_He2", &event.e1_He2);
	tree->SetBranchAddress("e2_He2", &event.e2_He2);
	tree->SetBranchAddress("e3_He2", &event.e3_He2);
	tree->SetBranchAddress("e4_He2", &event.e4_He2);
	tree->SetBranchAddress("e1_He3", &event.e1_He3);
	tree->SetBranchAddress("e2_He3", &event.e2_He3);
	tree->SetBranchAddress("e3_He3", &event.e3_He3);
	tree->SetBranchAddress("e4_He3", &event.e4_He3);
	tree->SetBranchAddress("et_He1", &event.et_He1);
	tree->SetBranchAddress("et_He2", &event.et_He2);
	tree->SetBranchAddress("et_He3", &event.et_He3);
	tree->SetBranchAddress("stop_He1", &event.stop_He1);
	tree->SetBranchAddress("stop_He2", &event.stop_He2);
	tree->SetBranchAddress("stop_He3", &event.stop_He3);
	tree->SetBranchAddress("is_14O", &event.is_14O);
	tree->SetBranchAddress("is_13N", &event.is_13N);
	tree->SetBranchAddress("is_12C", &event.is_12C);
	tree->SetBranchAddress("ppac_valid", &event.ppac_valid);
	tree->SetBranchAddress("target_x", &event.target_x);
	tree->SetBranchAddress("target_y", &event.target_y);
	tree->SetBranchAddress("dir_x", &event.dir_x);
	tree->SetBranchAddress("dir_y", &event.dir_y);
	tree->SetBranchAddress("t0d2_He1_x", &event.t0d2_He1_x);
	tree->SetBranchAddress("t0d2_He1_y", &event.t0d2_He1_y);
	tree->SetBranchAddress("t0d2_He1_z", &event.t0d2_He1_z);
	tree->SetBranchAddress("t0d2_He2_x", &event.t0d2_He2_x);
	tree->SetBranchAddress("t0d2_He2_y", &event.t0d2_He2_y);
	tree->SetBranchAddress("t0d2_He2_z", &event.t0d2_He2_z);
	tree->SetBranchAddress("t0d2_He3_x", &event.t0d2_He3_x);
	tree->SetBranchAddress("t0d2_He3_y", &event.t0d2_He3_y);
	tree->SetBranchAddress("t0d2_He3_z", &event.t0d2_He3_z);
	tree->SetBranchAddress("theta_beam", &event.theta_beam);
	tree->SetBranchAddress("phi_beam", &event.phi_beam);
	tree->SetBranchAddress("theta_He1", &event.theta_He1);
	tree->SetBranchAddress("theta_He2", &event.theta_He2);
	tree->SetBranchAddress("theta_He3", &event.theta_He3);
	tree->SetBranchAddress("opening_He1_He2", &event.opening_He1_He2);
	tree->SetBranchAddress("opening_He1_He3", &event.opening_He1_He3);
	tree->SetBranchAddress("opening_He2_He3", &event.opening_He2_He3);
	tree->SetBranchAddress("e_sum_He1", &event.e_sum_He1);
	tree->SetBranchAddress("e_sum_He2", &event.e_sum_He2);
	tree->SetBranchAddress("e_sum_He3", &event.e_sum_He3);
}

void SetupOutputThreeAlpha(TTree *tree, ThreeAlphaEvent &event) {
	tree->Branch("run_number", &event.run_number, "run_number/I");
	tree->Branch("entry", &event.entry, "entry/L");
	tree->Branch("e1_He1", &event.e1_He1, "e1_He1/D");
	tree->Branch("e2_He1", &event.e2_He1, "e2_He1/D");
	tree->Branch("e3_He1", &event.e3_He1, "e3_He1/D");
	tree->Branch("e4_He1", &event.e4_He1, "e4_He1/D");
	tree->Branch("e1_He2", &event.e1_He2, "e1_He2/D");
	tree->Branch("e2_He2", &event.e2_He2, "e2_He2/D");
	tree->Branch("e3_He2", &event.e3_He2, "e3_He2/D");
	tree->Branch("e4_He2", &event.e4_He2, "e4_He2/D");
	tree->Branch("e1_He3", &event.e1_He3, "e1_He3/D");
	tree->Branch("e2_He3", &event.e2_He3, "e2_He3/D");
	tree->Branch("e3_He3", &event.e3_He3, "e3_He3/D");
	tree->Branch("e4_He3", &event.e4_He3, "e4_He3/D");
	tree->Branch("et_He1", &event.et_He1, "et_He1/D");
	tree->Branch("et_He2", &event.et_He2, "et_He2/D");
	tree->Branch("et_He3", &event.et_He3, "et_He3/D");
	tree->Branch("stop_He1", &event.stop_He1, "stop_He1/O");
	tree->Branch("stop_He2", &event.stop_He2, "stop_He2/O");
	tree->Branch("stop_He3", &event.stop_He3, "stop_He3/O");
	tree->Branch("is_14O", &event.is_14O, "is_14O/O");
	tree->Branch("is_13N", &event.is_13N, "is_13N/O");
	tree->Branch("is_12C", &event.is_12C, "is_12C/O");
	tree->Branch("ppac_valid", &event.ppac_valid, "ppac_valid/O");
	tree->Branch("target_x", &event.target_x, "target_x/D");
	tree->Branch("target_y", &event.target_y, "target_y/D");
	tree->Branch("dir_x", &event.dir_x, "dir_x/D");
	tree->Branch("dir_y", &event.dir_y, "dir_y/D");
	tree->Branch("t0d2_He1_x", &event.t0d2_He1_x, "t0d2_He1_x/D");
	tree->Branch("t0d2_He1_y", &event.t0d2_He1_y, "t0d2_He1_y/D");
	tree->Branch("t0d2_He1_z", &event.t0d2_He1_z, "t0d2_He1_z/D");
	tree->Branch("t0d2_He2_x", &event.t0d2_He2_x, "t0d2_He2_x/D");
	tree->Branch("t0d2_He2_y", &event.t0d2_He2_y, "t0d2_He2_y/D");
	tree->Branch("t0d2_He2_z", &event.t0d2_He2_z, "t0d2_He2_z/D");
	tree->Branch("t0d2_He3_x", &event.t0d2_He3_x, "t0d2_He3_x/D");
	tree->Branch("t0d2_He3_y", &event.t0d2_He3_y, "t0d2_He3_y/D");
	tree->Branch("t0d2_He3_z", &event.t0d2_He3_z, "t0d2_He3_z/D");
	tree->Branch("theta_beam", &event.theta_beam, "theta_beam/D");
	tree->Branch("phi_beam", &event.phi_beam, "phi_beam/D");
	tree->Branch("theta_He1", &event.theta_He1, "theta_He1/D");
	tree->Branch("theta_He2", &event.theta_He2, "theta_He2/D");
	tree->Branch("theta_He3", &event.theta_He3, "theta_He3/D");
	tree->Branch("opening_He1_He2", &event.opening_He1_He2, "opening_He1_He2/D");
	tree->Branch("opening_He1_He3", &event.opening_He1_He3, "opening_He1_He3/D");
	tree->Branch("opening_He2_He3", &event.opening_He2_He3, "opening_He2_He3/D");
	tree->Branch("e_sum_He1", &event.e_sum_He1, "e_sum_He1/D");
	tree->Branch("e_sum_He2", &event.e_sum_He2, "e_sum_He2/D");
	tree->Branch("e_sum_He3", &event.e_sum_He3, "e_sum_He3/D");
}

} // namespace brill