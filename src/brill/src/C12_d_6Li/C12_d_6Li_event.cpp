#include "include/C12_d_6Li/C12_d_6Li_event.h"

namespace brill {

void SetupInputC12D6Li(TTree *tree, C12D6LiEvent &event) {
	tree->SetBranchAddress("run_number", &event.run_number);
	tree->SetBranchAddress("entry", &event.entry);
	tree->SetBranchAddress("e1_6Li", &event.e1_6Li);
	tree->SetBranchAddress("e2_6Li", &event.e2_6Li);
	tree->SetBranchAddress("e1_4He1", &event.e1_4He1);
	tree->SetBranchAddress("e2_4He1", &event.e2_4He1);
	tree->SetBranchAddress("e3_4He1", &event.e3_4He1);
	tree->SetBranchAddress("e4_4He1", &event.e4_4He1);
	tree->SetBranchAddress("e1_4He2", &event.e1_4He2);
	tree->SetBranchAddress("e2_4He2", &event.e2_4He2);
	tree->SetBranchAddress("e3_4He2", &event.e3_4He2);
	tree->SetBranchAddress("e4_4He2", &event.e4_4He2);
	tree->SetBranchAddress("e5", &event.e5);
	tree->SetBranchAddress("is_14O", &event.is_14O);
	tree->SetBranchAddress("is_13N", &event.is_13N);
	tree->SetBranchAddress("is_12C", &event.is_12C);
	tree->SetBranchAddress("ppac_valid", &event.ppac_valid);
	tree->SetBranchAddress("target_x", &event.target_x);
	tree->SetBranchAddress("target_y", &event.target_y);
	tree->SetBranchAddress("dir_x", &event.dir_x);
	tree->SetBranchAddress("dir_y", &event.dir_y);
	tree->SetBranchAddress("t0d2_6Li_x", &event.t0d2_6Li_x);
	tree->SetBranchAddress("t0d2_6Li_y", &event.t0d2_6Li_y);
	tree->SetBranchAddress("t0d2_6Li_z", &event.t0d2_6Li_z);
	tree->SetBranchAddress("t0d2_4He1_x", &event.t0d2_4He1_x);
	tree->SetBranchAddress("t0d2_4He1_y", &event.t0d2_4He1_y);
	tree->SetBranchAddress("t0d2_4He1_z", &event.t0d2_4He1_z);
	tree->SetBranchAddress("t0d2_4He2_x", &event.t0d2_4He2_x);
	tree->SetBranchAddress("t0d2_4He2_y", &event.t0d2_4He2_y);
	tree->SetBranchAddress("t0d2_4He2_z", &event.t0d2_4He2_z);
	tree->SetBranchAddress("theta_beam", &event.theta_beam);
	tree->SetBranchAddress("phi_beam", &event.phi_beam);
	tree->SetBranchAddress("theta_6Li", &event.theta_6Li);
	tree->SetBranchAddress("theta_4He1", &event.theta_4He1);
	tree->SetBranchAddress("theta_4He2", &event.theta_4He2);
	tree->SetBranchAddress("opening_6Li_4He1", &event.opening_6Li_4He1);
	tree->SetBranchAddress("opening_6Li_4He2", &event.opening_6Li_4He2);
	tree->SetBranchAddress("opening_4He1_4He2", &event.opening_4He1_4He2);
}

void SetupOutputC12D6Li(TTree *tree, C12D6LiEvent &event) {
	tree->Branch("run_number", &event.run_number, "run_number/I");
	tree->Branch("entry", &event.entry, "entry/L");
	tree->Branch("e1_6Li", &event.e1_6Li, "e1_6Li/D");
	tree->Branch("e2_6Li", &event.e2_6Li, "e2_6Li/D");
	tree->Branch("e1_4He1", &event.e1_4He1, "e1_4He1/D");
	tree->Branch("e2_4He1", &event.e2_4He1, "e2_4He1/D");
	tree->Branch("e3_4He1", &event.e3_4He1, "e3_4He1/D");
	tree->Branch("e4_4He1", &event.e4_4He1, "e4_4He1/D");
	tree->Branch("e1_4He2", &event.e1_4He2, "e1_4He2/D");
	tree->Branch("e2_4He2", &event.e2_4He2, "e2_4He2/D");
	tree->Branch("e3_4He2", &event.e3_4He2, "e3_4He2/D");
	tree->Branch("e4_4He2", &event.e4_4He2, "e4_4He2/D");
	tree->Branch("e5", &event.e5, "e5/D");
	tree->Branch("is_14O", &event.is_14O, "is_14O/O");
	tree->Branch("is_13N", &event.is_13N, "is_13N/O");
	tree->Branch("is_12C", &event.is_12C, "is_12C/O");
	tree->Branch("ppac_valid", &event.ppac_valid, "ppac_valid/O");
	tree->Branch("target_x", &event.target_x, "target_x/D");
	tree->Branch("target_y", &event.target_y, "target_y/D");
	tree->Branch("dir_x", &event.dir_x, "dir_x/D");
	tree->Branch("dir_y", &event.dir_y, "dir_y/D");
	tree->Branch("t0d2_6Li_x", &event.t0d2_6Li_x, "t0d2_6Li_x/D");
	tree->Branch("t0d2_6Li_y", &event.t0d2_6Li_y, "t0d2_6Li_y/D");
	tree->Branch("t0d2_6Li_z", &event.t0d2_6Li_z, "t0d2_6Li_z/D");
	tree->Branch("t0d2_4He1_x", &event.t0d2_4He1_x, "t0d2_4He1_x/D");
	tree->Branch("t0d2_4He1_y", &event.t0d2_4He1_y, "t0d2_4He1_y/D");
	tree->Branch("t0d2_4He1_z", &event.t0d2_4He1_z, "t0d2_4He1_z/D");
	tree->Branch("t0d2_4He2_x", &event.t0d2_4He2_x, "t0d2_4He2_x/D");
	tree->Branch("t0d2_4He2_y", &event.t0d2_4He2_y, "t0d2_4He2_y/D");
	tree->Branch("t0d2_4He2_z", &event.t0d2_4He2_z, "t0d2_4He2_z/D");
	tree->Branch("theta_beam", &event.theta_beam, "theta_beam/D");
	tree->Branch("phi_beam", &event.phi_beam, "phi_beam/D");
	tree->Branch("theta_6Li", &event.theta_6Li, "theta_6Li/D");
	tree->Branch("theta_4He1", &event.theta_4He1, "theta_4He1/D");
	tree->Branch("theta_4He2", &event.theta_4He2, "theta_4He2/D");
	tree->Branch("opening_6Li_4He1", &event.opening_6Li_4He1, "opening_6Li_4He1/D");
	tree->Branch("opening_6Li_4He2", &event.opening_6Li_4He2, "opening_6Li_4He2/D");
	tree->Branch("opening_4He1_4He2", &event.opening_4He1_4He2, "opening_4He1_4He2/D");
}

} // namespace brill