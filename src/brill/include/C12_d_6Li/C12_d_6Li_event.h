#pragma once

#include <TTree.h>

namespace brill {

struct C12D6LiEvent {
	int run_number = 0;
	Long64_t entry = 0;

	double e1_6Li = 0.0;
	double e2_6Li = 0.0;

	double e1_4He1 = 0.0;
	double e2_4He1 = 0.0;
	double e3_4He1 = 0.0;
	double e4_4He1 = 0.0;

	double e1_4He2 = 0.0;
	double e2_4He2 = 0.0;
	double e3_4He2 = 0.0;
	double e4_4He2 = 0.0;

	double e5 = 0.0;

	bool is_14O = false;
	bool is_13N = false;
	bool is_12C = false;

	bool ppac_valid = false;
	double target_x = 0.0;
	double target_y = 0.0;
	double dir_x = 0.0;
	double dir_y = 0.0;

	double t0d2_6Li_x = 0.0;
	double t0d2_6Li_y = 0.0;
	double t0d2_6Li_z = 0.0;
	double t0d2_4He1_x = 0.0;
	double t0d2_4He1_y = 0.0;
	double t0d2_4He1_z = 0.0;
	double t0d2_4He2_x = 0.0;
	double t0d2_4He2_y = 0.0;
	double t0d2_4He2_z = 0.0;

	double theta_beam = 0.0;
	double phi_beam = 0.0;
	double theta_6Li = 0.0;
	double theta_4He1 = 0.0;
	double theta_4He2 = 0.0;
	double opening_6Li_4He1 = 0.0;
	double opening_6Li_4He2 = 0.0;
	double opening_4He1_4He2 = 0.0;
};

void SetupInputC12D6Li(TTree *tree, C12D6LiEvent &event);
void SetupOutputC12D6Li(TTree *tree, C12D6LiEvent &event);

} // namespace brill