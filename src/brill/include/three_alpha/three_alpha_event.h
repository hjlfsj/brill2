#pragma once

#include <TTree.h>

namespace brill {

struct ThreeAlphaEvent {
	int run_number = 0;
	Long64_t entry = 0;

	double e1_He1 = 0.0;
	double e2_He1 = 0.0;
	double e3_He1 = 0.0;
	double e4_He1 = 0.0;

	double e1_He2 = 0.0;
	double e2_He2 = 0.0;
	double e3_He2 = 0.0;
	double e4_He2 = 0.0;

	double e1_He3 = 0.0;
	double e2_He3 = 0.0;
	double e3_He3 = 0.0;
	double e4_He3 = 0.0;

	double et_He1 = 0.0;
	double et_He2 = 0.0;
	double et_He3 = 0.0;

	bool stop_He1 = false;
	bool stop_He2 = false;
	bool stop_He3 = false;

	bool is_14O = false;
	bool is_13N = false;
	bool is_12C = false;

	bool ppac_valid = false;
	double target_x = 0.0;
	double target_y = 0.0;
	double dir_x = 0.0;
	double dir_y = 0.0;

	double t0d2_He1_x = 0.0;
	double t0d2_He1_y = 0.0;
	double t0d2_He1_z = 0.0;
	double t0d2_He2_x = 0.0;
	double t0d2_He2_y = 0.0;
	double t0d2_He2_z = 0.0;
	double t0d2_He3_x = 0.0;
	double t0d2_He3_y = 0.0;
	double t0d2_He3_z = 0.0;

	double theta_beam = 0.0;
	double phi_beam = 0.0;
	double theta_He1 = 0.0;
	double theta_He2 = 0.0;
	double theta_He3 = 0.0;
	double opening_He1_He2 = 0.0;
	double opening_He1_He3 = 0.0;
	double opening_He2_He3 = 0.0;

	double e_sum_He1 = 0.0;
	double e_sum_He2 = 0.0;
	double e_sum_He3 = 0.0;
};

void SetupInputThreeAlpha(TTree *tree, ThreeAlphaEvent &event);
void SetupOutputThreeAlpha(TTree *tree, ThreeAlphaEvent &event);

} // namespace brill