#pragma once

#include <TH1D.h>
#include <TTree.h>

namespace brill {

struct BeamSortResult {
	double mean_14O = 0.0;
	double sigma_14O = 0.0;
	double x_low_14O = 0.0;
	double x_high_14O = 0.0;

	double mean_13N = 0.0;
	double sigma_13N = 0.0;
	double x_low_13N = 0.0;
	double x_high_13N = 0.0;

	double mean_12C = 0.0;
	double sigma_12C = 0.0;
	double x_low_12C = 0.0;
	double x_high_12C = 0.0;
};

int SortBeamTOF(TH1D *h_tof, BeamSortResult &result);

void SetupInputSortBeamTree(TTree *tree, bool &v_14O, bool &v_13N, bool &v_12C);

void SetupOutputSortBeamTree(TTree *tree, bool &v_14O, bool &v_13N, bool &v_12C);

void ClassifyBeam(
	double tof,
	bool valid,
	const BeamSortResult &result,
	bool &v_14O,
	bool &v_13N,
	bool &v_12C
);

} // namespace brill