#include "include/event/beam/beam_sort.h"

#include <TF1.h>
#include <TSpectrum.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace brill {

int SortBeamTOF(TH1D *h_tof, BeamSortResult &result) {
	if (!h_tof || h_tof->GetEntries() == 0) return 1;

	TSpectrum spectrum(10);
	int n_peaks = spectrum.Search(h_tof, 2, "", 0.10);

	if (n_peaks < 3) {
		fprintf(stderr, "Error: SortBeamTOF found only %d peaks, need at least 3.\n", n_peaks);
		return 2;
	}

	std::vector<int> peak_indices(n_peaks);
	for (int i = 0; i < n_peaks; ++i) {
		peak_indices[i] = i;
	}

	std::sort(peak_indices.begin(), peak_indices.end(), [&](int a, int b) {
		return spectrum.GetPositionY()[a] > spectrum.GetPositionY()[b];
	});

	std::vector<int> top3 = {peak_indices[0], peak_indices[1], peak_indices[2]};
	std::sort(top3.begin(), top3.end(), [&](int a, int b) {
		return spectrum.GetPositionX()[a] < spectrum.GetPositionX()[b];
	});

	double pos[3] = {
		spectrum.GetPositionX()[top3[0]],
		spectrum.GetPositionX()[top3[1]],
		spectrum.GetPositionX()[top3[2]]
	};
	double hgt[3] = {
		spectrum.GetPositionY()[top3[0]],
		spectrum.GetPositionY()[top3[1]],
		spectrum.GetPositionY()[top3[2]]
	};

	printf("Peaks found (left to right):\n");
	printf("  14O: pos=%.2f, height=%.0f\n", pos[0], hgt[0]);
	printf("  13N: pos=%.2f, height=%.0f\n", pos[1], hgt[1]);
	printf("  12C: pos=%.2f, height=%.0f\n", pos[2], hgt[2]);

	double *means[3] = {&result.mean_14O, &result.mean_13N, &result.mean_12C};
	double *sigmas[3] = {&result.sigma_14O, &result.sigma_13N, &result.sigma_12C};
	double *x_lows[3] = {&result.x_low_14O, &result.x_low_13N, &result.x_low_12C};
	double *x_highs[3] = {&result.x_high_14O, &result.x_high_13N, &result.x_high_12C};
	const char *labels[3] = {"14O", "13N", "12C"};

	for (int i = 0; i < 3; ++i) {
		double fit_min = pos[i] - 3.0;
		double fit_max = pos[i] + 3.0;

		TString func_name = TString::Format("gaus_%s", labels[i]);
		TF1 *gaus = new TF1(func_name, "gaus", fit_min, fit_max);
		gaus->SetParameter(0, h_tof->GetBinContent(h_tof->FindBin(pos[i])));
		gaus->SetParameter(1, pos[i]);
		gaus->SetParameter(2, 5.0);

		h_tof->Fit(gaus, "RQ0");

		*means[i] = gaus->GetParameter(1);
		*sigmas[i] = std::abs(gaus->GetParameter(2));
		*x_lows[i] = *means[i] - 5.0 * (*sigmas[i]);
		*x_highs[i] = *means[i] + 5.0 * (*sigmas[i]);

		printf("  %s fit: mean=%.2f, sigma=%.2f, 5sigma range: [%.2f, %.2f]\n",
			labels[i], *means[i], *sigmas[i], *x_lows[i], *x_highs[i]);
	}

	for (int i = 0; i < 2; ++i) {
		if (*x_highs[i] > *x_lows[i+1]) {
			double mid = (*means[i] + *means[i+1]) / 2.0;
			printf("  Overlap: %s and %s 5sigma ranges overlap, adjusted boundary to %.2f\n",
				labels[i], labels[i+1], mid);
			*x_highs[i] = mid;
			*x_lows[i+1] = mid;
		} else {
			printf("  No overlap: %s and %s, gap=%.2f\n",
				labels[i], labels[i+1], *x_lows[i+1] - *x_highs[i]);
		}
	}

	return 0;
}

void SetupInputSortBeamTree(TTree *tree, bool &v_14O, bool &v_13N, bool &v_12C) {
	tree->SetBranchAddress("14O_valid", &v_14O);
	tree->SetBranchAddress("13N_valid", &v_13N);
	tree->SetBranchAddress("12C_valid", &v_12C);
}

void SetupOutputSortBeamTree(TTree *tree, bool &v_14O, bool &v_13N, bool &v_12C) {
	tree->Branch("14O_valid", &v_14O, "14O_valid/O");
	tree->Branch("13N_valid", &v_13N, "13N_valid/O");
	tree->Branch("12C_valid", &v_12C, "12C_valid/O");
}

void ClassifyBeam(
	double tof,
	bool valid,
	const BeamSortResult &result,
	bool &v_14O,
	bool &v_13N,
	bool &v_12C
) {
	v_14O = false;
	v_13N = false;
	v_12C = false;

	if (!valid) return;

	if (tof >= result.x_low_14O && tof <= result.x_high_14O) {
		v_14O = true;
	} else if (tof >= result.x_low_13N && tof <= result.x_high_13N) {
		v_13N = true;
	} else if (tof >= result.x_low_12C && tof <= result.x_high_12C) {
		v_12C = true;
	}
}

} // namespace brill