#include "include/event/ppac/adjust_ppac.h"

#include <cmath>
#include <iostream>
#include <vector>

#include <TFile.h>
#include <TTree.h>

namespace brill {

static void LinearFit(
	const double *z, const double *pos, int n,
	double &intercept, double &slope
) {
	double sum_z = 0, sum_pos = 0, sum_zz = 0, sum_zpos = 0;
	for (int i = 0; i < n; ++i) {
		sum_z += z[i];
		sum_pos += pos[i];
		sum_zz += z[i] * z[i];
		sum_zpos += z[i] * pos[i];
	}
	double denom = n * sum_zz - sum_z * sum_z;
	if (std::fabs(denom) < 1e-12) {
		intercept = 0;
		slope = 0;
		return;
	}
	slope = (n * sum_zpos - sum_z * sum_pos) / denom;
	intercept = (sum_pos - slope * sum_z) / n;
}

static double ResidualSum(
	const double *z, const double *pos, int n,
	double intercept, double slope
) {
	double sum = 0;
	for (int i = 0; i < n; ++i) {
		double diff = pos[i] - (intercept + slope * z[i]);
		sum += diff * diff;
	}
	return sum;
}

static double ComputeAvgResidual(
	const std::vector<double> ev_z[3],
	const std::vector<double> ev_pos[3],
	int count,
	double d2, double d3,
	double *var_out = nullptr
) {
	double sum_residual = 0;
	double sum_sq = 0;
	for (int ev = 0; ev < count; ++ev) {
		double z[3] = {ev_z[0][ev], ev_z[1][ev], ev_z[2][ev]};
		double pos[3] = {
			ev_pos[0][ev],
			ev_pos[1][ev] + d2,
			ev_pos[2][ev] + d3
		};
		double intercept, slope;
		LinearFit(z, pos, 3, intercept, slope);
		double res = ResidualSum(z, pos, 3, intercept, slope);
		sum_residual += res;
		sum_sq += res * res;
	}
	double mean = sum_residual / count;
	if (var_out) {
		double mean_sq = sum_sq / count;
		*var_out = mean_sq - mean * mean;
		if (*var_out < 0) *var_out = 0;
	}
	return mean;
}

int AdjustPpacPosition(
	const std::string &track_path,
	const PpacConfig &ppac_config,
	int n,
	PpacPositionAdjust &result
) {
	TFile f(track_path.c_str(), "READ");
	if (!f.IsOpen()) {
		std::cerr << "Error: Cannot open " << track_path << "\n";
		return 1;
	}
	TTree *tree = (TTree*)f.Get("tree");
	if (!tree) {
		std::cerr << "Error: No tree in " << track_path << "\n";
		return 1;
	}

	PpacTrackEvent event;
	SetupInput(tree, event);

	const int kMaxEvents = 1e6;
	long long total = tree->GetEntriesFast();

	std::vector<double> x_events_z[3];
	std::vector<double> x_events_pos[3];
	std::vector<double> y_events_z[3];
	std::vector<double> y_events_pos[3];

	int x_count = 0, y_count = 0;
	for (long long entry = 0; entry < total && x_count < kMaxEvents && y_count < kMaxEvents; ++entry) {
		tree->GetEntry(entry);

		if (event.x_used_ppac == 7 && x_count < kMaxEvents) {
			for (int i = 0; i < 3; ++i) {
				x_events_z[i].push_back(ppac_config.z_x_mm[i]);
				x_events_pos[i].push_back(event.ppac_x[i]);
			}
			++x_count;
		}

		if (event.y_used_ppac == 7 && y_count < kMaxEvents) {
			for (int i = 0; i < 3; ++i) {
				y_events_z[i].push_back(ppac_config.z_y_mm[i]);
				y_events_pos[i].push_back(event.ppac_y[i]);
			}
			++y_count;
		}
	}

	std::cout << "X events: " << x_count << "  Y events: " << y_count << "\n";

	const double kCoarseStep = 0.5;
	const double kFineStep = 0.1;

	auto ScanDirection = [&](
		const std::vector<double> ev_z[3],
		const std::vector<double> ev_pos[3],
		int count,
		double &baseline_residual,
		double &best_d2,
		double &best_d3,
		double &best_residual,
		double &dsigma2,
		double &dsigma3,
		double &improvement
	) {
		dsigma2 = 0;
		dsigma3 = 0;
		improvement = 0;

		baseline_residual = ComputeAvgResidual(ev_z, ev_pos, count, 0.0, 0.0);
		std::cout << "  Baseline residual: " << baseline_residual << "\n";

		// Regularization: penalize large offsets to break the
		// 3-point PPAC degeneracy. The penalty is ~1% of the baseline
		// residual at the edge of the search range.
		double lambda = baseline_residual / (n * n * 100.0);

		best_residual = baseline_residual;
		best_d2 = 0;
		best_d3 = 0;

		// Stage 1: coarse grid search (0.5mm)
		int coarse_steps = static_cast<int>(2 * n / kCoarseStep) + 1;
		for (int i2 = 0; i2 < coarse_steps; ++i2) {
			double d2 = -n + i2 * kCoarseStep;
			for (int i3 = 0; i3 < coarse_steps; ++i3) {
				double d3 = -n + i3 * kCoarseStep;
				double avg_residual = ComputeAvgResidual(ev_z, ev_pos, count, d2, d3);
				double reg = lambda * (d2 * d2 + d3 * d3);
				double total = avg_residual + reg;
				if (total < best_residual) {
					best_residual = total;
					best_d2 = d2;
					best_d3 = d3;
				}
			}
		}
		std::cout << "  Coarse best: d2=" << best_d2 << "  d3=" << best_d3
			<< "  residual=" << ComputeAvgResidual(ev_z, ev_pos, count, best_d2, best_d3) << "\n";

		// Stage 2: fine grid search (0.1mm) around coarse best
		double fine_lo2 = best_d2 - kCoarseStep;
		double fine_hi2 = best_d2 + kCoarseStep;
		double fine_lo3 = best_d3 - kCoarseStep;
		int fine_steps = static_cast<int>((fine_hi2 - fine_lo2) / kFineStep) + 1;
		for (int i2 = 0; i2 < fine_steps; ++i2) {
			double d2 = fine_lo2 + i2 * kFineStep;
			for (int i3 = 0; i3 < fine_steps; ++i3) {
				double d3 = fine_lo3 + i3 * kFineStep;
				double avg_residual = ComputeAvgResidual(ev_z, ev_pos, count, d2, d3);
				double reg = lambda * (d2 * d2 + d3 * d3);
				double total = avg_residual + reg;
				if (total < best_residual) {
					best_residual = total;
					best_d2 = d2;
					best_d3 = d3;
				}
			}
		}
		best_residual = ComputeAvgResidual(ev_z, ev_pos, count, best_d2, best_d3);
		std::cout << "  Fine best:   d2=" << best_d2 << "  d3=" << best_d3
			<< "  residual=" << best_residual << "\n";

		// Stage 3: parabolic interpolation for sub-grid accuracy
		// R(d2,d3) = A*d2^2 + B*d3^2 + C*d2*d3 + D*d2 + E*d3 + F
		double h = kFineStep;
		double R3[3][3];
		for (int i = 0; i < 3; ++i) {
			for (int j = 0; j < 3; ++j) {
				double d2 = best_d2 + (i - 1) * h;
				double d3 = best_d3 + (j - 1) * h;
				double reg = lambda * (d2 * d2 + d3 * d3);
				R3[i][j] = ComputeAvgResidual(ev_z, ev_pos, count, d2, d3) + reg;
			}
		}
		double D = (R3[2][1] - R3[0][1]) / (2 * h);
		double E = (R3[1][2] - R3[1][0]) / (2 * h);
		double A = (R3[2][1] - 2 * R3[1][1] + R3[0][1]) / (2 * h * h);
		double B = (R3[1][2] - 2 * R3[1][1] + R3[1][0]) / (2 * h * h);
		double C = (R3[2][2] - R3[2][0] - R3[0][2] + R3[0][0]) / (4 * h * h);

		double det = 4 * A * B - C * C;
		if (std::fabs(det) > 1e-12 && A > 0 && B > 0) {
			double d2_corr = (C * E - 2 * B * D) / det;
			double d3_corr = (C * D - 2 * A * E) / det;
			if (std::fabs(d2_corr) < h && std::fabs(d3_corr) < h) {
				best_d2 += d2_corr;
				best_d3 += d3_corr;
				double var_r;
				best_residual = ComputeAvgResidual(ev_z, ev_pos, count, best_d2, best_d3, &var_r);
				std::cout << "  Interp best: d2=" << best_d2 << "  d3=" << best_d3
					<< "  residual=" << best_residual << "\n";
			}
		}

		// Stage 4: compute Hessian with larger step size for reliable uncertainty
		double h_hess = kCoarseStep;
		double H3[3][3];
		for (int i = 0; i < 3; ++i) {
			for (int j = 0; j < 3; ++j) {
				double d2 = best_d2 + (i - 1) * h_hess;
				double d3 = best_d3 + (j - 1) * h_hess;
				double reg = lambda * (d2 * d2 + d3 * d3);
				H3[i][j] = ComputeAvgResidual(ev_z, ev_pos, count, d2, d3) + reg;
			}
		}
		double Ah = (H3[2][1] - 2 * H3[1][1] + H3[0][1]) / (2 * h_hess * h_hess);
		double Bh = (H3[1][2] - 2 * H3[1][1] + H3[1][0]) / (2 * h_hess * h_hess);
		double Ch = (H3[2][2] - H3[2][0] - H3[0][2] + H3[0][0]) / (4 * h_hess * h_hess);
		double det_h = 4 * Ah * Bh - Ch * Ch;

		improvement = (baseline_residual - best_residual) / baseline_residual * 100;

		if (det_h > 1e-12 && Ah > 0 && Bh > 0) {
			double scale = best_residual / count;
			dsigma2 = std::sqrt(scale * 2 * Bh / det_h);
			dsigma3 = std::sqrt(scale * 2 * Ah / det_h);
			std::cout << "  Uncertainty: sigma_d2=" << dsigma2 << "  sigma_d3=" << dsigma3
				<< "  improvement=" << improvement << "%\n";
		} else {
			dsigma2 = -1;
			dsigma3 = -1;
			std::cout << "  Uncertainty: degeneracy broken by regularization\n";
			std::cout << "  improvement=" << improvement << "%\n";
		}
	};

	ScanDirection(x_events_z, x_events_pos, x_count,
		result.baseline_residual_x, result.dx2, result.dx3, result.best_residual_x,
		result.dsigma_x2, result.dsigma_x3, result.improvement_x);

	ScanDirection(y_events_z, y_events_pos, y_count,
		result.baseline_residual_y, result.dy2, result.dy3, result.best_residual_y,
		result.dsigma_y2, result.dsigma_y3, result.improvement_y);

	f.Close();
	return 0;
}

} // namespace brill