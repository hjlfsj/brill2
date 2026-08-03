#include "include/event/t0/adjust_t0.h"

#include <cmath>
#include <iostream>

#include <TFile.h>
#include <TFitResult.h>
#include <TH1D.h>
#include <TROOT.h>
#include <TString.h>
#include <TTree.h>

#include "include/event/ppac/ppac_track_event.h"
#include "include/event/t0/dssd_match_event.h"
#include "include/utils.h"

namespace brill {

int AdjustT0Step1(
	const T0AdjustConfig &adjust_config,
	const std::vector<SquareDetectorConfig> &detectors,
	std::vector<T0AdjustResult> &results
) {
	results.clear();
	int ndet = (int)detectors.size();

	TFile track_file(adjust_config.track_path.c_str(), "read");
	TTree *track_tree = (TTree*)track_file.Get("tree");
	if (!track_tree) {
		std::cerr << "Error: Get tree from " << adjust_config.track_path << " failed.\n";
		return 1;
	}
	PpacTrackEvent track_event;
	SetupInput(track_tree, track_event);

	std::vector<TFile*> match_files(ndet, nullptr);
	std::vector<TTree*> match_trees(ndet, nullptr);
	std::vector<DssdMatchEvent> match_events(ndet);

	for (int i = 0; i < ndet; ++i) {
		TString match_path = TString::Format(
			"%s/%s_%s%04d.root",
			adjust_config.match_dir.c_str(),
			detectors[i].name.c_str(),
			TriggerInfix(adjust_config.trigger).c_str(),
			adjust_config.run
		);
		match_files[i] = new TFile(match_path, "read");
		match_trees[i] = (TTree*)match_files[i]->Get("tree");
		if (!match_trees[i]) {
			std::cerr << "Error: Get tree from " << match_path << " failed.\n";
			return 1;
		}
		SetupInput(match_trees[i], match_events[i]);
	}

	TH1D *h_res_x[3] = {nullptr, nullptr, nullptr};
	TH1D *h_res_y[3] = {nullptr, nullptr, nullptr};
	TH1D *h_energy[3] = {nullptr, nullptr, nullptr};
	for (int i = 0; i < ndet; ++i) {
		h_res_x[i] = new TH1D(
			TString::Format("h_res_x_%s", detectors[i].name.c_str()),
			TString::Format("%s X residual;residual_x (mm);Counts", detectors[i].name.c_str()),
			200, -10.0, 10.0);
		h_res_y[i] = new TH1D(
			TString::Format("h_res_y_%s", detectors[i].name.c_str()),
			TString::Format("%s Y residual;residual_y (mm);Counts", detectors[i].name.c_str()),
			200, -10.0, 10.0);
		h_energy[i] = new TH1D(
			TString::Format("h_energy_%s", detectors[i].name.c_str()),
			TString::Format("%s energy;ADC;Counts", detectors[i].name.c_str()),
			1000, 0.0, 65536.0);
		h_res_x[i]->SetDirectory(nullptr);
		h_res_y[i]->SetDirectory(nullptr);
		h_energy[i]->SetDirectory(nullptr);
	}

	long long total = match_trees[0]->GetEntries();
	long long n_track_valid = 0;
	long long n_all_num1 = 0;
	long long n_energy_pass = 0;
	int num_selected = 0;
	int sample_count = 0;

	for (long long entry = 0; entry < total; ++entry) {
		track_tree->GetEntry(entry);
		if (!track_event.valid) continue;
		++n_track_valid;

		bool all_ok = true;
		for (int i = 0; i < ndet; ++i) {
			match_trees[i]->GetEntry(entry);
			if (match_events[i].num != 1) { all_ok = false; break; }
		}
		if (!all_ok) continue;
		++n_all_num1;

		if (sample_count < 20) {
			std::cout << "  Entry " << entry << " energies:";
			for (int i = 0; i < ndet; ++i)
				std::cout << " " << detectors[i].name << "=" << match_events[i].energy[0];
			std::cout << "\n";
			++sample_count;
		}

		for (int i = 0; i < ndet; ++i) {
			h_energy[i]->Fill(match_events[i].energy[0]);
		}

		for (int i = 0; i < ndet; ++i) {
			const auto &cut = adjust_config.energy_cuts[i];
			if (cut.max > 0.0) {
				double e = match_events[i].energy[0];
				if (e <= cut.min || e >= cut.max) { all_ok = false; break; }
			}
		}
		if (!all_ok) continue;
		++n_energy_pass;

		for (int i = 0; i < ndet; ++i) {
			const auto &det = detectors[i];
			double z = det.z_mm;
			double predicted_x = track_event.target_x + track_event.dir_x * z;
			double predicted_y = track_event.target_y + track_event.dir_y * z;

			double match_x = match_events[i].x[0];
			double match_y = match_events[i].y[0];

			h_res_x[i]->Fill(match_x - predicted_x);
			h_res_y[i]->Fill(match_y - predicted_y);
		}
		++num_selected;
	}

	std::cout << "Total events:    " << total << "\n";
	std::cout << "Track valid:     " << n_track_valid << "\n";
	std::cout << "All num==1:      " << n_all_num1 << "\n";
	std::cout << "Energy pass:     " << n_energy_pass << "\n";
	std::cout << "Selected events: " << num_selected << "\n";

	if (num_selected < 10) {
		std::cerr << "Error: Too few events (" << num_selected << ") for fitting.\n";
		return 1;
	}

	results.resize(ndet);
	for (int i = 0; i < ndet; ++i) {
		results[i].detector_name = detectors[i].name;
		results[i].num_events = num_selected;

		double hx_entries = h_res_x[i]->GetEntries();
		double hy_entries = h_res_y[i]->GetEntries();
		std::cout << detectors[i].name << " h_res_x entries: " << hx_entries
			<< "  h_res_y entries: " << hy_entries << "\n";

		TFitResultPtr fit_x = h_res_x[i]->Fit("gaus", "QNS");
		if (!fit_x.Get() || !fit_x->IsValid()) {
			fit_x = h_res_x[i]->Fit("gaus", "QNS", "", h_res_x[i]->GetMean() - h_res_x[i]->GetRMS(),
				h_res_x[i]->GetMean() + h_res_x[i]->GetRMS());
		}

		TFitResultPtr fit_y = h_res_y[i]->Fit("gaus", "QNS");
		if (!fit_y.Get() || !fit_y->IsValid()) {
			fit_y = h_res_y[i]->Fit("gaus", "QNS", "", h_res_y[i]->GetMean() - h_res_y[i]->GetRMS(),
				h_res_y[i]->GetMean() + h_res_y[i]->GetRMS());
		}

		if (fit_x.Get() && fit_x->IsValid()) {
			results[i].dx = -fit_x->Parameter(1);
			results[i].dsigma_x = fit_x->Error(1);
			results[i].residual_x_mean = fit_x->Parameter(1);
			results[i].residual_x_sigma = fit_x->Parameter(2);
		} else {
			std::cerr << "Warning: X fit failed for " << detectors[i].name << ".\n";
			results[i].dx = 0.0;
			results[i].dsigma_x = -1.0;
		}

		if (fit_y.Get() && fit_y->IsValid()) {
			results[i].dy = -fit_y->Parameter(1);
			results[i].dsigma_y = fit_y->Error(1);
			results[i].residual_y_mean = fit_y->Parameter(1);
			results[i].residual_y_sigma = fit_y->Parameter(2);
		} else {
			std::cerr << "Warning: Y fit failed for " << detectors[i].name << ".\n";
			results[i].dy = 0.0;
			results[i].dsigma_y = -1.0;
		}
	}

	TString output_path = TString::Format(
		"%s/t0_offset_%s%04d.root",
		adjust_config.output_dir.c_str(),
		TriggerInfix(adjust_config.trigger).c_str(),
		adjust_config.run
	);

	{
		TFile opf(output_path, "recreate");
		for (int i = 0; i < ndet; ++i) {
			h_res_x[i]->Write();
			h_res_y[i]->Write();
			h_energy[i]->Write();
		}
		opf.Close();
	}
	gROOT->cd();

	std::cout << "Residual histograms saved to " << output_path << "\n";

	return 0;
}

int AdjustT0Step2(
	const T0AdjustConfig &adjust_config,
	const std::vector<SquareDetectorConfig> &detectors,
	std::vector<T0AdjustResult> &results
) {
	if (detectors.size() != 4) {
		std::cerr << "Error: AdjustT0Step2 requires exactly 4 detectors (d1-d4).\n";
		return 1;
	}
	results.clear();

	std::vector<TFile*> match_files(4, nullptr);
	std::vector<TTree*> match_trees(4, nullptr);
	std::vector<DssdMatchEvent> match_events(4);

	for (int i = 0; i < 4; ++i) {
		TString match_path = TString::Format(
			"%s/%s_%s%04d.root",
			adjust_config.match_dir.c_str(),
			detectors[i].name.c_str(),
			TriggerInfix(adjust_config.trigger).c_str(),
			adjust_config.run
		);
		match_files[i] = new TFile(match_path, "read");
		match_trees[i] = (TTree*)match_files[i]->Get("tree");
		if (!match_trees[i]) {
			std::cerr << "Error: Get tree from " << match_path << " failed.\n";
			return 1;
		}
		SetupInput(match_trees[i], match_events[i]);
	}

	TH1D *h_res_x = new TH1D(
		"h_res_x_t0d4", "T0D4 X residual;residual_x (mm);Counts",
		20, -10.0, 10.0);
	TH1D *h_res_y = new TH1D(
		"h_res_y_t0d4", "T0D4 Y residual;residual_y (mm);Counts",
		20, -10.0, 10.0);
	h_res_x->SetDirectory(nullptr);
	h_res_y->SetDirectory(nullptr);

	const int kMaxEvents = 10000;
	long long total = match_trees[0]->GetEntries();
	int num_selected = 0;

	double z1 = detectors[0].z_mm;
	double z2 = detectors[1].z_mm;
	double z3 = detectors[2].z_mm;
	double z4 = detectors[3].z_mm;

	double mean_z = (z1 + z2 + z3) / 3.0;
	double szz = (z1 - mean_z) * (z1 - mean_z)
		+ (z2 - mean_z) * (z2 - mean_z)
		+ (z3 - mean_z) * (z3 - mean_z);

	for (long long entry = 0; entry < total && num_selected < kMaxEvents; ++entry) {
		bool all_ok = true;
		for (int i = 0; i < 4; ++i) {
			match_trees[i]->GetEntry(entry);
			if (match_events[i].num != 1) { all_ok = false; break; }
		}
		if (!all_ok) continue;

		double x1 = match_events[0].x[0];
		double x2 = match_events[1].x[0];
		double x3 = match_events[2].x[0];
		double x4 = match_events[3].x[0];
		double mean_x = (x1 + x2 + x3) / 3.0;
		double szx = (z1 - mean_z) * (x1 - mean_x)
			+ (z2 - mean_z) * (x2 - mean_x)
			+ (z3 - mean_z) * (x3 - mean_x);
		double b_x = szx / szz;
		double a_x = mean_x - b_x * mean_z;
		double predicted_x = a_x + b_x * z4;
		h_res_x->Fill(x4 - predicted_x);

		double y1 = match_events[0].y[0];
		double y2 = match_events[1].y[0];
		double y3 = match_events[2].y[0];
		double y4 = match_events[3].y[0];
		double mean_y = (y1 + y2 + y3) / 3.0;
		double szy = (z1 - mean_z) * (y1 - mean_y)
			+ (z2 - mean_z) * (y2 - mean_y)
			+ (z3 - mean_z) * (y3 - mean_y);
		double b_y = szy / szz;
		double a_y = mean_y - b_y * mean_z;
		double predicted_y = a_y + b_y * z4;
		h_res_y->Fill(y4 - predicted_y);

		++num_selected;
	}

	std::cout << "Total events:    " << total << "\n";
	std::cout << "Selected events: " << num_selected << "\n";

	if (num_selected < 10) {
		std::cerr << "Error: Too few events (" << num_selected << ") for fitting.\n";
		return 1;
	}

	T0AdjustResult result;
	result.detector_name = "t0d4";
	result.num_events = num_selected;

	TFitResultPtr fit_x = h_res_x->Fit("gaus", "QNS");
	if (!fit_x.Get() || !fit_x->IsValid()) {
		fit_x = h_res_x->Fit("gaus", "QNS", "",
			h_res_x->GetMean() - h_res_x->GetRMS(),
			h_res_x->GetMean() + h_res_x->GetRMS());
	}
	if (fit_x.Get() && fit_x->IsValid()) {
		result.dx = -fit_x->Parameter(1);
		result.dsigma_x = fit_x->Error(1);
		result.residual_x_mean = fit_x->Parameter(1);
		result.residual_x_sigma = fit_x->Parameter(2);
	} else {
		std::cerr << "Warning: X fit failed for t0d4.\n";
		result.dx = 0.0;
		result.dsigma_x = -1.0;
	}

	TFitResultPtr fit_y = h_res_y->Fit("gaus", "QNS");
	if (!fit_y.Get() || !fit_y->IsValid()) {
		fit_y = h_res_y->Fit("gaus", "QNS", "",
			h_res_y->GetMean() - h_res_y->GetRMS(),
			h_res_y->GetMean() + h_res_y->GetRMS());
	}
	if (fit_y.Get() && fit_y->IsValid()) {
		result.dy = -fit_y->Parameter(1);
		result.dsigma_y = fit_y->Error(1);
		result.residual_y_mean = fit_y->Parameter(1);
		result.residual_y_sigma = fit_y->Parameter(2);
	} else {
		std::cerr << "Warning: Y fit failed for t0d4.\n";
		result.dy = 0.0;
		result.dsigma_y = -1.0;
	}

	results.push_back(result);

	TString output_path = TString::Format(
		"%s/t0_offset_d4_%s%04d.root",
		adjust_config.output_dir.c_str(),
		TriggerInfix(adjust_config.trigger).c_str(),
		adjust_config.run
	);

	{
		TFile opf(output_path, "recreate");
		h_res_x->Write();
		h_res_y->Write();
		opf.Close();
	}
	gROOT->cd();

	std::cout << "Residual histograms saved to " << output_path << "\n";

	return 0;
}

} // namespace brill