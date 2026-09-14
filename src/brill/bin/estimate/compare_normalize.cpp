#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMultiGraph.h>
#include <TString.h>

#include "external/cxxopts.hpp"
#include "include/config.h"
#include "include/t0/dssd.h"
#include "include/utils.h"

void PrintUsage(const cxxopts::Options &options) {
	std::cout << options.help() << "\n";
}

int main(int argc, char **argv) {
	cxxopts::Options options(
		"compare_normalize",
		"Compare DSSD normalize parameters across different runs."
	);
	options.add_options()
		("h,help", "Print help information.")
		(
			"c,config",
			"Config file path.",
			cxxopts::value<std::string>()->default_value("config.toml"),
			"file"
		)
		(
			"t,trigger",
			"Trigger type (default: use config value).",
			cxxopts::value<std::string>(),
			"trigger"
		);

	auto result = options.parse(argc, argv);
	if (result.count("help")) {
		PrintUsage(options);
		return 0;
	}

	brill::AppConfig config;
	if (brill::LoadConfig(result["config"].as<std::string>(), config)) {
		return 1;
	}
	if (result.count("trigger")) {
		config.trigger = result["trigger"].as<std::string>();
	}

	if (config.normalize.runs.empty()) {
		std::cerr << "Error: No normalize.runs configured in config.toml.\n";
		return 1;
	}

	const std::vector<std::string> detectors = {"t0d1", "t0d2", "t0d3", "t0d4"};
	const int d_offset[] = {0, 32, 96, 128};
	const int colors[] = {kBlack, kRed, kBlue, kGreen+2, kMagenta, kCyan, kOrange+1, kViolet};
	const int n_colors = sizeof(colors) / sizeof(colors[0]);

	std::string normalize_dir = brill::JoinPath(config.workspace, config.paths.normalize);
	std::string output_dir = brill::JoinPath(config.workspace, config.paths.estimate);
	TString output_path = TString::Format(
		"%s/compare_normalize.root",
		output_dir.c_str()
	);
	TFile opf(output_path, "recreate");

	TMultiGraph mg_front_p0;
	TMultiGraph mg_front_p1;
	TMultiGraph mg_back_p0;
	TMultiGraph mg_back_p1;
	TMultiGraph mg_front_fwhm;
	TMultiGraph mg_front_mean;
	TMultiGraph mg_back_fwhm;
	TMultiGraph mg_back_mean;

	mg_front_p0.SetTitle("Front p0 (offset);Global strip;Offset");
	mg_front_p1.SetTitle("Front p1 (scale);Global strip;Scale");
	mg_back_p0.SetTitle("Back p0 (offset);Global strip;Offset");
	mg_back_p1.SetTitle("Back p1 (scale);Global strip;Scale");
	mg_front_fwhm.SetTitle("Front FWHM;Global strip;FWHM");
	mg_front_mean.SetTitle("Front Peak;Global strip;Peak");
	mg_back_fwhm.SetTitle("Back FWHM;Global strip;FWHM");
	mg_back_mean.SetTitle("Back Peak;Global strip;Peak");

	TLegend leg_front_p0(0.82, 0.78, 0.95, 0.95);
	TLegend leg_front_p1(0.82, 0.78, 0.95, 0.95);
	TLegend leg_back_p0(0.82, 0.78, 0.95, 0.95);
	TLegend leg_back_p1(0.82, 0.78, 0.95, 0.95);
	TLegend leg_front_fwhm(0.82, 0.78, 0.95, 0.95);
	TLegend leg_front_mean(0.82, 0.78, 0.95, 0.95);
	TLegend leg_back_fwhm(0.82, 0.78, 0.95, 0.95);
	TLegend leg_back_mean(0.82, 0.78, 0.95, 0.95);

	for (size_t det_idx = 0; det_idx < detectors.size(); ++det_idx) {
		const std::string &detector_name = detectors[det_idx];
		const brill::SquareDetectorConfig *detector =
			brill::FindDetectorConfig(config, detector_name);
		if (!detector) {
			std::cerr << "Warning: Detector " << detector_name
				<< " not found in config, skipping.\n";
			continue;
		}

		const int front_strips = detector->front_strips;
		const int back_strips = detector->back_strips;
		const int offset = d_offset[det_idx];

		for (size_t run_idx = 0; run_idx < config.normalize.runs.size(); ++run_idx) {
			int normalize_run = config.normalize.runs[run_idx].second;

			TString front_path = TString::Format(
				"%s/%s_front_%s%04d.txt",
				normalize_dir.c_str(),
				detector_name.c_str(),
				brill::TriggerInfix(config.trigger).c_str(),
				normalize_run
			);
			TString back_path = TString::Format(
				"%s/%s_back_%s%04d.txt",
				normalize_dir.c_str(),
				detector_name.c_str(),
				brill::TriggerInfix(config.trigger).c_str(),
				normalize_run
			);

			brill::DssdNormalizeParameters parameters;
			parameters.front_strips = front_strips;
			parameters.back_strips = back_strips;

			std::cout << "Reading " << detector_name << " run " << normalize_run
				<< ":\n  " << front_path << "\n  " << back_path << "\n";

			if (brill::ReadDssdNormalizeParameters(
				front_path.Data(), back_path.Data(), parameters
			)) {
				std::cerr << "Warning: Failed to read parameters for "
					<< detector_name << " run " << normalize_run << ", skipping.\n";
				continue;
			}

			int run_color = colors[run_idx % n_colors];
			TString run_label = TString::Format("run %d", normalize_run);

			TGraph *g_front_p0 = new TGraph(front_strips);
			TGraph *g_front_p1 = new TGraph(front_strips);
			TGraph *g_back_p0 = new TGraph(back_strips);
			TGraph *g_back_p1 = new TGraph(back_strips);
			TGraph *g_front_fwhm = new TGraph();
			TGraph *g_front_mean = new TGraph();
			TGraph *g_back_fwhm = new TGraph();
			TGraph *g_back_mean = new TGraph();

			g_front_p0->SetName(TString::Format("g_%s_fp0_r%d", detector_name.c_str(), normalize_run));
			g_front_p1->SetName(TString::Format("g_%s_fp1_r%d", detector_name.c_str(), normalize_run));
			g_back_p0->SetName(TString::Format("g_%s_bp0_r%d", detector_name.c_str(), normalize_run));
			g_back_p1->SetName(TString::Format("g_%s_bp1_r%d", detector_name.c_str(), normalize_run));
			g_front_fwhm->SetName(TString::Format("g_%s_ffwhm_r%d", detector_name.c_str(), normalize_run));
			g_front_mean->SetName(TString::Format("g_%s_fmean_r%d", detector_name.c_str(), normalize_run));
			g_back_fwhm->SetName(TString::Format("g_%s_bfwhm_r%d", detector_name.c_str(), normalize_run));
			g_back_mean->SetName(TString::Format("g_%s_bmean_r%d", detector_name.c_str(), normalize_run));

			TString full_label = TString::Format("%s %s", detector_name.c_str(), run_label.Data());
			g_front_p0->SetTitle(full_label);
			g_front_p1->SetTitle(full_label);
			g_back_p0->SetTitle(full_label);
			g_back_p1->SetTitle(full_label);
			g_front_fwhm->SetTitle(full_label);
			g_front_mean->SetTitle(full_label);
			g_back_fwhm->SetTitle(full_label);
			g_back_mean->SetTitle(full_label);

			for (auto *g : {g_front_p0, g_front_p1, g_back_p0, g_back_p1,
			                g_front_fwhm, g_front_mean, g_back_fwhm, g_back_mean}) {
				g->SetLineColor(run_color);
				g->SetMarkerColor(run_color);
				g->SetMarkerStyle(20 + int(run_idx));
			}

			for (int i = 0; i < front_strips; ++i) {
				g_front_p0->SetPoint(i, offset + i, parameters.front_p0[i]);
				g_front_p1->SetPoint(i, offset + i, parameters.front_p1[i]);
			}
			for (int i = 0; i < back_strips; ++i) {
				g_back_p0->SetPoint(i, offset + i, parameters.back_p0[i]);
				g_back_p1->SetPoint(i, offset + i, parameters.back_p1[i]);
			}

			int end_run = normalize_run + 19;
			if (normalize_run == 97) end_run = 115;

			TString root_path = TString::Format(
				"%s/%s_%s%04d_%04d.root",
				normalize_dir.c_str(),
				detector_name.c_str(),
				brill::TriggerInfix(config.trigger).c_str(),
				normalize_run,
				end_run
			);

			std::cout << "  " << root_path << "\n";
			TFile rf(root_path, "read");
			if (rf.IsZombie()) {
				std::cerr << "Warning: Cannot open ROOT file " << root_path
					<< ", skipping residual analysis.\n";
			} else {
				int front_failed = 0, back_failed = 0;

				for (int i = 0; i < front_strips; ++i) {
					TString hname = TString::Format("resf%d", i);
					TH1D *h = dynamic_cast<TH1D*>(rf.Get(hname));
					if (!h) {
						front_failed++;
						continue;
					}
					TF1 gaus(TString::Format("gf%d_%zu_%zu", i, det_idx, run_idx), "gaus", -500, 500);
					h->Fit(&gaus, "QR");
					double mean = gaus.GetParameter(1);
					double sigma = gaus.GetParameter(2);
					if (sigma > 0 && sigma < 1000) {
						g_front_fwhm->SetPoint(g_front_fwhm->GetN(), offset + i, 2.355 * sigma);
						g_front_mean->SetPoint(g_front_mean->GetN(), offset + i, mean);
					} else {
						front_failed++;
					}
				}

				for (int i = 0; i < back_strips; ++i) {
					TString hname = TString::Format("resb%d", i);
					TH1D *h = dynamic_cast<TH1D*>(rf.Get(hname));
					if (!h) {
						back_failed++;
						continue;
					}
					TF1 gaus(TString::Format("gb%d_%zu_%zu", i, det_idx, run_idx), "gaus", -500, 500);
					h->Fit(&gaus, "QR");
					double mean = gaus.GetParameter(1);
					double sigma = gaus.GetParameter(2);
					if (sigma > 0 && sigma < 1000) {
						g_back_fwhm->SetPoint(g_back_fwhm->GetN(), offset + i, 2.355 * sigma);
						g_back_mean->SetPoint(g_back_mean->GetN(), offset + i, mean);
					} else {
						back_failed++;
					}
				}

				rf.Close();
				opf.cd();

				if (front_failed > 0) {
					std::cerr << "Warning: " << detector_name << " run " << normalize_run
						<< " front: " << front_failed << "/" << front_strips
						<< " strips failed Gaussian fit.\n";
				}
				if (back_failed > 0) {
					std::cerr << "Warning: " << detector_name << " run " << normalize_run
						<< " back: " << back_failed << "/" << back_strips
						<< " strips failed Gaussian fit.\n";
				}
			}

			mg_front_p0.Add(g_front_p0, "lp");
			mg_front_p1.Add(g_front_p1, "lp");
			mg_back_p0.Add(g_back_p0, "lp");
			mg_back_p1.Add(g_back_p1, "lp");
			mg_front_fwhm.Add(g_front_fwhm, "lp");
			mg_front_mean.Add(g_front_mean, "lp");
			mg_back_fwhm.Add(g_back_fwhm, "lp");
			mg_back_mean.Add(g_back_mean, "lp");

			if (det_idx == 0) {
				leg_front_p0.AddEntry(g_front_p0, run_label, "lp");
				leg_front_p1.AddEntry(g_front_p1, run_label, "lp");
				leg_back_p0.AddEntry(g_back_p0, run_label, "lp");
				leg_back_p1.AddEntry(g_back_p1, run_label, "lp");
				leg_front_fwhm.AddEntry(g_front_fwhm, run_label, "lp");
				leg_front_mean.AddEntry(g_front_mean, run_label, "lp");
				leg_back_fwhm.AddEntry(g_back_fwhm, run_label, "lp");
				leg_back_mean.AddEntry(g_back_mean, run_label, "lp");
			}
		}
	}

	TCanvas c_front_p0("c_front_p0", "Front p0 (offset)", 1600, 800);
	mg_front_p0.Draw("a lp");
	double ymin = mg_front_p0.GetHistogram()->GetMinimum();
	double ymax = mg_front_p0.GetHistogram()->GetMaximum();
	TLine vl1(31.5, ymin, 31.5, ymax);
	TLine vl2(95.5, ymin, 95.5, ymax);
	TLine vl3(127.5, ymin, 127.5, ymax);
	vl1.SetLineStyle(2); vl2.SetLineStyle(2); vl3.SetLineStyle(2);
	vl1.SetLineColor(kGray+2); vl2.SetLineColor(kGray+2); vl3.SetLineColor(kGray+2);
	vl1.Draw(); vl2.Draw(); vl3.Draw();
	TLatex().DrawLatex(15, ymin, "t0d1");
	TLatex().DrawLatex(63, ymin, "t0d2");
	TLatex().DrawLatex(111, ymin, "t0d3");
	TLatex().DrawLatex(143, ymin, "t0d4");
	leg_front_p0.Draw();
	c_front_p0.Write();

	TCanvas c_front_p1("c_front_p1", "Front p1 (scale)", 1600, 800);
	mg_front_p1.Draw("a lp");
	ymin = mg_front_p1.GetHistogram()->GetMinimum();
	ymax = mg_front_p1.GetHistogram()->GetMaximum();
	TLine vl1b(31.5, ymin, 31.5, ymax);
	TLine vl2b(95.5, ymin, 95.5, ymax);
	TLine vl3b(127.5, ymin, 127.5, ymax);
	vl1b.SetLineStyle(2); vl2b.SetLineStyle(2); vl3b.SetLineStyle(2);
	vl1b.SetLineColor(kGray+2); vl2b.SetLineColor(kGray+2); vl3b.SetLineColor(kGray+2);
	vl1b.Draw(); vl2b.Draw(); vl3b.Draw();
	TLatex().DrawLatex(15, ymin, "t0d1");
	TLatex().DrawLatex(63, ymin, "t0d2");
	TLatex().DrawLatex(111, ymin, "t0d3");
	TLatex().DrawLatex(143, ymin, "t0d4");
	leg_front_p1.Draw();
	c_front_p1.Write();

	TCanvas c_back_p0("c_back_p0", "Back p0 (offset)", 1600, 800);
	mg_back_p0.Draw("a lp");
	ymin = mg_back_p0.GetHistogram()->GetMinimum();
	ymax = mg_back_p0.GetHistogram()->GetMaximum();
	TLine vl1c(31.5, ymin, 31.5, ymax);
	TLine vl2c(95.5, ymin, 95.5, ymax);
	TLine vl3c(127.5, ymin, 127.5, ymax);
	vl1c.SetLineStyle(2); vl2c.SetLineStyle(2); vl3c.SetLineStyle(2);
	vl1c.SetLineColor(kGray+2); vl2c.SetLineColor(kGray+2); vl3c.SetLineColor(kGray+2);
	vl1c.Draw(); vl2c.Draw(); vl3c.Draw();
	TLatex().DrawLatex(15, ymin, "t0d1");
	TLatex().DrawLatex(63, ymin, "t0d2");
	TLatex().DrawLatex(111, ymin, "t0d3");
	TLatex().DrawLatex(143, ymin, "t0d4");
	leg_back_p0.Draw();
	c_back_p0.Write();

	TCanvas c_back_p1("c_back_p1", "Back p1 (scale)", 1600, 800);
	mg_back_p1.Draw("a lp");
	ymin = mg_back_p1.GetHistogram()->GetMinimum();
	ymax = mg_back_p1.GetHistogram()->GetMaximum();
	TLine vl1d(31.5, ymin, 31.5, ymax);
	TLine vl2d(95.5, ymin, 95.5, ymax);
	TLine vl3d(127.5, ymin, 127.5, ymax);
	vl1d.SetLineStyle(2); vl2d.SetLineStyle(2); vl3d.SetLineStyle(2);
	vl1d.SetLineColor(kGray+2); vl2d.SetLineColor(kGray+2); vl3d.SetLineColor(kGray+2);
	vl1d.Draw(); vl2d.Draw(); vl3d.Draw();
	TLatex().DrawLatex(15, ymin, "t0d1");
	TLatex().DrawLatex(63, ymin, "t0d2");
	TLatex().DrawLatex(111, ymin, "t0d3");
	TLatex().DrawLatex(143, ymin, "t0d4");
	leg_back_p1.Draw();
	c_back_p1.Write();

	TCanvas c_front_fwhm("c_front_fwhm", "Front FWHM", 1600, 800);
	mg_front_fwhm.Draw("a lp");
	ymin = mg_front_fwhm.GetHistogram()->GetMinimum();
	ymax = mg_front_fwhm.GetHistogram()->GetMaximum();
	TLine vf1(31.5, ymin, 31.5, ymax);
	TLine vf2(95.5, ymin, 95.5, ymax);
	TLine vf3(127.5, ymin, 127.5, ymax);
	vf1.SetLineStyle(2); vf2.SetLineStyle(2); vf3.SetLineStyle(2);
	vf1.SetLineColor(kGray+2); vf2.SetLineColor(kGray+2); vf3.SetLineColor(kGray+2);
	vf1.Draw(); vf2.Draw(); vf3.Draw();
	TLatex().DrawLatex(15, ymin, "t0d1");
	TLatex().DrawLatex(63, ymin, "t0d2");
	TLatex().DrawLatex(111, ymin, "t0d3");
	TLatex().DrawLatex(143, ymin, "t0d4");
	leg_front_fwhm.Draw();
	c_front_fwhm.Write();

	TCanvas c_front_mean("c_front_mean", "Front Peak", 1600, 800);
	mg_front_mean.Draw("a lp");
	ymin = mg_front_mean.GetHistogram()->GetMinimum();
	ymax = mg_front_mean.GetHistogram()->GetMaximum();
	TLine vm1(31.5, ymin, 31.5, ymax);
	TLine vm2(95.5, ymin, 95.5, ymax);
	TLine vm3(127.5, ymin, 127.5, ymax);
	vm1.SetLineStyle(2); vm2.SetLineStyle(2); vm3.SetLineStyle(2);
	vm1.SetLineColor(kGray+2); vm2.SetLineColor(kGray+2); vm3.SetLineColor(kGray+2);
	vm1.Draw(); vm2.Draw(); vm3.Draw();
	TLatex().DrawLatex(15, ymin, "t0d1");
	TLatex().DrawLatex(63, ymin, "t0d2");
	TLatex().DrawLatex(111, ymin, "t0d3");
	TLatex().DrawLatex(143, ymin, "t0d4");
	leg_front_mean.Draw();
	c_front_mean.Write();

	TCanvas c_back_fwhm("c_back_fwhm", "Back FWHM", 1600, 800);
	mg_back_fwhm.Draw("a lp");
	ymin = mg_back_fwhm.GetHistogram()->GetMinimum();
	ymax = mg_back_fwhm.GetHistogram()->GetMaximum();
	TLine vb1(31.5, ymin, 31.5, ymax);
	TLine vb2(95.5, ymin, 95.5, ymax);
	TLine vb3(127.5, ymin, 127.5, ymax);
	vb1.SetLineStyle(2); vb2.SetLineStyle(2); vb3.SetLineStyle(2);
	vb1.SetLineColor(kGray+2); vb2.SetLineColor(kGray+2); vb3.SetLineColor(kGray+2);
	vb1.Draw(); vb2.Draw(); vb3.Draw();
	TLatex().DrawLatex(15, ymin, "t0d1");
	TLatex().DrawLatex(63, ymin, "t0d2");
	TLatex().DrawLatex(111, ymin, "t0d3");
	TLatex().DrawLatex(143, ymin, "t0d4");
	leg_back_fwhm.Draw();
	c_back_fwhm.Write();

	TCanvas c_back_mean("c_back_mean", "Back Peak", 1600, 800);
	mg_back_mean.Draw("a lp");
	ymin = mg_back_mean.GetHistogram()->GetMinimum();
	ymax = mg_back_mean.GetHistogram()->GetMaximum();
	TLine vn1(31.5, ymin, 31.5, ymax);
	TLine vn2(95.5, ymin, 95.5, ymax);
	TLine vn3(127.5, ymin, 127.5, ymax);
	vn1.SetLineStyle(2); vn2.SetLineStyle(2); vn3.SetLineStyle(2);
	vn1.SetLineColor(kGray+2); vn2.SetLineColor(kGray+2); vn3.SetLineColor(kGray+2);
	vn1.Draw(); vn2.Draw(); vn3.Draw();
	TLatex().DrawLatex(15, ymin, "t0d1");
	TLatex().DrawLatex(63, ymin, "t0d2");
	TLatex().DrawLatex(111, ymin, "t0d3");
	TLatex().DrawLatex(143, ymin, "t0d4");
	leg_back_mean.Draw();
	c_back_mean.Write();

	opf.Close();
	std::cout << "Output written to " << output_path << "\n";
	return 0;
}