#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TAxis.h>
#include <TFile.h>
#include <TGraph.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMultiGraph.h>
#include <TString.h>

#include "external/cxxopts.hpp"
#include "include/config.h"
#include "include/utils.h"

void PrintUsage(const cxxopts::Options &options) {
	std::cout << options.help() << "\n";
}

int main(int argc, char **argv) {
	cxxopts::Options options(
		"compare_calibration",
		"Compare T0 energy calibration parameters across different runs."
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

	if (config.calibration.runs.empty()) {
		std::cerr << "Error: No calibration.runs configured in config.toml.\n";
		return 1;
	}

	const std::vector<int> &runs = config.calibration.runs;

	const int kLayerCount = 5;
	const char *layer_names[kLayerCount] = {"t0d1", "t0d2", "t0d3", "t0d4", "t0s"};

	const int colors[] = {kBlack, kRed, kBlue, kGreen+2, kMagenta, kCyan, kOrange+1, kViolet};
	const int n_colors = sizeof(colors) / sizeof(colors[0]);

	std::string calibration_dir = brill::JoinPath(config.workspace, config.paths.calibration);
	std::string output_dir = brill::JoinPath(config.workspace, config.paths.estimate);
	TString output_path = TString::Format(
		"%s/compare_calibration.root",
		output_dir.c_str()
	);
	TFile opf(output_path, "recreate");

	TLegend leg_p0(0.75, 0.78, 0.90, 0.95);
	TLegend leg_p1(0.75, 0.78, 0.90, 0.95);

	TMultiGraph mg_p0;
	TMultiGraph mg_p1;
	mg_p0.SetTitle("Calibration p0 (offset);Detector;Offset");
	mg_p1.SetTitle("Calibration p1 (scale);Detector;Scale");

	std::vector<std::vector<double>> p0_data(runs.size(), std::vector<double>(kLayerCount));
	std::vector<std::vector<double>> p1_data(runs.size(), std::vector<double>(kLayerCount));
	std::vector<bool> load_ok(runs.size(), false);

	for (size_t run_idx = 0; run_idx < runs.size(); ++run_idx) {
		int run = runs[run_idx];
		TString file_path = TString::Format(
			"%s/t0_%04d.txt",
			calibration_dir.c_str(),
			run
		);

		std::cout << "Reading " << file_path << "\n";
		std::ifstream fin(file_path.Data());
		if (!fin.good()) {
			std::cerr << "Warning: Cannot open " << file_path << ", skipping.\n";
			continue;
		}

		std::string header;
		std::getline(fin, header);

		int index = -1;
		double p0 = 0.0;
		double p1 = 1.0;
		bool found[kLayerCount] = {false};

		while (fin >> index >> p0 >> p1) {
			if (index < 0 || index >= kLayerCount) continue;
			p0_data[run_idx][index] = p0;
			p1_data[run_idx][index] = p1;
			found[index] = true;
		}

		load_ok[run_idx] = true;
		for (int i = 0; i < kLayerCount; ++i) {
			if (!found[i]) {
				std::cerr << "Warning: Layer " << i << " not found in " << file_path << "\n";
				load_ok[run_idx] = false;
			}
		}

		std::cout << "  "
			<< "d1: (" << p0_data[run_idx][0] << ", " << p1_data[run_idx][0] << ")  "
			<< "d2: (" << p0_data[run_idx][1] << ", " << p1_data[run_idx][1] << ")  "
			<< "d3: (" << p0_data[run_idx][2] << ", " << p1_data[run_idx][2] << ")  "
			<< "d4: (" << p0_data[run_idx][3] << ", " << p1_data[run_idx][3] << ")  "
			<< "s:  (" << p0_data[run_idx][4] << ", " << p1_data[run_idx][4] << ")\n";
	}

	TGraph *g_p0 = new TGraph();
	TGraph *g_p1 = new TGraph();
	g_p0->SetName("g_p0_bar");
	g_p1->SetName("g_p1_bar");

	for (size_t run_idx = 0; run_idx < runs.size(); ++run_idx) {
		if (!load_ok[run_idx]) continue;

		int run = runs[run_idx];
		int run_color = colors[run_idx % n_colors];
		TString run_label = TString::Format("run %04d", run);

		TGraph *gp0 = new TGraph(kLayerCount);
		TGraph *gp1 = new TGraph(kLayerCount);
		gp0->SetName(TString::Format("g_p0_r%d", run));
		gp1->SetName(TString::Format("g_p1_r%d", run));
		gp0->SetTitle(run_label);
		gp1->SetTitle(run_label);

		gp0->SetLineColor(run_color);
		gp0->SetMarkerColor(run_color);
		gp0->SetMarkerStyle(20 + int(run_idx));
		gp1->SetLineColor(run_color);
		gp1->SetMarkerColor(run_color);
		gp1->SetMarkerStyle(20 + int(run_idx));

		for (int layer = 0; layer < kLayerCount; ++layer) {
			gp0->SetPoint(layer, layer, p0_data[run_idx][layer]);
			gp1->SetPoint(layer, layer, p1_data[run_idx][layer]);
			g_p0->SetPoint(g_p0->GetN(), layer + run_idx * 0.1, p0_data[run_idx][layer]);
			g_p1->SetPoint(g_p1->GetN(), layer + run_idx * 0.1, p1_data[run_idx][layer]);
		}

		mg_p0.Add(gp0, "lp");
		mg_p1.Add(gp1, "lp");
		leg_p0.AddEntry(gp0, run_label, "lp");
		leg_p1.AddEntry(gp1, run_label, "lp");
	}

	TCanvas c_p0("c_p0", "Calibration p0 (offset)", 1200, 600);
	mg_p0.Draw("a lp");
	mg_p0.GetXaxis()->Set(5, -0.5, 4.5);
	mg_p0.GetXaxis()->SetNdivisions(5);
	for (int i = 0; i < kLayerCount; ++i) {
		mg_p0.GetXaxis()->SetBinLabel(i + 1, layer_names[i]);
	}
	leg_p0.Draw();
	c_p0.Write();

	TCanvas c_p1("c_p1", "Calibration p1 (scale)", 1200, 600);
	mg_p1.Draw("a lp");
	mg_p1.GetXaxis()->Set(5, -0.5, 4.5);
	mg_p1.GetXaxis()->SetNdivisions(5);
	for (int i = 0; i < kLayerCount; ++i) {
		mg_p1.GetXaxis()->SetBinLabel(i + 1, layer_names[i]);
	}
	leg_p1.Draw();
	c_p1.Write();

	TCanvas c_both("c_both", "Calibration Parameters Summary", 1400, 600);
	c_both.Divide(2, 1);
	c_both.cd(1);
	mg_p0.Draw("a lp");
	mg_p0.GetXaxis()->Set(5, -0.5, 4.5);
	mg_p0.GetXaxis()->SetNdivisions(5);
	for (int i = 0; i < kLayerCount; ++i) {
		mg_p0.GetXaxis()->SetBinLabel(i + 1, layer_names[i]);
	}
	leg_p0.Draw();
	c_both.cd(2);
	mg_p1.Draw("a lp");
	mg_p1.GetXaxis()->Set(5, -0.5, 4.5);
	mg_p1.GetXaxis()->SetNdivisions(5);
	for (int i = 0; i < kLayerCount; ++i) {
		mg_p1.GetXaxis()->SetBinLabel(i + 1, layer_names[i]);
	}
	leg_p1.Draw();
	c_both.Write();

	opf.Close();
	std::cout << "Output written to " << output_path << "\n";
	return 0;
}