#include "external/cxxopts.hpp"
#include "include/config.h"
#include "include/event/ingot/beam_event.h"
#include "include/event/beam/beam_sort.h"
#include "include/utils.h"

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TLine.h>
#include <TLatex.h>
#include <TCanvas.h>

#include <filesystem>
#include <iostream>
#include <string>

int main(int argc, char **argv) {
	cxxopts::Options options("sort_beam", "Sort beam by TOF peaks");
	options.add_options()
		("h,help", "Print help information.")
		("r,run", "Run number.", cxxopts::value<int>())
		("t,trigger", "Trigger type.", cxxopts::value<std::string>())
		("c,config", "Config file path.", cxxopts::value<std::string>()->default_value("config.toml"))
	;

	auto result = options.parse(argc, argv);

	if (result.count("help")) {
		std::cout << options.help() << std::endl;
		return 0;
	}

	if (!result.count("run") || !result.count("trigger")) {
		std::cerr << "Error: --run and --trigger are required.\n";
		std::cout << options.help() << std::endl;
		return 1;
	}

	int run = result["run"].as<int>();
	std::string trigger = result["trigger"].as<std::string>();
	std::string config_path = result["config"].as<std::string>();

	brill::AppConfig config;
	if (brill::LoadConfig(config_path, config)) {
		std::cerr << "Error: Load config failed.\n";
		return 1;
	}

	if (brill::IsJumpRun(config, run)) {
		std::cerr << "Error: Run " << run << " is a jump run.\n";
		return 1;
	}

	std::string ingot_dir = brill::JoinPath(config.workspace, config.paths.ingot);
	std::string output_dir = brill::JoinPath(config.workspace, config.paths.beam);
	std::string trigger_infix = brill::TriggerInfix(trigger);

	TString beam_path = TString::Format("%s/beam_%s%04d.root",
		ingot_dir.c_str(), trigger_infix.c_str(), run);

	if (!std::filesystem::exists(beam_path.Data())) {
		std::cerr << "Error: beam file not found: " << beam_path << "\n";
		return 1;
	}

	TFile *beam_file = new TFile(beam_path, "read");
	TTree *beam_tree = (TTree*)beam_file->Get("tree");
	if (!beam_tree) {
		std::cerr << "Error: tree not found in " << beam_path << "\n";
		beam_file->Close();
		return 1;
	}

	Long64_t n_entries = beam_tree->GetEntries();
	printf("Run %d: %lld events\n", run, n_entries);

	TCanvas *c1 = new TCanvas("c1_beam", "Beam TOF", 1200, 800);
	c1->cd();
	beam_tree->Draw("tof>>h_tof(5000,-500,500)", "valid", "goff");
	TH1D *h_drawn = (TH1D*)gDirectory->Get("h_tof");
	h_drawn->SetDirectory(0);
	if (!h_drawn || h_drawn->GetEntries() == 0) {
		std::cerr << "Error: No valid beam events found.\n";
		beam_file->Close();
		return 1;
	}

	brill::BeamSortResult sort_result;
	if (brill::SortBeamTOF(h_drawn, sort_result)) {
		std::cerr << "Error: Beam sorting failed.\n";
		beam_file->Close();
		return 1;
	}

	TString output_path = TString::Format("%s/beam_%s%04d.root",
		output_dir.c_str(), trigger_infix.c_str(), run);

	TFile *output_file = new TFile(output_path, "recreate");
	if (!output_file || output_file->IsZombie()) {
		std::cerr << "Error: Create output file failed: " << output_path << "\n";
		beam_file->Close();
		return 1;
	}

	h_drawn->Draw();

	TLine *l_14O_low = new TLine(sort_result.x_low_14O, 0, sort_result.x_low_14O, h_drawn->GetMaximum());
	l_14O_low->SetLineColor(kRed);
	l_14O_low->SetLineStyle(7);
	l_14O_low->Draw();

	TLine *l_14O_high = new TLine(sort_result.x_high_14O, 0, sort_result.x_high_14O, h_drawn->GetMaximum());
	l_14O_high->SetLineColor(kRed);
	l_14O_high->SetLineStyle(7);
	l_14O_high->Draw();

	TLine *l_13N_low = new TLine(sort_result.x_low_13N, 0, sort_result.x_low_13N, h_drawn->GetMaximum());
	l_13N_low->SetLineColor(kBlue);
	l_13N_low->SetLineStyle(7);
	l_13N_low->Draw();

	TLine *l_13N_high = new TLine(sort_result.x_high_13N, 0, sort_result.x_high_13N, h_drawn->GetMaximum());
	l_13N_high->SetLineColor(kBlue);
	l_13N_high->SetLineStyle(7);
	l_13N_high->Draw();

	TLine *l_12C_low = new TLine(sort_result.x_low_12C, 0, sort_result.x_low_12C, h_drawn->GetMaximum());
	l_12C_low->SetLineColor(kGreen+2);
	l_12C_low->SetLineStyle(7);
	l_12C_low->Draw();

	TLine *l_12C_high = new TLine(sort_result.x_high_12C, 0, sort_result.x_high_12C, h_drawn->GetMaximum());
	l_12C_high->SetLineColor(kGreen+2);
	l_12C_high->SetLineStyle(7);
	l_12C_high->Draw();

	TLine *l_part_12 = new TLine(sort_result.x_high_14O, 0, sort_result.x_high_14O, h_drawn->GetMaximum());
	l_part_12->SetLineColor(kBlack);
	l_part_12->SetLineWidth(2);
	l_part_12->SetLineStyle(2);
	l_part_12->Draw();

	TLine *l_part_23 = new TLine(sort_result.x_high_13N, 0, sort_result.x_high_13N, h_drawn->GetMaximum());
	l_part_23->SetLineColor(kBlack);
	l_part_23->SetLineWidth(2);
	l_part_23->SetLineStyle(2);
	l_part_23->Draw();

	double label_y = h_drawn->GetMaximum() * 0.85;
	TLatex *tex_14O = new TLatex(sort_result.mean_14O, label_y, "14O");
	tex_14O->SetTextAlign(22);
	tex_14O->SetTextColor(kRed);
	tex_14O->SetTextSize(0.04);
	tex_14O->Draw();

	TLatex *tex_13N = new TLatex(sort_result.mean_13N, label_y, "13N");
	tex_13N->SetTextAlign(22);
	tex_13N->SetTextColor(kBlue);
	tex_13N->SetTextSize(0.04);
	tex_13N->Draw();

	TLatex *tex_12C = new TLatex(sort_result.mean_12C, label_y, "12C");
	tex_12C->SetTextAlign(22);
	tex_12C->SetTextColor(kGreen+2);
	tex_12C->SetTextSize(0.04);
	tex_12C->Draw();

	c1->Write();

	TTree *output_tree = new TTree("tree", "sorted beam");
	bool v_14O = false;
	bool v_13N = false;
	bool v_12C = false;
	brill::SetupOutputSortBeamTree(output_tree, v_14O, v_13N, v_12C);

	brill::BeamEvent beam_event;
	brill::SetupInput(beam_tree, beam_event);

	printf("Run %d: %lld events", run, n_entries);
	fflush(stdout);

	Long64_t last_pct = 0;
	for (Long64_t i = 0; i < n_entries; ++i) {
		if (i * 100ll / n_entries > last_pct) {
			last_pct = i * 100ll / n_entries;
			printf("\b\b\b\b%3lld%%", last_pct);
			fflush(stdout);
		}
		beam_tree->GetEntry(i);

		brill::ClassifyBeam(beam_event.tof, beam_event.valid, sort_result, v_14O, v_13N, v_12C);

		output_tree->Fill();
	}
	printf("\b\b\b\b100%%\n");

	output_file->cd();
	output_tree->Write();
	h_drawn->Write();

	output_file->Close();
	beam_file->Close();

	printf("Output written to %s\n", output_path.Data());
	return 0;
}