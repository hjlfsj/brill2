#include "include/config.h"
#include "include/d_6Li/d_6Li_event.h"
#include "include/Lise++/e_theta.h"
#include "include/physics/kinematics.h"
#include "include/rebuild/rebuild_d_6Li.h"
#include "include/utils.h"
#include "external/cxxopts.hpp"

#include <TApplication.h>
#include <TCanvas.h>
#include <TCutG.h>
#include <TFile.h>
#include <TGClient.h>
#include <TGFileDialog.h>
#include <TGFrame.h>
#include <TGLayout.h>
#include <TGMenu.h>
#include <TGStatusBar.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TPad.h>
#include <TRootEmbeddedCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

struct BeamCanvases {
	TRootEmbeddedCanvas *embed = nullptr;
	TCanvas *canvas = nullptr;
	TH2D *h_e1_10C_e2_10C = nullptr;
	TH2D *h_e1_6Li_e2_6Li = nullptr;
	TH2D *h_e2_10C_e3_10C = nullptr;
	TH2D *h_e3_10C_e4_10C = nullptr;
};

struct AnalysisCanvas {
	TCanvas *canvas = nullptr;
	TH2D *h_E_6Li_E_10C = nullptr;
	TH2D *h_10C_e_theta = nullptr;
	TH2D *h_6Li_e_theta = nullptr;
	TH2D *h_theta_theta = nullptr;
};

struct GUIContext {
	TGMainFrame *main_frame = nullptr;
	TGStatusBar *status_bar = nullptr;

	BeamCanvases bcs;
	AnalysisCanvas ac;

	std::string current_file;
	std::string config_path;
	std::string d_Li6_dir;
	std::string assets_dir;
	std::vector<brill::D6LiEvent> all_events;
};

static GUIContext g_ctx;
static volatile int g_menu_action = 0;

static void RebuildHistograms() {
	auto &bc = g_ctx.bcs;
	if (bc.h_e1_10C_e2_10C) delete bc.h_e1_10C_e2_10C;
	if (bc.h_e1_6Li_e2_6Li) delete bc.h_e1_6Li_e2_6Li;
	if (bc.h_e2_10C_e3_10C) delete bc.h_e2_10C_e3_10C;
	if (bc.h_e3_10C_e4_10C) delete bc.h_e3_10C_e4_10C;

	bc.h_e1_10C_e2_10C = new TH2D(
		"h_e1_10C_e2_10C",
		"E1_10C vs E2_10C (all);E2_10C (MeV);E1_10C (MeV)",
		1000, 0, 400, 1000, 0, 25);
	bc.h_e1_10C_e2_10C->SetDirectory(0);
	bc.h_e1_6Li_e2_6Li = new TH2D(
		"h_e1_6Li_e2_6Li",
		"E1_6Li vs E2_6Li (all);E2_6Li (MeV);E1_6Li (MeV)",
		1000, 0, 400, 1000, 0, 25);
	bc.h_e1_6Li_e2_6Li->SetDirectory(0);
	bc.h_e2_10C_e3_10C = new TH2D(
		"h_e2_10C_e3_10C",
		"E2_10C vs E3_10C (all);E3_10C (MeV);E2_10C (MeV)",
		1000, 0, 350, 1000, 0, 400);
	bc.h_e2_10C_e3_10C->SetDirectory(0);
	bc.h_e3_10C_e4_10C = new TH2D(
		"h_e3_10C_e4_10C",
		"E3_10C vs E4_10C (all);E4_10C (MeV);E3_10C (MeV)",
		1000, 0, 250, 1000, 0, 300);
	bc.h_e3_10C_e4_10C->SetDirectory(0);
}

static void RebuildAnalysisHistograms() {
	auto &ac = g_ctx.ac;
	if (ac.h_E_6Li_E_10C) delete ac.h_E_6Li_E_10C;
	if (ac.h_10C_e_theta) delete ac.h_10C_e_theta;
	if (ac.h_6Li_e_theta) delete ac.h_6Li_e_theta;
	if (ac.h_theta_theta) delete ac.h_theta_theta;

	ac.h_E_6Li_E_10C = new TH2D(
		"h_E_6Li_E_10C",
		"E_{6Li} vs E_{10C};E_{10C} (MeV);E_{6Li} (MeV)",
		500, 0, 500, 500, 0, 200);
	ac.h_E_6Li_E_10C->SetDirectory(0);
	ac.h_10C_e_theta = new TH2D(
		"h_10C_e_theta",
		"^{10}C E-#theta;#theta_{10C} (deg);E_{10C} (MeV)",
		180, 0, 180, 500, 0, 500);
	ac.h_10C_e_theta->SetDirectory(0);
	ac.h_6Li_e_theta = new TH2D(
		"h_6Li_e_theta",
		"^{6}Li E-#theta;#theta_{6Li} (deg);E_{6Li} (MeV)",
		180, 0, 180, 500, 0, 500);
	ac.h_6Li_e_theta->SetDirectory(0);
	ac.h_theta_theta = new TH2D(
		"h_theta_theta",
		"#theta_{10C} vs #theta_{6Li};#theta_{10C} (deg);#theta_{6Li} (deg)",
		180, 0, 180, 180, 0, 180);
	ac.h_theta_theta->SetDirectory(0);
}

static void FillHistograms() {
	RebuildHistograms();

	if (g_ctx.all_events.empty()) return;

	int total_events = (int)g_ctx.all_events.size();
	int passed = 0;

	for (int i = 0; i < total_events; ++i) {
		if (i % 100 == 0 || i == total_events - 1) {
			printf("\r  Processing %d/%d events, passed=%d...",
				i + 1, total_events, passed);
			fflush(stdout);
			g_ctx.status_bar->SetText(
				TString::Format("Processing %d/%d events, passed=%d...",
					i + 1, total_events, passed));
			gSystem->ProcessEvents();
		}

		const auto &ev = g_ctx.all_events[i];

		passed++;
		g_ctx.bcs.h_e1_10C_e2_10C->Fill(ev.e2_10C, ev.e1_10C);
		g_ctx.bcs.h_e1_6Li_e2_6Li->Fill(ev.e2_6Li, ev.e1_6Li);
		g_ctx.bcs.h_e2_10C_e3_10C->Fill(ev.e3_10C, ev.e2_10C);
		g_ctx.bcs.h_e3_10C_e4_10C->Fill(ev.e4_10C, ev.e3_10C);
	}

	printf("\r  Done: %d events, passed=%d        \n", total_events, passed);
}

static void FillAnalysisHistograms(TCutG *cut) {
	RebuildAnalysisHistograms();

	if (g_ctx.all_events.empty()) return;

	int total_events = (int)g_ctx.all_events.size();
	int passed = 0;

	for (int i = 0; i < total_events; ++i) {
		if (i % 100 == 0 || i == total_events - 1) {
			printf("\r  Analysis %d/%d events, passed=%d...",
				i + 1, total_events, passed);
			fflush(stdout);
			g_ctx.status_bar->SetText(
				TString::Format("Analysis %d/%d events, passed=%d...",
					i + 1, total_events, passed));
			gSystem->ProcessEvents();
		}

		const auto &ev = g_ctx.all_events[i];

		if (!ev.ppac_valid) continue;
		if (!cut->IsInside(ev.e2_6Li, ev.e1_6Li)) continue;

		passed++;

		brill::D6LiKinematics kin = brill::ComputeKinematics(ev);
		g_ctx.ac.h_E_6Li_E_10C->Fill(kin.E_10C, kin.E_6Li);
		g_ctx.ac.h_10C_e_theta->Fill(ev.theta_10C, kin.E_10C);
		g_ctx.ac.h_6Li_e_theta->Fill(ev.theta_6Li, kin.E_6Li);
		g_ctx.ac.h_theta_theta->Fill(ev.theta_10C, ev.theta_6Li);
	}

	printf("\r  Analysis done: %d events, passed=%d (%d%%)        \n",
		total_events, passed,
		total_events > 0 ? (int)(passed * 100.0 / total_events) : 0);
	g_ctx.status_bar->SetText(
		TString::Format("Analysis: %d passed (%d%%). File: %s",
			passed,
			total_events > 0 ? (int)(passed * 100.0 / total_events) : 0,
			g_ctx.current_file.c_str()));
}

static void DrawHistograms() {
	auto &bc = g_ctx.bcs;
	if (!bc.canvas) return;
	bc.canvas->Clear();
	bc.canvas->Divide(2, 2);
	bc.canvas->cd(1);
	bc.h_e1_10C_e2_10C->Draw("colz");
	bc.canvas->cd(2);
	bc.h_e1_6Li_e2_6Li->Draw("colz");
	bc.canvas->cd(3);
	bc.h_e2_10C_e3_10C->Draw("colz");
	bc.canvas->cd(4);
	bc.h_e3_10C_e4_10C->Draw("colz");
	bc.canvas->Modified();
	bc.canvas->Update();
}

static void DrawAnalysisHistograms() {
	auto &ac = g_ctx.ac;
	if (!ac.canvas) return;
	ac.canvas->Clear();
	ac.canvas->Divide(2, 2);

	TGraph *ref_10C = brill::LoadEThetaCurve(
		brill::JoinPath(g_ctx.assets_dir, "14O_d_6Li_0+_e_theta.txt"), "10C");
	TGraph *ref_6Li = brill::LoadEThetaCurve(
		brill::JoinPath(g_ctx.assets_dir, "14O_d_6Li_0+_e_theta.txt"), "6Li");
	TGraph *ref_theta = brill::LoadThetaThetaCurve(
		brill::JoinPath(g_ctx.assets_dir, "14O_d_6Li_0+_theta_theta.txt"));

	ac.canvas->cd(1);
	ac.h_E_6Li_E_10C->Draw("colz");

	ac.canvas->cd(2);
	ac.h_10C_e_theta->Draw("colz");
	ref_10C->Draw("L same");

	ac.canvas->cd(3);
	ac.h_6Li_e_theta->Draw("colz");
	ref_6Li->Draw("L same");

	ac.canvas->cd(4);
	ac.h_theta_theta->Draw("colz");
	ref_theta->Draw("L same");

	ac.canvas->Modified();
	ac.canvas->Update();
}

void OnFileOpen() {
	TGFileInfo fi;
	const char *filetypes[] = {"ROOT files", "*.root", nullptr, nullptr};
	fi.fFileTypes = filetypes;
	fi.fIniDir = StrDup(g_ctx.d_Li6_dir.c_str());
	new TGFileDialog(gClient->GetRoot(), g_ctx.main_frame, kFDOpen, &fi);

	if (!fi.fFilename) return;

	std::string filepath = fi.fFilename;
	g_ctx.current_file = filepath;
	printf("Opening file: %s\n", filepath.c_str());

	TFile *input_file = new TFile(filepath.c_str(), "read");
	if (!input_file || input_file->IsZombie()) {
		std::cerr << "Error: Cannot open " << filepath << "\n";
		return;
	}

	TTree *input_tree = (TTree*)input_file->Get("tree");
	if (!input_tree) {
		std::cerr << "Error: No tree in " << filepath << "\n";
		input_file->Close();
		return;
	}

	brill::D6LiEvent ev;
	brill::SetupInputD6Li(input_tree, ev);

	g_ctx.all_events.clear();
	Long64_t n_entries = input_tree->GetEntries();
	printf("  Loading %lld events from d_Li6 file...\n", n_entries);
	for (Long64_t i = 0; i < n_entries; ++i) {
		input_tree->GetEntry(i);
		g_ctx.all_events.push_back(ev);
	}
	input_file->Close();

	printf("  Loaded %zu events\n", g_ctx.all_events.size());

	g_ctx.status_bar->SetText("Filling histograms...");
	gSystem->ProcessEvents();

	FillHistograms();
	DrawHistograms();

	std::string cut_path = "src/brill/Cut/cal_d1_d2_6Li_cut.C";
	printf("  Loading cut: %s\n", cut_path.c_str());
	TCutG *cut = brill::LoadCutGFromFile(cut_path);

	FillAnalysisHistograms(cut);
	DrawAnalysisHistograms();

	printf("  Done.\n");
}

int main(int argc, char **argv) {
	cxxopts::Options options("GUI_d_Li6", "GUI for d+6Li physics analysis");
	options.add_options()
		("c,config", "Config file path.", cxxopts::value<std::string>()->default_value("config.toml"))
		("h,help", "Print help information.");
	auto result = options.parse(argc, argv);

	if (result.count("help")) {
		std::cout << options.help() << "\n";
		return 0;
	}

	g_ctx.config_path = result["config"].as<std::string>();
	brill::AppConfig config;
	if (brill::LoadConfig(g_ctx.config_path, config)) {
		std::cerr << "Error: Load config failed.\n";
		return 1;
	}
	g_ctx.d_Li6_dir = brill::JoinPath(config.workspace, config.paths.d_Li6);
	g_ctx.assets_dir = "assets";

	TApplication app("GUI_d_Li6", &argc, argv);
	gStyle->SetPalette(kRainBow);

	gInterpreter->Declare(
		TString::Format("volatile int &g_menu_action = *((volatile int*)%lu);",
			(unsigned long)&g_menu_action).Data()
	);
	gInterpreter->Declare("void HandleMenuSlot(Int_t id) { g_menu_action = (int)id; }");

	TGMainFrame *main_frame = new TGMainFrame(gClient->GetRoot(), 1200, 900);
	main_frame->SetWindowName("GUI_d_Li6");
	g_ctx.main_frame = main_frame;

	TGPopupMenu *menu_file = new TGPopupMenu(gClient->GetRoot());
	menu_file->AddEntry("&Open...", 1);
	menu_file->AddSeparator();
	menu_file->AddEntry("&Quit", 2);
	menu_file->Connect("Activated(Int_t)", nullptr, nullptr, "HandleMenuSlot(Int_t)");

	TGMenuBar *menu_bar = new TGMenuBar(main_frame, 1, 1, kHorizontalFrame);
	menu_bar->AddPopup("&File", menu_file, new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 0, 0, 0));
	main_frame->AddFrame(menu_bar, new TGLayoutHints(kLHintsTop | kLHintsExpandX));

	TRootEmbeddedCanvas *embed = new TRootEmbeddedCanvas("embed_main", main_frame, 1200, 800);
	main_frame->AddFrame(embed, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));
	g_ctx.bcs.embed = embed;
	g_ctx.bcs.canvas = embed->GetCanvas();

	g_ctx.ac.canvas = new TCanvas("canvas_analysis", "d+6Li Analysis", 1200, 800);

	TGStatusBar *status_bar = new TGStatusBar(main_frame, 1, 1);
	main_frame->AddFrame(status_bar, new TGLayoutHints(kLHintsBottom | kLHintsExpandX));
	g_ctx.status_bar = status_bar;
	status_bar->SetText("Ready. File > Open to load d_Li6 data.");

	RebuildHistograms();
	DrawHistograms();

	main_frame->MapSubwindows();
	main_frame->Resize(main_frame->GetDefaultSize());
	main_frame->MapWindow();

	while (true) {
		gSystem->DispatchOneEvent(kFALSE);
		if (g_menu_action == 1) {
			OnFileOpen();
		} else if (g_menu_action == 2) {
			break;
		}
		g_menu_action = 0;
	}

	return 0;
}