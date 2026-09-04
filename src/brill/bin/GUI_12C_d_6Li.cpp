#include "include/config.h"
#include "include/C12_d_6Li/C12_d_6Li_event.h"
#include "include/Lise++/e_theta.h"
#include "include/utils.h"
#include "external/cxxopts.hpp"

#include <TApplication.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGClient.h>
#include <TGFileDialog.h>
#include <TGFrame.h>
#include <TGLayout.h>
#include <TGMenu.h>
#include <TGButton.h>
#include <TGLabel.h>
#include <TGNumberEntry.h>
#include <TGStatusBar.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TPad.h>
#include <TROOT.h>
#include <TRootEmbeddedCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

struct MainCanvases {
	TRootEmbeddedCanvas *embed = nullptr;
	TCanvas *canvas = nullptr;
	TH2D *h_e1_6Li_e2_6Li = nullptr;
	TH2D *h_e1_4He1_e2_4He1 = nullptr;
	TH2D *h_e1_4He2_e2_4He2 = nullptr;
	TH2D *h_e2_4He1_e3_4He1 = nullptr;
	TH2D *h_e2_4He2_e3_4He2 = nullptr;
	TH2D *h_e3_4He1_e4_4He1 = nullptr;
};

struct SecondCanvases {
	TCanvas *canvas = nullptr;
	TH2D *h_e3_4He2_e4_4He2 = nullptr;
	TH2D *h_e4sum_e5 = nullptr;
	TH2D *h_6Li_e_theta = nullptr;
	TH2D *h_6Li_e_theta_filtered = nullptr;
};

struct GUIContext {
	TGMainFrame *main_frame = nullptr;
	TGStatusBar *status_bar = nullptr;

	TGTextButton *btn_all = nullptr;
	TGTextButton *btn_14O = nullptr;
	TGTextButton *btn_13N = nullptr;
	TGTextButton *btn_12C = nullptr;
	TGTextButton *btn_draw = nullptr;

	TGNumberEntry *entry_run_min = nullptr;
	TGNumberEntry *entry_run_max = nullptr;

	MainCanvases mcs;
	SecondCanvases scs;

	std::string current_file;
	std::string config_path;
	std::string c12_d_6Li_dir;
	std::string assets_dir;
	std::vector<brill::C12D6LiEvent> all_events;
};

static GUIContext g_ctx;
static volatile int g_menu_action = 0;
static volatile int g_beam_filter = 0;
static volatile bool g_beam_changed = false;
static volatile bool g_redraw = false;

static const char *BeamLabel() {
	switch (g_beam_filter) {
		case 1: return "14O";
		case 2: return "13N";
		case 3: return "12C";
		default: return "All";
	}
}

static bool PassBeamFilter(const brill::C12D6LiEvent &ev) {
	switch (g_beam_filter) {
		case 1: return ev.is_14O;
		case 2: return ev.is_13N;
		case 3: return ev.is_12C;
		default: return true;
	}
}

static void RebuildMainHistograms() {
	auto &mc = g_ctx.mcs;
	if (mc.h_e1_6Li_e2_6Li) delete mc.h_e1_6Li_e2_6Li;
	if (mc.h_e1_4He1_e2_4He1) delete mc.h_e1_4He1_e2_4He1;
	if (mc.h_e1_4He2_e2_4He2) delete mc.h_e1_4He2_e2_4He2;
	if (mc.h_e2_4He1_e3_4He1) delete mc.h_e2_4He1_e3_4He1;
	if (mc.h_e2_4He2_e3_4He2) delete mc.h_e2_4He2_e3_4He2;
	if (mc.h_e3_4He1_e4_4He1) delete mc.h_e3_4He1_e4_4He1;

	TString label = TString::Format(" (%s)", BeamLabel());

	mc.h_e1_6Li_e2_6Li = new TH2D(
		"h_e1_6Li_e2_6Li",
		TString::Format("E1_6Li vs E2_6Li%s;E2_6Li (MeV);E1_6Li (MeV)", label.Data()),
		400, 0, 400, 400, 0, 25);
	mc.h_e1_6Li_e2_6Li->SetDirectory(0);
	mc.h_e1_4He1_e2_4He1 = new TH2D(
		"h_e1_4He1_e2_4He1",
		TString::Format("E1_4He1 vs E2_4He1%s;E2_4He1 (MeV);E1_4He1 (MeV)", label.Data()),
		400, 0, 400, 400, 0, 25);
	mc.h_e1_4He1_e2_4He1->SetDirectory(0);
	mc.h_e1_4He2_e2_4He2 = new TH2D(
		"h_e1_4He2_e2_4He2",
		TString::Format("E1_4He2 vs E2_4He2%s;E2_4He2 (MeV);E1_4He2 (MeV)", label.Data()),
		400, 0, 400, 400, 0, 25);
	mc.h_e1_4He2_e2_4He2->SetDirectory(0);
	mc.h_e2_4He1_e3_4He1 = new TH2D(
		"h_e2_4He1_e3_4He1",
		TString::Format("E2_4He1 vs E3_4He1%s;E3_4He1 (MeV);E2_4He1 (MeV)", label.Data()),
		400, 0, 350, 400, 0, 400);
	mc.h_e2_4He1_e3_4He1->SetDirectory(0);
	mc.h_e2_4He2_e3_4He2 = new TH2D(
		"h_e2_4He2_e3_4He2",
		TString::Format("E2_4He2 vs E3_4He2%s;E3_4He2 (MeV);E2_4He2 (MeV)", label.Data()),
		400, 0, 350, 400, 0, 400);
	mc.h_e2_4He2_e3_4He2->SetDirectory(0);
	mc.h_e3_4He1_e4_4He1 = new TH2D(
		"h_e3_4He1_e4_4He1",
		TString::Format("E3_4He1 vs E4_4He1%s;E4_4He1 (MeV);E3_4He1 (MeV)", label.Data()),
		400, 0, 250, 400, 0, 300);
	mc.h_e3_4He1_e4_4He1->SetDirectory(0);
}

static void RebuildSecondHistograms() {
	auto &sc = g_ctx.scs;
	if (sc.h_e3_4He2_e4_4He2) delete sc.h_e3_4He2_e4_4He2;
	if (sc.h_e4sum_e5) delete sc.h_e4sum_e5;
	if (sc.h_6Li_e_theta) delete sc.h_6Li_e_theta;

	TString label = TString::Format(" (%s)", BeamLabel());

	sc.h_e3_4He2_e4_4He2 = new TH2D(
		"h_e3_4He2_e4_4He2",
		TString::Format("E3_4He2 vs E4_4He2%s;E4_4He2 (MeV);E3_4He2 (MeV)", label.Data()),
		400, 0, 250, 400, 0, 300);
	sc.h_e3_4He2_e4_4He2->SetDirectory(0);
	sc.h_e4sum_e5 = new TH2D(
		"h_e4sum_e5",
		TString::Format("E4_4He1+E4_4He2 vs E5%s;E5 (MeV);E4_4He1+E4_4He2 (MeV)", label.Data()),
		400, 0, 500, 400, 0, 500);
	sc.h_e4sum_e5->SetDirectory(0);
	sc.h_6Li_e_theta = new TH2D(
		"h_6Li_e_theta",
		TString::Format("^{6}Li E-#theta%s;#theta_{6Li} (deg);E_{6Li} (MeV)", label.Data()),
		30, 0, 30, 200, 0, 400);
	sc.h_6Li_e_theta->SetDirectory(0);
}

static void FillAllHistograms() {
	RebuildMainHistograms();

	if (g_ctx.all_events.empty()) return;

	int total = (int)g_ctx.all_events.size();
	int passed = 0;
	for (int i = 0; i < total; ++i) {
		if (i % 100 == 0 || i == total - 1) {
			printf("\r  Filling all: %d/%d, passed=%d...", i + 1, total, passed);
			fflush(stdout);
			g_ctx.status_bar->SetText(
				TString::Format("Filling all: %d/%d, passed=%d", i + 1, total, passed));
			gSystem->ProcessEvents();
		}
		const auto &ev = g_ctx.all_events[i];
		if (!PassBeamFilter(ev)) continue;
		passed++;
		g_ctx.mcs.h_e1_6Li_e2_6Li->Fill(ev.e2_6Li, ev.e1_6Li);
		g_ctx.mcs.h_e1_4He1_e2_4He1->Fill(ev.e2_4He1, ev.e1_4He1);
		g_ctx.mcs.h_e1_4He2_e2_4He2->Fill(ev.e2_4He2, ev.e1_4He2);
		g_ctx.mcs.h_e2_4He1_e3_4He1->Fill(ev.e3_4He1, ev.e2_4He1);
		g_ctx.mcs.h_e2_4He2_e3_4He2->Fill(ev.e3_4He2, ev.e2_4He2);
		g_ctx.mcs.h_e3_4He1_e4_4He1->Fill(ev.e4_4He1, ev.e3_4He1);
	}
	printf("\r  All events: %d total, passed=%d        \n", total, passed);
}

static void FillSecondHistograms() {
	RebuildSecondHistograms();

	if (g_ctx.all_events.empty()) return;

	int total = (int)g_ctx.all_events.size();
	int passed = 0;

	for (int i = 0; i < total; ++i) {
		if (i % 100 == 0 || i == total - 1) {
			printf("\r  Filtering: %d/%d, passed=%d...", i + 1, total, passed);
			fflush(stdout);
			g_ctx.status_bar->SetText(
				TString::Format("Filtering: %d/%d, passed=%d", i + 1, total, passed));
			gSystem->ProcessEvents();
		}

		const auto &ev = g_ctx.all_events[i];

		if (!PassBeamFilter(ev)) continue;
		passed++;

		auto &sc = g_ctx.scs;
		sc.h_e3_4He2_e4_4He2->Fill(ev.e4_4He2, ev.e3_4He2);
		sc.h_e4sum_e5->Fill(ev.e5, ev.e4_4He1 + ev.e4_4He2);
		double T_6Li = ev.e1_6Li + ev.e2_6Li;
		sc.h_6Li_e_theta->Fill(ev.theta_6Li, T_6Li);
	}

	printf("\r  Filtering done: %d total, passed=%d        \n", total, passed);
	g_ctx.status_bar->SetText(
		TString::Format("Beam=%s: %d passed. File: %s",
			BeamLabel(), passed, g_ctx.current_file.c_str()));
}

static void RebuildFilteredSecondHistograms() {
	auto &sc = g_ctx.scs;
	if (sc.h_6Li_e_theta_filtered) delete sc.h_6Li_e_theta_filtered;

	TString label = TString::Format(" (%s)", BeamLabel());

	sc.h_6Li_e_theta_filtered = new TH2D(
		"h_6Li_e_theta_filtered",
		TString::Format("^{6}Li E-#theta (4He cut)%s;#theta_{6Li} (deg);E_{6Li} (MeV)", label.Data()),
		30, 0, 30, 200, 0, 400);
	sc.h_6Li_e_theta_filtered->SetDirectory(0);
}

static void FillFilteredSecondHistograms(TCutG *cut) {
	RebuildFilteredSecondHistograms();

	if (g_ctx.all_events.empty()) return;

	int total = (int)g_ctx.all_events.size();
	int passed = 0, n_beam = 0, n_ppac = 0, n_cut1 = 0, n_cut2 = 0;

	for (int i = 0; i < total; ++i) {
		if (i % 100 == 0 || i == total - 1) {
			printf("\r  Filtering (4He cut): %d/%d, passed=%d...", i + 1, total, passed);
			fflush(stdout);
			g_ctx.status_bar->SetText(
				TString::Format("Filtering (4He cut): %d/%d, passed=%d", i + 1, total, passed));
			gSystem->ProcessEvents();
		}

		const auto &ev = g_ctx.all_events[i];

		if (!PassBeamFilter(ev)) { n_beam++; continue; }
		if (!ev.ppac_valid) { n_ppac++; continue; }
		if (!cut->IsInside(ev.e4_4He1, ev.e3_4He1)) { n_cut1++; continue; }
		if (!cut->IsInside(ev.e4_4He2, ev.e3_4He2)) { n_cut2++; continue; }

		passed++;

		const char *beam = "?";
		if (ev.is_14O) beam = "14O";
		else if (ev.is_13N) beam = "13N";
		else if (ev.is_12C) beam = "12C";
		printf("run=%d  entry=%lld  beam=%s\n", ev.run_number, ev.entry, beam);

		double T_6Li = ev.e1_6Li + ev.e2_6Li;
		g_ctx.scs.h_6Li_e_theta_filtered->Fill(ev.theta_6Li, T_6Li);
	}

	printf("\r  Filtering (4He cut) done: %d total, beam=%d ppac=%d cut1=%d cut2=%d passed=%d\n",
		total, n_beam, n_ppac, n_cut1, n_cut2, passed);
}

static void DrawMainHistograms() {
	auto &mc = g_ctx.mcs;
	if (!mc.canvas) return;
	mc.canvas->Clear();
	mc.canvas->Divide(3, 2);
	mc.canvas->cd(1);
	mc.h_e1_6Li_e2_6Li->Draw("colz");
	mc.canvas->cd(2);
	mc.h_e1_4He1_e2_4He1->Draw("colz");
	mc.canvas->cd(3);
	mc.h_e1_4He2_e2_4He2->Draw("colz");
	mc.canvas->cd(4);
	mc.h_e2_4He1_e3_4He1->Draw("colz");
	mc.canvas->cd(5);
	mc.h_e2_4He2_e3_4He2->Draw("colz");
	mc.canvas->cd(6);
	mc.h_e3_4He1_e4_4He1->Draw("colz");
	mc.canvas->Modified();
	mc.canvas->Update();
}

static void DrawSecondHistograms() {
	auto &sc = g_ctx.scs;
	if (!sc.canvas) return;
	sc.canvas->Clear();
	sc.canvas->Divide(3, 2);

	TGraph *ref_6Li_0p = brill::LoadEThetaCurve(
		brill::JoinPath(g_ctx.assets_dir, "12C_d_6Li_0+_e_theta.txt"), "6Li");
	TGraph *ref_6Li_2p = brill::LoadEThetaCurve(
		brill::JoinPath(g_ctx.assets_dir, "12C_d_6Li_2+_e_theta.txt"), "6Li");
	TGraph *ref_6Li_4p = brill::LoadEThetaCurve(
		brill::JoinPath(g_ctx.assets_dir, "12C_d_6Li_4+_e_theta.txt"), "6Li");
	if (ref_6Li_2p) { ref_6Li_2p->SetLineColor(kGreen+2); ref_6Li_2p->SetLineStyle(kDashed); }
	if (ref_6Li_4p) { ref_6Li_4p->SetLineColor(kMagenta); ref_6Li_4p->SetLineStyle(kDotted); }

	sc.canvas->cd(1);
	sc.h_e3_4He2_e4_4He2->Draw("colz");
	sc.canvas->cd(2);
	sc.h_e4sum_e5->Draw("colz");
	sc.canvas->cd(3);
	sc.h_6Li_e_theta->Draw("colz");
	ref_6Li_0p->Draw("L same");
	ref_6Li_2p->Draw("L same");
	ref_6Li_4p->Draw("L same");

	sc.canvas->cd(6);
	if (sc.h_6Li_e_theta_filtered) {
		sc.h_6Li_e_theta_filtered->Draw("colz");
		ref_6Li_0p->Draw("L same");
		ref_6Li_2p->Draw("L same");
		ref_6Li_4p->Draw("L same");
	}

	sc.canvas->Modified();
	sc.canvas->Update();
}

static void OnBeamChanged() {
	if (g_ctx.current_file.empty()) return;
	FillAllHistograms();
	DrawMainHistograms();
	FillSecondHistograms();

	TString cut_path = "src/brill/Cut/cal_d3_d4_all_4He_cut.C";
	gROOT->ProcessLine(TString::Format(".x %s", cut_path.Data()));
	TCutG *cut = (TCutG*)gROOT->FindObject("cal_d3_d4_all_4He_cut");
	if (cut) {
		FillFilteredSecondHistograms(cut);
	}
	DrawSecondHistograms();
}

static void OnBeamChanged();
static void RebuildFromFile(const std::string &filepath);

void OnFileOpen() {
	TGFileInfo fi;
	const char *filetypes[] = {"ROOT files", "*.root", nullptr, nullptr};
	fi.fFileTypes = filetypes;
	fi.fIniDir = StrDup(g_ctx.c12_d_6Li_dir.c_str());
	new TGFileDialog(gClient->GetRoot(), g_ctx.main_frame, kFDOpen, &fi);

	if (!fi.fFilename) return;

	g_ctx.current_file = fi.fFilename;
	RebuildFromFile(g_ctx.current_file);
}

void RebuildFromFile(const std::string &filepath) {
	printf("Loading file: %s\n", filepath.c_str());

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

	brill::C12D6LiEvent ev;
	brill::SetupInputC12D6Li(input_tree, ev);

	g_ctx.all_events.clear();
	Long64_t n_entries = input_tree->GetEntries();

	int run_min = g_ctx.entry_run_min ? g_ctx.entry_run_min->GetIntNumber() : 0;
	int run_max = g_ctx.entry_run_max ? g_ctx.entry_run_max->GetIntNumber() : 0;
	bool apply_run_filter = (run_min > 0 || run_max > 0);

	printf("  Loading %lld events...\n", n_entries);
	if (apply_run_filter) {
		printf("  Run filter: %d - %d\n", run_min > 0 ? run_min : 0, run_max > 0 ? run_max : 99999);
	}
	int skipped = 0;
	for (Long64_t i = 0; i < n_entries; ++i) {
		input_tree->GetEntry(i);
		if (apply_run_filter) {
			if (run_min > 0 && ev.run_number < run_min) { skipped++; continue; }
			if (run_max > 0 && ev.run_number > run_max) { skipped++; continue; }
		}
		g_ctx.all_events.push_back(ev);
	}
	input_file->Close();
	printf("  Loaded %zu events", g_ctx.all_events.size());
	if (apply_run_filter) printf(" (skipped %d)", skipped);
	printf("\n");

	g_ctx.status_bar->SetText("Filling all histograms...");
	gSystem->ProcessEvents();

	FillAllHistograms();
	DrawMainHistograms();

	FillSecondHistograms();

	TString cut_path = "src/brill/Cut/cal_d3_d4_all_4He_cut.C";
	gROOT->ProcessLine(TString::Format(".x %s", cut_path.Data()));
	TCutG *cut = (TCutG*)gROOT->FindObject("cal_d3_d4_all_4He_cut");
	if (cut) {
		FillFilteredSecondHistograms(cut);
	}

	DrawSecondHistograms();

	printf("  Done.\n");
}

int main(int argc, char **argv) {
	cxxopts::Options options("GUI_12C_d_6Li", "GUI for 12C+d->6Li+2alpha physics analysis");
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
	g_ctx.c12_d_6Li_dir = brill::JoinPath(config.workspace, config.paths.c12_d_6Li);
	g_ctx.assets_dir = "assets";

	TApplication app("GUI_12C_d_6Li", &argc, argv);
	gStyle->SetPalette(kRainBow);

	gInterpreter->Declare(
		TString::Format("volatile int &g_menu_action = *((volatile int*)%lu);",
			(unsigned long)&g_menu_action).Data()
	);
	gInterpreter->Declare("void HandleMenuSlot(Int_t id) { g_menu_action = (int)id; }");
	gInterpreter->Declare(
		TString::Format("volatile int &g_beam_filter = *((volatile int*)%lu);",
			(unsigned long)&g_beam_filter).Data()
	);
	gInterpreter->Declare(
		TString::Format("volatile bool &g_beam_changed = *((volatile bool*)%lu);",
			(unsigned long)&g_beam_changed).Data()
	);
	gInterpreter->Declare(
		TString::Format("volatile bool &g_redraw = *((volatile bool*)%lu);",
			(unsigned long)&g_redraw).Data()
	);

	TGMainFrame *main_frame = new TGMainFrame(gClient->GetRoot(), 1200, 900);
	main_frame->SetWindowName("GUI_12C_d_6Li");
	g_ctx.main_frame = main_frame;

	TGPopupMenu *menu_file = new TGPopupMenu(gClient->GetRoot());
	menu_file->AddEntry("&Open...", 1);
	menu_file->AddSeparator();
	menu_file->AddEntry("&Quit", 2);
	menu_file->Connect("Activated(Int_t)", nullptr, nullptr, "HandleMenuSlot(Int_t)");

	TGMenuBar *menu_bar = new TGMenuBar(main_frame, 1, 1, kHorizontalFrame);
	menu_bar->AddPopup("&File", menu_file, new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 0, 0, 0));
	main_frame->AddFrame(menu_bar, new TGLayoutHints(kLHintsTop | kLHintsExpandX));

	TGHorizontalFrame *beam_frame = new TGHorizontalFrame(main_frame, 400, 30);
	TGLabel *beam_label = new TGLabel(beam_frame, "Beam: ");
	beam_frame->AddFrame(beam_label, new TGLayoutHints(kLHintsCenterY, 4, 2, 2, 2));

	g_ctx.btn_all = new TGTextButton(beam_frame, "All");
	g_ctx.btn_14O = new TGTextButton(beam_frame, "14O");
	g_ctx.btn_13N = new TGTextButton(beam_frame, "13N");
	g_ctx.btn_12C = new TGTextButton(beam_frame, "12C");

	g_ctx.btn_all->SetCommand("g_beam_filter = 0; g_beam_changed = true;");
	g_ctx.btn_14O->SetCommand("g_beam_filter = 1; g_beam_changed = true;");
	g_ctx.btn_13N->SetCommand("g_beam_filter = 2; g_beam_changed = true;");
	g_ctx.btn_12C->SetCommand("g_beam_filter = 3; g_beam_changed = true;");

	beam_frame->AddFrame(g_ctx.btn_all, new TGLayoutHints(kLHintsCenterY, 4, 2, 2, 2));
	beam_frame->AddFrame(g_ctx.btn_14O, new TGLayoutHints(kLHintsCenterY, 4, 2, 2, 2));
	beam_frame->AddFrame(g_ctx.btn_13N, new TGLayoutHints(kLHintsCenterY, 4, 2, 2, 2));
	beam_frame->AddFrame(g_ctx.btn_12C, new TGLayoutHints(kLHintsCenterY, 4, 2, 2, 2));

	TGLabel *run_label = new TGLabel(beam_frame, "  Run: ");
	beam_frame->AddFrame(run_label, new TGLayoutHints(kLHintsCenterY, 10, 2, 2, 2));

	TGNumberEntry *entry_run_min = new TGNumberEntry(beam_frame, 0, 5, -1,
		TGNumberFormat::kNESInteger,
		TGNumberFormat::kNEANonNegative,
		TGNumberFormat::kNELLimitMinMax, 0, 99999);
	beam_frame->AddFrame(entry_run_min, new TGLayoutHints(kLHintsCenterY, 2, 2, 2, 2));
	g_ctx.entry_run_min = entry_run_min;

	TGLabel *dash_label = new TGLabel(beam_frame, "-");
	beam_frame->AddFrame(dash_label, new TGLayoutHints(kLHintsCenterY, 2, 2, 2, 2));

	TGNumberEntry *entry_run_max = new TGNumberEntry(beam_frame, 0, 5, -1,
		TGNumberFormat::kNESInteger,
		TGNumberFormat::kNEANonNegative,
		TGNumberFormat::kNELLimitMinMax, 0, 99999);
	beam_frame->AddFrame(entry_run_max, new TGLayoutHints(kLHintsCenterY, 2, 2, 2, 2));
	g_ctx.entry_run_max = entry_run_max;

	TGTextButton *btn_draw = new TGTextButton(beam_frame, "Draw");
	btn_draw->SetCommand("g_redraw = true;");
	beam_frame->AddFrame(btn_draw, new TGLayoutHints(kLHintsCenterY, 10, 2, 2, 2));
	g_ctx.btn_draw = btn_draw;

	main_frame->AddFrame(beam_frame, new TGLayoutHints(kLHintsTop | kLHintsLeft, 4, 4, 2, 2));

	TRootEmbeddedCanvas *embed = new TRootEmbeddedCanvas("embed_main", main_frame, 1200, 700);
	main_frame->AddFrame(embed, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));
	g_ctx.mcs.embed = embed;
	g_ctx.mcs.canvas = embed->GetCanvas();

	g_ctx.scs.canvas = new TCanvas("canvas_second", "12C+d->6Li+2alpha Second", 1200, 800);

	TGStatusBar *status_bar = new TGStatusBar(main_frame, 1, 1);
	main_frame->AddFrame(status_bar, new TGLayoutHints(kLHintsBottom | kLHintsExpandX));
	g_ctx.status_bar = status_bar;
	status_bar->SetText("Ready. File > Open to load 12C+d->6Li+2alpha data.");

	RebuildMainHistograms();
	DrawMainHistograms();

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
		if (g_beam_changed) {
			g_beam_changed = false;
			OnBeamChanged();
		}
		if (g_redraw) {
			g_redraw = false;
			if (!g_ctx.current_file.empty()) {
				RebuildFromFile(g_ctx.current_file);
			}
		}
		g_menu_action = 0;
	}

	return 0;
}