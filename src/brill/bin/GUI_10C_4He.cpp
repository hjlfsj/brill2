#include "include/config.h"
#include "include/10C+4He/10C_4He_event.h"
#include "include/physics/kinematics.h"
#include "include/rebuild/nuclear_data.h"
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
#include <TGButton.h>
#include <TGLabel.h>
#include <TGNumberEntry.h>
#include <TGStatusBar.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TPad.h>
#include <TROOT.h>
#include <TRootEmbeddedCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

struct BeamCanvases {
	TRootEmbeddedCanvas *embed = nullptr;
	TCanvas *canvas = nullptr;
	TH2D *h_e1_10C_e2_10C = nullptr;
	TH2D *h_e2_10C_e3_10C = nullptr;
	TH2D *h_e1_4He_e2_4He = nullptr;
	TH2D *h_e2_4He_e3_4He = nullptr;
	TH2D *h_e3_4He_e4_4He = nullptr;
	TH2D *h_e4_4He_e5_4He = nullptr;
};

struct FilteredCanvases {
	TCanvas *canvas = nullptr;
	TH2D *h_e1_10C_e2_10C = nullptr;
	TH2D *h_e2_10C_e3_10C = nullptr;
	TH2D *h_e1_4He_e2_4He = nullptr;
	TH2D *h_e2_4He_e3_4He = nullptr;
	TH2D *h_e3_4He_e4_4He = nullptr;
	TH2D *h_e4_4He_e5_4He = nullptr;
};

struct AnalysisCanvas {
	TCanvas *canvas = nullptr;
	TH2D *h_10C_e_theta = nullptr;
	TH2D *h_4He_e_theta = nullptr;
	TH2D *h_theta_theta = nullptr;
	TH1D *h_excitation = nullptr;
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

	BeamCanvases bcs;
	FilteredCanvases fcs;
	AnalysisCanvas ac;

	std::string current_file;
	std::string config_path;
	std::string c10_he4_dir;
	std::vector<brill::C10He4Event> all_events;
};

static GUIContext g_ctx;
static volatile int g_menu_action = 0;
static volatile int g_beam_filter = 0;
static volatile bool g_beam_changed = false;
static volatile bool g_redraw = false;

static constexpr double kEbreak = 10.12;

static const char *BeamLabel() {
	switch (g_beam_filter) {
		case 1: return "14O";
		case 2: return "13N";
		case 3: return "12C";
		default: return "All";
	}
}

static bool PassBeamFilter(const brill::C10He4Event &ev) {
	switch (g_beam_filter) {
		case 1: return ev.is_14O;
		case 2: return ev.is_13N;
		case 3: return ev.is_12C;
		default: return true;
	}
}

static void RebuildBeamHistograms() {
	auto &bc = g_ctx.bcs;
	if (bc.h_e1_10C_e2_10C) delete bc.h_e1_10C_e2_10C;
	if (bc.h_e2_10C_e3_10C) delete bc.h_e2_10C_e3_10C;
	if (bc.h_e1_4He_e2_4He) delete bc.h_e1_4He_e2_4He;
	if (bc.h_e2_4He_e3_4He) delete bc.h_e2_4He_e3_4He;
	if (bc.h_e3_4He_e4_4He) delete bc.h_e3_4He_e4_4He;
	if (bc.h_e4_4He_e5_4He) delete bc.h_e4_4He_e5_4He;

	bc.h_e1_10C_e2_10C = new TH2D(
		"h_e1_10C_e2_10C",
		"E1_10C vs E2_10C (all);E2_10C (MeV);E1_10C (MeV)",
		400, 0, 400, 400, 0, 25);
	bc.h_e1_10C_e2_10C->SetDirectory(0);
	bc.h_e2_10C_e3_10C = new TH2D(
		"h_e2_10C_e3_10C",
		"E2_10C vs E3_10C (all);E3_10C (MeV);E2_10C (MeV)",
		400, 0, 350, 400, 0, 400);
	bc.h_e2_10C_e3_10C->SetDirectory(0);
	bc.h_e1_4He_e2_4He = new TH2D(
		"h_e1_4He_e2_4He",
		"E1_4He vs E2_4He (all);E2_4He (MeV);E1_4He (MeV)",
		400, 0, 400, 400, 0, 25);
	bc.h_e1_4He_e2_4He->SetDirectory(0);
	bc.h_e2_4He_e3_4He = new TH2D(
		"h_e2_4He_e3_4He",
		"E2_4He vs E3_4He (all);E3_4He (MeV);E2_4He (MeV)",
		400, 0, 350, 400, 0, 400);
	bc.h_e2_4He_e3_4He->SetDirectory(0);
	bc.h_e3_4He_e4_4He = new TH2D(
		"h_e3_4He_e4_4He",
		"E3_4He vs E4_4He (all);E4_4He (MeV);E3_4He (MeV)",
		400, 0, 250, 400, 0, 300);
	bc.h_e3_4He_e4_4He->SetDirectory(0);
	bc.h_e4_4He_e5_4He = new TH2D(
		"h_e4_4He_e5_4He",
		"E4_4He vs E5_4He (all);E5_4He (MeV);E4_4He (MeV)",
		400, 0, 250, 400, 0, 100);
	bc.h_e4_4He_e5_4He->SetDirectory(0);
}

static void RebuildFilteredHistograms() {
	auto &fc = g_ctx.fcs;
	if (fc.h_e1_10C_e2_10C) delete fc.h_e1_10C_e2_10C;
	if (fc.h_e2_10C_e3_10C) delete fc.h_e2_10C_e3_10C;
	if (fc.h_e1_4He_e2_4He) delete fc.h_e1_4He_e2_4He;
	if (fc.h_e2_4He_e3_4He) delete fc.h_e2_4He_e3_4He;
	if (fc.h_e3_4He_e4_4He) delete fc.h_e3_4He_e4_4He;
	if (fc.h_e4_4He_e5_4He) delete fc.h_e4_4He_e5_4He;

	TString label = TString::Format(" (%s)", BeamLabel());

	fc.h_e1_10C_e2_10C = new TH2D(
		"h_e1_10C_e2_10C_f",
		TString::Format("E1_10C vs E2_10C%s;E2_10C (MeV);E1_10C (MeV)", label.Data()),
		400, 0, 400, 400, 0, 25);
	fc.h_e1_10C_e2_10C->SetDirectory(0);
	fc.h_e2_10C_e3_10C = new TH2D(
		"h_e2_10C_e3_10C_f",
		TString::Format("E2_10C vs E3_10C%s;E3_10C (MeV);E2_10C (MeV)", label.Data()),
		400, 0, 350, 400, 0, 400);
	fc.h_e2_10C_e3_10C->SetDirectory(0);
	fc.h_e1_4He_e2_4He = new TH2D(
		"h_e1_4He_e2_4He_f",
		TString::Format("E1_4He vs E2_4He%s;E2_4He (MeV);E1_4He (MeV)", label.Data()),
		400, 0, 400, 400, 0, 25);
	fc.h_e1_4He_e2_4He->SetDirectory(0);
	fc.h_e2_4He_e3_4He = new TH2D(
		"h_e2_4He_e3_4He_f",
		TString::Format("E2_4He vs E3_4He%s;E3_4He (MeV);E2_4He (MeV)", label.Data()),
		400, 0, 350, 400, 0, 400);
	fc.h_e2_4He_e3_4He->SetDirectory(0);
	fc.h_e3_4He_e4_4He = new TH2D(
		"h_e3_4He_e4_4He_f",
		TString::Format("E3_4He vs E4_4He%s;E4_4He (MeV);E3_4He (MeV)", label.Data()),
		400, 0, 250, 400, 0, 300);
	fc.h_e3_4He_e4_4He->SetDirectory(0);
	fc.h_e4_4He_e5_4He = new TH2D(
		"h_e4_4He_e5_4He_f",
		TString::Format("E4_4He vs E5_4He%s;E5_4He (MeV);E4_4He (MeV)", label.Data()),
		400, 0, 250, 400, 0, 100);
	fc.h_e4_4He_e5_4He->SetDirectory(0);
}

static void RebuildAnalysisHistograms() {
	auto &ac = g_ctx.ac;
	if (ac.h_10C_e_theta) delete ac.h_10C_e_theta;
	if (ac.h_4He_e_theta) delete ac.h_4He_e_theta;
	if (ac.h_theta_theta) delete ac.h_theta_theta;
	if (ac.h_excitation) delete ac.h_excitation;

	TString label = TString::Format(" (%s)", BeamLabel());

	ac.h_10C_e_theta = new TH2D(
		"h_10C_e_theta",
		TString::Format("^{10}C E-#theta%s;#theta_{10C} (deg);E_{10C} (MeV)", label.Data()),
		20, 0, 20, 200, 220, 420);
	ac.h_10C_e_theta->SetDirectory(0);
	ac.h_4He_e_theta = new TH2D(
		"h_4He_e_theta",
		TString::Format("^{4}He E-#theta%s;#theta_{4He} (deg);E_{4He} (MeV)", label.Data()),
		25, 0, 25, 200, 0, 200);
	ac.h_4He_e_theta->SetDirectory(0);
	ac.h_theta_theta = new TH2D(
		"h_theta_theta",
		TString::Format("#theta_{10C} vs #theta_{4He}%s;#theta_{10C} (deg);#theta_{4He} (deg)", label.Data()),
		20, 0, 20, 30, 0, 30);
	ac.h_theta_theta->SetDirectory(0);
	ac.h_excitation = new TH1D(
		"h_excitation",
		TString::Format("^{14}O Excitation Energy%s;E_{x} (MeV);Counts", label.Data()),
		200, 5, 25);
	ac.h_excitation->SetDirectory(0);
}

static void FillAllHistograms() {
	RebuildBeamHistograms();

	if (g_ctx.all_events.empty()) return;

	int total = (int)g_ctx.all_events.size();
	for (int i = 0; i < total; ++i) {
		if (i % 100 == 0 || i == total - 1) {
			printf("\r  Filling all: %d/%d...", i + 1, total);
			fflush(stdout);
			g_ctx.status_bar->SetText(
				TString::Format("Filling all: %d/%d...", i + 1, total));
			gSystem->ProcessEvents();
		}
		const auto &ev = g_ctx.all_events[i];
		g_ctx.bcs.h_e1_10C_e2_10C->Fill(ev.e2_10C, ev.e1_10C);
		g_ctx.bcs.h_e2_10C_e3_10C->Fill(ev.e3_10C, ev.e2_10C);
		g_ctx.bcs.h_e1_4He_e2_4He->Fill(ev.e2_4He, ev.e1_4He);
		g_ctx.bcs.h_e2_4He_e3_4He->Fill(ev.e3_4He, ev.e2_4He);
		g_ctx.bcs.h_e3_4He_e4_4He->Fill(ev.e4_4He, ev.e3_4He);
		g_ctx.bcs.h_e4_4He_e5_4He->Fill(ev.e5_4He, ev.e4_4He);
	}
	printf("\r  All events: %d filled.        \n", total);
}

static void FillFilteredAndAnalysis(TCutG *cut) {
	RebuildFilteredHistograms();
	RebuildAnalysisHistograms();

	if (g_ctx.all_events.empty()) return;

	double m_10C = brill::GetMass(6, 10);
	double m_4He = brill::GetMass(2, 4);

	int total = (int)g_ctx.all_events.size();
	int passed = 0;
	int analyzed = 0;

	printf("  --- Selected events (run, entry) ---\n");
	for (int i = 0; i < total; ++i) {
		if (i % 100 == 0 || i == total - 1) {
			printf("\r  Filtering: %d/%d, passed=%d, analyzed=%d...",
				i + 1, total, passed, analyzed);
			fflush(stdout);
			g_ctx.status_bar->SetText(
				TString::Format("Filtering: %d/%d, passed=%d, analyzed=%d",
					i + 1, total, passed, analyzed));
			gSystem->ProcessEvents();
		}

		const auto &ev = g_ctx.all_events[i];

		if (!cut->IsInside(ev.e5_4He, ev.e4_4He)) continue;
		if (!PassBeamFilter(ev)) continue;

		passed++;
		auto &fc = g_ctx.fcs;
		fc.h_e1_10C_e2_10C->Fill(ev.e2_10C, ev.e1_10C);
		fc.h_e2_10C_e3_10C->Fill(ev.e3_10C, ev.e2_10C);
		fc.h_e1_4He_e2_4He->Fill(ev.e2_4He, ev.e1_4He);
		fc.h_e2_4He_e3_4He->Fill(ev.e3_4He, ev.e2_4He);
		fc.h_e3_4He_e4_4He->Fill(ev.e4_4He, ev.e3_4He);
		fc.h_e4_4He_e5_4He->Fill(ev.e5_4He, ev.e4_4He);

		if (!ev.ppac_valid) continue;

		printf("    run=%d  entry=%lld\n", ev.run_number, ev.entry);
		analyzed++;

		double T_10C = ev.e1_10C + ev.e2_10C + ev.e3_10C;
		double T_4He = ev.e1_4He + ev.e2_4He + ev.e3_4He + ev.e4_4He + ev.e5_4He;

		double E_10C = T_10C + m_10C;
		double E_4He = T_4He + m_4He;

		double p_10C = std::sqrt(E_10C * E_10C - m_10C * m_10C);
		double p_4He = std::sqrt(E_4He * E_4He - m_4He * m_4He);

		double cos_open = std::cos(ev.opening_angle * M_PI / 180.0);

		double M_inv_sq = m_10C * m_10C + m_4He * m_4He
			+ 2.0 * (E_10C * E_4He - p_10C * p_4He * cos_open);
		double M_inv = std::sqrt(M_inv_sq);

		double E_rel = M_inv - (m_10C + m_4He);
		double E_x = E_rel + kEbreak;

		auto &ac = g_ctx.ac;
		ac.h_10C_e_theta->Fill(ev.theta_10C, T_10C);
		ac.h_4He_e_theta->Fill(ev.theta_4He, T_4He);
		ac.h_theta_theta->Fill(ev.theta_10C, ev.theta_4He);
		ac.h_excitation->Fill(E_x);
	}

	printf("\r  Filtering done: %d total, passed=%d, analyzed=%d        \n",
		total, passed, analyzed);
	g_ctx.status_bar->SetText(
		TString::Format("Beam=%s: %d passed, %d analyzed. File: %s",
			BeamLabel(), passed, analyzed, g_ctx.current_file.c_str()));
}

static void DrawBeamHistograms() {
	auto &bc = g_ctx.bcs;
	if (!bc.canvas) return;
	bc.canvas->Clear();
	bc.canvas->Divide(3, 2);
	bc.canvas->cd(1);
	bc.h_e1_10C_e2_10C->Draw("colz");
	bc.canvas->cd(2);
	bc.h_e2_10C_e3_10C->Draw("colz");
	bc.canvas->cd(3);
	bc.h_e1_4He_e2_4He->Draw("colz");
	bc.canvas->cd(4);
	bc.h_e2_4He_e3_4He->Draw("colz");
	bc.canvas->cd(5);
	bc.h_e3_4He_e4_4He->Draw("colz");
	bc.canvas->cd(6);
	bc.h_e4_4He_e5_4He->Draw("colz");
	bc.canvas->Modified();
	bc.canvas->Update();
}

static void DrawFilteredHistograms() {
	auto &fc = g_ctx.fcs;
	if (!fc.canvas) return;
	fc.canvas->Clear();
	fc.canvas->Divide(3, 2);
	fc.canvas->cd(1);
	fc.h_e1_10C_e2_10C->Draw("colz");
	fc.canvas->cd(2);
	fc.h_e2_10C_e3_10C->Draw("colz");
	fc.canvas->cd(3);
	fc.h_e1_4He_e2_4He->Draw("colz");
	fc.canvas->cd(4);
	fc.h_e2_4He_e3_4He->Draw("colz");
	fc.canvas->cd(5);
	fc.h_e3_4He_e4_4He->Draw("colz");
	fc.canvas->cd(6);
	fc.h_e4_4He_e5_4He->Draw("colz");
	fc.canvas->Modified();
	fc.canvas->Update();
}

static void DrawAnalysisHistograms() {
	auto &ac = g_ctx.ac;
	if (!ac.canvas) return;
	ac.canvas->Clear();
	ac.canvas->Divide(2, 2);
	ac.canvas->cd(1);
	ac.h_10C_e_theta->Draw("colz");
	ac.canvas->cd(2);
	ac.h_4He_e_theta->Draw("colz");
	ac.canvas->cd(3);
	ac.h_theta_theta->Draw("colz");
	ac.canvas->cd(4);
	ac.h_excitation->Draw();
	ac.canvas->Modified();
	ac.canvas->Update();
}

static void OnBeamChanged() {
	if (g_ctx.current_file.empty()) return;

	TString cut_path = "src/brill/Cut/cal_d4_s1_stop_4He_cut.C";
	gROOT->ProcessLine(TString::Format(".x %s", cut_path.Data()));
	TCutG *cut = (TCutG*)gROOT->FindObject("cal_d4_s1_stop_4He_cut");
	if (!cut) {
		std::cerr << "Error: Load cut failed.\n";
		return;
	}

	FillFilteredAndAnalysis(cut);
	DrawFilteredHistograms();
	DrawAnalysisHistograms();
}

static void OnBeamChanged();
static void RebuildFromFile(const std::string &filepath);

void OnFileOpen() {
	TGFileInfo fi;
	const char *filetypes[] = {"ROOT files", "*.root", nullptr, nullptr};
	fi.fFileTypes = filetypes;
	fi.fIniDir = StrDup(g_ctx.c10_he4_dir.c_str());
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

	brill::C10He4Event ev;
	brill::SetupInputC10He4(input_tree, ev);

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
	DrawBeamHistograms();

	TString cut_path = "src/brill/Cut/cal_d4_s1_stop_4He_cut.C";
	printf("  Loading cut: %s\n", cut_path.Data());
	gROOT->ProcessLine(TString::Format(".x %s", cut_path.Data()));
	TCutG *cut = (TCutG*)gROOT->FindObject("cal_d4_s1_stop_4He_cut");
	if (!cut) {
		std::cerr << "Error: Load cut failed.\n";
		return;
	}

	FillFilteredAndAnalysis(cut);
	DrawFilteredHistograms();
	DrawAnalysisHistograms();

	printf("  Done.\n");
}

int main(int argc, char **argv) {
	cxxopts::Options options("GUI_10C_4He", "GUI for 10C+4He physics analysis");
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
	g_ctx.c10_he4_dir = brill::JoinPath(config.workspace, config.paths.c10_he4);

	TApplication app("GUI_10C_4He", &argc, argv);
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
	main_frame->SetWindowName("GUI_10C_4He");
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
	g_ctx.bcs.embed = embed;
	g_ctx.bcs.canvas = embed->GetCanvas();

	g_ctx.fcs.canvas = new TCanvas("canvas_filtered", "10C+4He Filtered", 1200, 800);
	g_ctx.ac.canvas = new TCanvas("canvas_analysis", "10C+4He Analysis", 1200, 800);

	TGStatusBar *status_bar = new TGStatusBar(main_frame, 1, 1);
	main_frame->AddFrame(status_bar, new TGLayoutHints(kLHintsBottom | kLHintsExpandX));
	g_ctx.status_bar = status_bar;
	status_bar->SetText("Ready. File > Open to load 10C+4He data.");

	RebuildBeamHistograms();
	DrawBeamHistograms();

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