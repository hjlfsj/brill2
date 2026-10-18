#include "include/config.h"
#include "include/three_alpha/three_alpha_event.h"
#include "include/rebuild/nuclear_data.h"
#include "include/utils.h"
#include "external/cxxopts.hpp"

#include <TApplication.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGButton.h>
#include <TGClient.h>
#include <TGFileDialog.h>
#include <TGFrame.h>
#include <TGLabel.h>
#include <TGLayout.h>
#include <TGMenu.h>
#include <TGNumberEntry.h>
#include <TGStatusBar.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLine.h>
#include <TLorentzVector.h>
#include <TPad.h>
#include <TRootEmbeddedCanvas.h>
#include <TSystem.h>
#include <TString.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

struct PIDCanvas {
	TCanvas *canvas = nullptr;
	TH2D *h_e1_e2 = nullptr;
	TH2D *h_e2_e3 = nullptr;
	TH2D *h_e3_e4 = nullptr;
	TH2D *h_e_theta = nullptr;
};

struct AnalysisPlots {
	TH1D *h_ex_8Be = nullptr;
	TH1D *h_ex_12C = nullptr;
	TH2D *h_ex_8Be_ex_12C = nullptr;
	TH1D *h_ex_12C_gated = nullptr;
	TLine *line_x1 = nullptr;
	TLine *line_x2 = nullptr;
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
	TGNumberEntry *entry_n1 = nullptr;
	TGNumberEntry *entry_n2 = nullptr;
	TGNumberEntry *entry_x1 = nullptr;
	TGNumberEntry *entry_x2 = nullptr;

	PIDCanvas pc_He1;
	PIDCanvas pc_He2;
	PIDCanvas pc_He3;

	TCanvas *analysis_canvas = nullptr;
	AnalysisPlots ap;

	std::string current_file;
	std::string config_path;
	std::string three_alpha_dir;
	std::string assets_dir;
	brill::AppConfig config;
	std::vector<brill::ThreeAlphaEvent> all_events;
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

static bool PassBeamFilter(const brill::ThreeAlphaEvent &ev) {
	switch (g_beam_filter) {
		case 1: return ev.is_14O;
		case 2: return ev.is_13N;
		case 3: return ev.is_12C;
		default: return true;
	}
}

static void RebuildPIDHistograms(PIDCanvas &pc, const char *name_prefix, const char *title_prefix) {
	delete pc.h_e1_e2; pc.h_e1_e2 = nullptr;
	delete pc.h_e2_e3; pc.h_e2_e3 = nullptr;
	delete pc.h_e3_e4; pc.h_e3_e4 = nullptr;
	delete pc.h_e_theta; pc.h_e_theta = nullptr;

	TString label = TString::Format(" (%s)", BeamLabel());

	pc.h_e1_e2 = new TH2D(
		TString::Format("h_%s_e1_e2", name_prefix),
		TString::Format("%s E1 vs E2%s;E2 (MeV);E1 (MeV)", title_prefix, label.Data()),
		1000, 0, 400, 1000, 0, 25);
	pc.h_e1_e2->SetDirectory(nullptr);

	pc.h_e2_e3 = new TH2D(
		TString::Format("h_%s_e2_e3", name_prefix),
		TString::Format("%s E2 vs E3%s;E3 (MeV);E2 (MeV)", title_prefix, label.Data()),
		1000, 0, 350, 1000, 0, 400);
	pc.h_e2_e3->SetDirectory(nullptr);

	pc.h_e3_e4 = new TH2D(
		TString::Format("h_%s_e3_e4", name_prefix),
		TString::Format("%s E3 vs E4%s;E4 (MeV);E3 (MeV)", title_prefix, label.Data()),
		1000, 0, 250, 1000, 0, 300);
	pc.h_e3_e4->SetDirectory(nullptr);

	pc.h_e_theta = new TH2D(
		TString::Format("h_%s_e_theta", name_prefix),
		TString::Format("%s E_{true} vs #theta%s;#theta (deg);E_{true} (MeV)", title_prefix, label.Data()),
		50, 0, 50, 100, 0, 200);
	pc.h_e_theta->SetDirectory(nullptr);
}

static void RebuildAnalysisHistograms() {
	auto &ap = g_ctx.ap;
	delete ap.h_ex_8Be; ap.h_ex_8Be = nullptr;
	delete ap.h_ex_12C; ap.h_ex_12C = nullptr;
	delete ap.h_ex_8Be_ex_12C; ap.h_ex_8Be_ex_12C = nullptr;
	delete ap.h_ex_12C_gated; ap.h_ex_12C_gated = nullptr;
	delete ap.line_x1; ap.line_x1 = nullptr;
	delete ap.line_x2; ap.line_x2 = nullptr;

	int n1 = g_ctx.entry_n1 ? g_ctx.entry_n1->GetIntNumber() : 300;
	int n2 = g_ctx.entry_n2 ? g_ctx.entry_n2->GetIntNumber() : 300;

	ap.h_ex_8Be = new TH1D("h_ex_8Be", "^{8}Be Excitation;E_{x} (MeV);Counts",
		n1, -1, 29);
	ap.h_ex_8Be->SetDirectory(nullptr);

	ap.h_ex_12C = new TH1D("h_ex_12C", "^{12}C Excitation;E_{x} (MeV);Counts",
		n2, -1, 29);
	ap.h_ex_12C->SetDirectory(nullptr);

	ap.h_ex_8Be_ex_12C = new TH2D("h_ex_8Be_ex_12C",
		"^{12}C vs ^{8}Be Excitation;E_{x}(^{12}C) (MeV);E_{x}(^{8}Be) (MeV)",
		n2, -1, 29, n1, -1, 29);
	ap.h_ex_8Be_ex_12C->SetDirectory(nullptr);

	ap.h_ex_12C_gated = new TH1D("h_ex_12C_gated",
		"^{12}C Excitation (^{8}Be gated);E_{x} (MeV);Counts",
		n2, -1, 29);
	ap.h_ex_12C_gated->SetDirectory(nullptr);

	double x1v = g_ctx.entry_x1 ? g_ctx.entry_x1->GetNumber() : -0.5;
	double x2v = g_ctx.entry_x2 ? g_ctx.entry_x2->GetNumber() : 0.5;
	ap.line_x1 = new TLine(x1v, 0, x1v, 1);
	ap.line_x1->SetLineColor(kRed);
	ap.line_x1->SetLineWidth(2);
	ap.line_x2 = new TLine(x2v, 0, x2v, 1);
	ap.line_x2->SetLineColor(kRed);
	ap.line_x2->SetLineWidth(2);
}

static void FillHistograms() {
	RebuildPIDHistograms(g_ctx.pc_He1, "He1", "He1");
	RebuildPIDHistograms(g_ctx.pc_He2, "He2", "He2");
	RebuildPIDHistograms(g_ctx.pc_He3, "He3", "He3");
	RebuildAnalysisHistograms();

	if (g_ctx.all_events.empty()) return;

	int total_events = (int)g_ctx.all_events.size();
	int passed = 0;
	double x1 = g_ctx.entry_x1 ? g_ctx.entry_x1->GetNumber() : -0.5;
	double x2 = g_ctx.entry_x2 ? g_ctx.entry_x2->GetNumber() : 0.5;
	if (x1 > x2) std::swap(x1, x2);

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
		if (!PassBeamFilter(ev)) continue;
		if (!ev.ppac_valid) continue;

		passed++;

		g_ctx.pc_He1.h_e1_e2->Fill(ev.e2_He1, ev.e1_He1);
		g_ctx.pc_He1.h_e2_e3->Fill(ev.e3_He1, ev.e2_He1);
		g_ctx.pc_He1.h_e3_e4->Fill(ev.e4_He1, ev.e3_He1);
		g_ctx.pc_He1.h_e_theta->Fill(ev.theta_He1, ev.et_He1);

		g_ctx.pc_He2.h_e1_e2->Fill(ev.e2_He2, ev.e1_He2);
		g_ctx.pc_He2.h_e2_e3->Fill(ev.e3_He2, ev.e2_He2);
		g_ctx.pc_He2.h_e3_e4->Fill(ev.e4_He2, ev.e3_He2);
		g_ctx.pc_He2.h_e_theta->Fill(ev.theta_He2, ev.et_He2);

		g_ctx.pc_He3.h_e1_e2->Fill(ev.e2_He3, ev.e1_He3);
		g_ctx.pc_He3.h_e2_e3->Fill(ev.e3_He3, ev.e2_He3);
		g_ctx.pc_He3.h_e3_e4->Fill(ev.e4_He3, ev.e3_He3);
		g_ctx.pc_He3.h_e_theta->Fill(ev.theta_He3, ev.et_He3);

		const double m_4He = brill::GetMass(2, 4);
		const double m_8Be_gs = brill::GetMass(4, 8);
		const double m_12C_gs = brill::GetMass(6, 12);

		double dx[3], dy[3], dz[3], norm[3];
		dx[0] = ev.t0d2_He1_x - ev.target_x;
		dy[0] = ev.t0d2_He1_y - ev.target_y;
		dz[0] = ev.t0d2_He1_z;
		dx[1] = ev.t0d2_He2_x - ev.target_x;
		dy[1] = ev.t0d2_He2_y - ev.target_y;
		dz[1] = ev.t0d2_He2_z;
		dx[2] = ev.t0d2_He3_x - ev.target_x;
		dy[2] = ev.t0d2_He3_y - ev.target_y;
		dz[2] = ev.t0d2_He3_z;

		double et[3] = {ev.et_He1, ev.et_He2, ev.et_He3};
		TLorentzVector P_alpha[3];

		for (int p = 0; p < 3; ++p) {
			norm[p] = std::sqrt(dx[p]*dx[p] + dy[p]*dy[p] + dz[p]*dz[p]);
			if (norm[p] < 1e-6) norm[p] = 1e-6;
			double E_tot = et[p] + m_4He;
			double momentum = std::sqrt(E_tot * E_tot - m_4He * m_4He);
			P_alpha[p] = TLorentzVector(
				momentum * dx[p] / norm[p],
				momentum * dy[p] / norm[p],
				momentum * dz[p] / norm[p],
				E_tot);
		}

		double opening[3] = {
			ev.opening_He1_He2,
			ev.opening_He1_He3,
			ev.opening_He2_He3
		};

		int min_idx = 0;
		if (opening[1] < opening[min_idx]) min_idx = 1;
		if (opening[2] < opening[min_idx]) min_idx = 2;

		int a = 0, b = 1;
		if (min_idx == 0) { a = 0; b = 1; }
		else if (min_idx == 1) { a = 0; b = 2; }
		else { a = 1; b = 2; }

		TLorentzVector P_8Be = P_alpha[a] + P_alpha[b];
		double M_inv_8Be = P_8Be.M();
		double Ex_8Be = M_inv_8Be - m_8Be_gs;

		TLorentzVector P_12C = P_alpha[0] + P_alpha[1] + P_alpha[2];
		double M_inv_12C = P_12C.M();
		double Ex_12C = M_inv_12C - m_12C_gs;

		g_ctx.ap.h_ex_8Be->Fill(Ex_8Be);
		g_ctx.ap.h_ex_12C->Fill(Ex_12C);
		g_ctx.ap.h_ex_8Be_ex_12C->Fill(Ex_12C, Ex_8Be);

		if (Ex_8Be >= x1 && Ex_8Be <= x2) {
			g_ctx.ap.h_ex_12C_gated->Fill(Ex_12C);
		}
	}

	printf("\r  Done: %d events, passed=%d        \n", total_events, passed);
	g_ctx.status_bar->SetText(
		TString::Format("Processed %d events, %d passed.", total_events, passed));
}

static void DrawPIDCanvas(PIDCanvas &pc) {
	if (!pc.canvas) return;
	pc.canvas->Clear();
	pc.canvas->Divide(2, 2);
	pc.canvas->cd(1);
	if (pc.h_e1_e2) pc.h_e1_e2->Draw("colz");
	pc.canvas->cd(2);
	if (pc.h_e2_e3) pc.h_e2_e3->Draw("colz");
	pc.canvas->cd(3);
	if (pc.h_e3_e4) pc.h_e3_e4->Draw("colz");
	pc.canvas->cd(4);
	if (pc.h_e_theta) pc.h_e_theta->Draw("colz");
	pc.canvas->Modified();
	pc.canvas->Update();
}

static void DrawAnalysisCanvas() {
	auto &ap = g_ctx.ap;
	TCanvas *c = g_ctx.analysis_canvas;
	if (!c) return;

	c->Clear();
	c->Divide(2, 2);

	c->cd(1);
	if (ap.h_ex_8Be) {
		ap.h_ex_8Be->Draw();
		double ymax = ap.h_ex_8Be->GetMaximum() * 1.15;
		delete ap.line_x1;
		delete ap.line_x2;
		double x1 = g_ctx.entry_x1 ? g_ctx.entry_x1->GetNumber() : -0.5;
		double x2 = g_ctx.entry_x2 ? g_ctx.entry_x2->GetNumber() : 0.5;
		ap.line_x1 = new TLine(x1, 0, x1, ymax);
		ap.line_x1->SetLineColor(kRed);
		ap.line_x1->SetLineWidth(2);
		ap.line_x1->Draw();
		ap.line_x2 = new TLine(x2, 0, x2, ymax);
		ap.line_x2->SetLineColor(kRed);
		ap.line_x2->SetLineWidth(2);
		ap.line_x2->Draw();
	}

	c->cd(2);
	if (ap.h_ex_8Be_ex_12C) ap.h_ex_8Be_ex_12C->Draw("colz");

	c->cd(3);
	if (ap.h_ex_12C) ap.h_ex_12C->Draw();

	c->cd(4);
	if (ap.h_ex_12C_gated) ap.h_ex_12C_gated->Draw();

	c->Modified();
	c->Update();
}

static void OnDraw() {
	g_ctx.status_bar->SetText("Drawing...");
	gSystem->ProcessEvents();

	FillHistograms();

	DrawPIDCanvas(g_ctx.pc_He1);
	DrawPIDCanvas(g_ctx.pc_He2);
	DrawPIDCanvas(g_ctx.pc_He3);
	DrawAnalysisCanvas();

	g_ctx.status_bar->SetText(
		TString::Format("Draw done. File: %s", g_ctx.current_file.c_str()));
}

static int LoadAllEvents(const std::string &file_path) {
	TFile *f = TFile::Open(file_path.c_str(), "read");
	if (!f || f->IsZombie()) {
		std::cerr << "Error: Cannot open file " << file_path << "\n";
		return -1;
	}

	TTree *tree = (TTree*)f->Get("tree");
	if (!tree) {
		std::cerr << "Error: No tree found in " << file_path << "\n";
		f->Close();
		return -1;
	}

	Long64_t n_entries = tree->GetEntries();
	g_ctx.all_events.clear();
	g_ctx.all_events.resize(n_entries);

	brill::ThreeAlphaEvent ev;
	brill::SetupInputThreeAlpha(tree, ev);

	for (Long64_t i = 0; i < n_entries; ++i) {
		tree->GetEntry(i);
		g_ctx.all_events[i] = ev;
	}

	f->Close();
	g_ctx.current_file = file_path;
	printf("Loaded %zd events from %s\n", g_ctx.all_events.size(), file_path.c_str());
	fflush(stdout);
	return 0;
}

static void OnFileOpen() {
	TGFileInfo fi;
	const char *filetypes[] = {"ROOT files", "*.root", nullptr, nullptr};
	fi.fFileTypes = filetypes;
	fi.fIniDir = StrDup(g_ctx.three_alpha_dir.c_str());
	new TGFileDialog(gClient->GetRoot(), g_ctx.main_frame, kFDOpen, &fi);

	if (!fi.fFilename) return;

	std::string filepath = fi.fFilename;
	printf("Opening file: %s\n", filepath.c_str());
	fflush(stdout);

	if (LoadAllEvents(filepath) != 0) return;
	OnDraw();
}

static void ProcessMenu(int action) {
	switch (action) {
		case 3: { // Draw
			int run_min = g_ctx.entry_run_min->GetIntNumber();
			int run_max = g_ctx.entry_run_max->GetIntNumber();
			if (run_min > run_max) std::swap(run_min, run_max);

			TString file_path = TString::Format(
				"%s/extract_three_alpha_main_%04d_%04d.root",
				g_ctx.three_alpha_dir.c_str(), run_min, run_max);

			if (!std::filesystem::exists(file_path.Data())) {
				file_path = TString::Format(
					"%s/extract_three_alpha_t1_%04d_%04d.root",
					g_ctx.three_alpha_dir.c_str(), run_min, run_max);
			}

			if (!std::filesystem::exists(file_path.Data())) {
				std::cerr << "Error: Output file not found for runs "
					<< run_min << "-" << run_max << "\n";
				std::cerr << "  Tried: " << file_path << "\n";
				return;
			}

			if (LoadAllEvents(file_path.Data()) != 0) return;
			OnDraw();
			break;
		}
		case 4: // 14O
			g_beam_filter = (g_beam_filter == 1) ? 0 : 1;
			g_beam_changed = true;
			break;
		case 5: // 13N
			g_beam_filter = (g_beam_filter == 2) ? 0 : 2;
			g_beam_changed = true;
			break;
		case 6: // 12C
			g_beam_filter = (g_beam_filter == 3) ? 0 : 3;
			g_beam_changed = true;
			break;
		case 7: // All
			g_beam_filter = 0;
			g_beam_changed = true;
			break;
		default:
			break;
	}

	if (g_beam_changed) {
		g_beam_changed = false;
		const char *label = BeamLabel();
		g_ctx.btn_all->SetState(g_beam_filter == 0 ? kButtonDown : kButtonUp);
		g_ctx.btn_14O->SetState(g_beam_filter == 1 ? kButtonDown : kButtonUp);
		g_ctx.btn_13N->SetState(g_beam_filter == 2 ? kButtonDown : kButtonUp);
		g_ctx.btn_12C->SetState(g_beam_filter == 3 ? kButtonDown : kButtonUp);
		g_ctx.status_bar->SetText(
			TString::Format("Beam: %s. Ready.", label));
		g_redraw = true;
	}
}

int main(int argc, char **argv) {
	TApplication app("GUI_three_alpha", &argc, argv);

	cxxopts::Options options("GUI_three_alpha", "GUI for three-alpha analysis");
	options.add_options()
		("c,config", "Config file path.", cxxopts::value<std::string>()->default_value("config.toml"))
		("h,help", "Print help information.")
	;

	auto result = options.parse(argc, argv);
	if (result.count("help")) {
		std::cout << options.help() << std::endl;
		return 0;
	}

	g_ctx.config_path = result["config"].as<std::string>();
	if (brill::LoadConfig(g_ctx.config_path, g_ctx.config)) {
		std::cerr << "Error: Load config failed.\n";
		return 1;
	}

	g_ctx.three_alpha_dir = brill::JoinPath(g_ctx.config.workspace, g_ctx.config.paths.three_alpha);
	g_ctx.assets_dir = brill::JoinPath(g_ctx.config.workspace, g_ctx.config.assets);

	const int main_w = 1700;
	const int main_h = 900;

	g_ctx.main_frame = new TGMainFrame(gClient->GetRoot(), 50, 50, kMainFrame | kVerticalFrame);

	gInterpreter->Declare(
		TString::Format("volatile int &g_menu_action = *((volatile int*)%lu);",
			(unsigned long)&g_menu_action).Data());
	gInterpreter->Declare("void HandleMenuSlot(Int_t id) { g_menu_action = (int)id; }");
	gInterpreter->Declare(
		TString::Format("volatile bool &g_redraw = *((volatile bool*)%lu);",
			(unsigned long)&g_redraw).Data());
	gInterpreter->Declare("void HandleParamChange() { g_redraw = true; }");

	TGPopupMenu *menu_file = new TGPopupMenu(gClient->GetRoot());
	menu_file->AddEntry("&Open...", 1);
	menu_file->AddSeparator();
	menu_file->AddEntry("&Quit", 2);
	menu_file->Connect("Activated(Int_t)", nullptr, nullptr, "HandleMenuSlot(Int_t)");

	TGMenuBar *menu_bar = new TGMenuBar(g_ctx.main_frame, 1, 1, kHorizontalFrame);
	menu_bar->AddPopup("&File", menu_file, new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 0, 0, 0));
	g_ctx.main_frame->AddFrame(menu_bar, new TGLayoutHints(kLHintsTop | kLHintsExpandX));

	TGCompositeFrame *ctrl_frame = new TGHorizontalFrame(g_ctx.main_frame, main_w, 40);
	g_ctx.main_frame->AddFrame(ctrl_frame, new TGLayoutHints(kLHintsTop | kLHintsExpandX));

	g_ctx.btn_all = new TGTextButton(ctrl_frame, "All");
	g_ctx.btn_all->SetState(kButtonDown);
	g_ctx.btn_all->SetCommand("g_menu_action = 7;");
	ctrl_frame->AddFrame(g_ctx.btn_all, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 5, 2));

	g_ctx.btn_14O = new TGTextButton(ctrl_frame, "14O");
	g_ctx.btn_14O->SetCommand("g_menu_action = 4;");
	ctrl_frame->AddFrame(g_ctx.btn_14O, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 2));

	g_ctx.btn_13N = new TGTextButton(ctrl_frame, "13N");
	g_ctx.btn_13N->SetCommand("g_menu_action = 5;");
	ctrl_frame->AddFrame(g_ctx.btn_13N, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 2));

	g_ctx.btn_12C = new TGTextButton(ctrl_frame, "12C");
	g_ctx.btn_12C->SetCommand("g_menu_action = 6;");
	ctrl_frame->AddFrame(g_ctx.btn_12C, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 2));

	ctrl_frame->AddFrame(new TGLabel(ctrl_frame, "  Run:"),
		new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 10, 2));

	g_ctx.entry_run_min = new TGNumberEntry(ctrl_frame, 0, 5, -1,
		TGNumberFormat::kNESInteger, TGNumberFormat::kNEAPositive,
		TGNumberFormat::kNELLimitMinMax, 1, 9999);
	ctrl_frame->AddFrame(g_ctx.entry_run_min, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 0, 2));

	ctrl_frame->AddFrame(new TGLabel(ctrl_frame, "-"),
		new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 2));

	g_ctx.entry_run_max = new TGNumberEntry(ctrl_frame, 100, 5, -1,
		TGNumberFormat::kNESInteger, TGNumberFormat::kNEAPositive,
		TGNumberFormat::kNELLimitMinMax, 1, 9999);
	ctrl_frame->AddFrame(g_ctx.entry_run_max, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 10));

	g_ctx.btn_draw = new TGTextButton(ctrl_frame, "Draw");
	g_ctx.btn_draw->SetCommand("g_menu_action = 3;");
	ctrl_frame->AddFrame(g_ctx.btn_draw, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 10, 5));

	ctrl_frame->AddFrame(new TGLabel(ctrl_frame, "  n1:"),
		new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 10, 2));

	g_ctx.entry_n1 = new TGNumberEntry(ctrl_frame, 300, 5, -1,
		TGNumberFormat::kNESInteger, TGNumberFormat::kNEAPositive,
		TGNumberFormat::kNELLimitMinMax, 10, 10000);
	g_ctx.entry_n1->Connect("ValueSet(Long_t)", nullptr, nullptr, "HandleParamChange()");
	ctrl_frame->AddFrame(g_ctx.entry_n1, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 0, 5));

	ctrl_frame->AddFrame(new TGLabel(ctrl_frame, "n2:"),
		new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 5, 2));

	g_ctx.entry_n2 = new TGNumberEntry(ctrl_frame, 300, 5, -1,
		TGNumberFormat::kNESInteger, TGNumberFormat::kNEAPositive,
		TGNumberFormat::kNELLimitMinMax, 10, 10000);
	g_ctx.entry_n2->Connect("ValueSet(Long_t)", nullptr, nullptr, "HandleParamChange()");
	ctrl_frame->AddFrame(g_ctx.entry_n2, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 0, 5));

	ctrl_frame->AddFrame(new TGLabel(ctrl_frame, "x1:"),
		new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 10, 2));

	g_ctx.entry_x1 = new TGNumberEntry(ctrl_frame, -0.5, 6, -1,
		TGNumberFormat::kNESReal, TGNumberFormat::kNEAAnyNumber,
		TGNumberFormat::kNELLimitMinMax, -30, 30);
	g_ctx.entry_x1->Connect("ValueSet(Long_t)", nullptr, nullptr, "HandleParamChange()");
	ctrl_frame->AddFrame(g_ctx.entry_x1, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 0, 5));

	ctrl_frame->AddFrame(new TGLabel(ctrl_frame, "x2:"),
		new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 5, 2));

	g_ctx.entry_x2 = new TGNumberEntry(ctrl_frame, 0.5, 6, -1,
		TGNumberFormat::kNESReal, TGNumberFormat::kNEAAnyNumber,
		TGNumberFormat::kNELLimitMinMax, -30, 30);
	g_ctx.entry_x2->Connect("ValueSet(Long_t)", nullptr, nullptr, "HandleParamChange()");
	ctrl_frame->AddFrame(g_ctx.entry_x2, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 0, 5));

	g_ctx.main_frame->AddFrame(ctrl_frame, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 2, 2, 2, 0));

	TRootEmbeddedCanvas *ec_main = new TRootEmbeddedCanvas("ec_main", g_ctx.main_frame, main_w, main_h);
	g_ctx.main_frame->AddFrame(ec_main, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
	g_ctx.analysis_canvas = ec_main->GetCanvas();
	g_ctx.analysis_canvas->SetName("c_Analysis");

	TCanvas *c1 = new TCanvas("c_He1", "He1 PID", main_w, main_h);
	g_ctx.pc_He1.canvas = c1;
	c1->SetWindowPosition(50 + main_w + 20, 50);

	TCanvas *c2 = new TCanvas("c_He2", "He2 PID", main_w, main_h);
	g_ctx.pc_He2.canvas = c2;
	c2->SetWindowPosition(50 + main_w + 20, 50 + main_h + 60);

	TCanvas *c3 = new TCanvas("c_He3", "He3 PID", main_w, main_h);
	g_ctx.pc_He3.canvas = c3;
	c3->SetWindowPosition(50 + main_w + 20 + main_w + 20, 50);

	g_ctx.status_bar = new TGStatusBar(g_ctx.main_frame, 100, 20);
	g_ctx.main_frame->AddFrame(g_ctx.status_bar, new TGLayoutHints(kLHintsBottom | kLHintsExpandX));

	g_ctx.main_frame->SetWindowName("Three-Alpha Analysis GUI");
	g_ctx.main_frame->MapSubwindows();
	g_ctx.main_frame->Resize(main_w, main_h);
	g_ctx.main_frame->MapWindow();

	g_ctx.status_bar->SetText("Ready. Select run range and click Draw.");

	while (true) {
		gSystem->DispatchOneEvent(kFALSE);

		if (g_menu_action == 1) {
			OnFileOpen();
		} else if (g_menu_action == 2) {
			break;
		} else if (g_menu_action) {
			ProcessMenu(g_menu_action);
		}

		g_menu_action = 0;

		if (g_redraw) {
			g_redraw = false;
			if (!g_ctx.current_file.empty()) {
				OnDraw();
			} else {
				g_ctx.status_bar->SetText(
					TString::Format("Beam: %s. Click Draw to load data.", BeamLabel()));
			}
		}
	}

	return 0;
}