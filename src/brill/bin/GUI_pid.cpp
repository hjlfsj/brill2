#include "include/config.h"
#include "include/d_6Li/extract.h"
#include "include/event/ingot/silicon_event.h"
#include "include/event/t0/dssd_match_event.h"
#include "include/utils.h"
#include "external/cxxopts.hpp"

#include <TApplication.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGButton.h>
#include <TGClient.h>
#include <TGComboBox.h>
#include <TGFrame.h>
#include <TGLayout.h>
#include <TGLabel.h>
#include <TGMenu.h>
#include <TGNumberEntry.h>
#include <TGStatusBar.h>
#include <TGTextEntry.h>
#include <TH2D.h>
#include <TPad.h>
#include <TRootEmbeddedCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <cstdio>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

struct PidCanvases {
	TRootEmbeddedCanvas *embed = nullptr;
	TCanvas *canvas = nullptr;
	TH2D *h_d1_d2 = nullptr;
	TH2D *h_d2_d3 = nullptr;
	TH2D *h_d3_d4 = nullptr;
	TH2D *h_d4_t0s = nullptr;
};

struct GUIContext {
	TGMainFrame *main_frame = nullptr;
	TGStatusBar *status_bar = nullptr;

	PidCanvases main_cv;
	PidCanvases sec_cv;

	TGComboBox *trigger_combo = nullptr;
	TGNumberEntryField *run_start = nullptr;
	TGNumberEntryField *run_end = nullptr;
	TGNumberEntry *hit_d1 = nullptr;
	TGNumberEntry *hit_d2 = nullptr;
	TGNumberEntry *hit_d3 = nullptr;
	TGNumberEntry *hit_d4 = nullptr;
	TGTextButton *draw_btn = nullptr;

	std::string config_path;
	std::string match_dir;
	std::string ingot_dir;
	brill::D6LiCalibration calib;
	brill::AppConfig config;
};

static GUIContext g_ctx;
static volatile int g_menu_action = 0;

static void RebuildMainHistograms() {
	auto &cv = g_ctx.main_cv;
	if (cv.h_d1_d2) delete cv.h_d1_d2;
	if (cv.h_d2_d3) delete cv.h_d2_d3;
	if (cv.h_d3_d4) delete cv.h_d3_d4;
	if (cv.h_d4_t0s) delete cv.h_d4_t0s;

	cv.h_d1_d2 = new TH2D("h_d1_d2", "PID D1-D2;D2 Energy (MeV);D1 Energy (MeV)",
		1000, 0, 400, 1000, 0, 200);
	cv.h_d1_d2->SetDirectory(0);
	cv.h_d2_d3 = new TH2D("h_d2_d3", "PID D2-D3;D3 Energy (MeV);D2 Energy (MeV)",
		1000, 0, 400, 1000, 0, 400);
	cv.h_d2_d3->SetDirectory(0);
	cv.h_d3_d4 = new TH2D("h_d3_d4", "PID D3-D4;D4 Energy (MeV);D3 Energy (MeV)",
		1000, 0, 300, 1000, 0, 300);
	cv.h_d3_d4->SetDirectory(0);
	cv.h_d4_t0s = new TH2D("h_d4_t0s", "PID D4-T0S;T0S Energy (MeV);D4 Energy (MeV)",
		1000, 0, 300, 1000, 0, 300);
	cv.h_d4_t0s->SetDirectory(0);
}

static void RebuildSecHistograms() {
	auto &cv = g_ctx.sec_cv;
	if (cv.h_d1_d2) delete cv.h_d1_d2;
	if (cv.h_d2_d3) delete cv.h_d2_d3;
	if (cv.h_d3_d4) delete cv.h_d3_d4;
	if (cv.h_d4_t0s) delete cv.h_d4_t0s;

	cv.h_d1_d2 = new TH2D("h_d1_d2_sec", "PID D1-D2 (sec);D2 Energy (MeV);D1 Energy (MeV)",
		1000, 0, 400, 1000, 0, 200);
	cv.h_d1_d2->SetDirectory(0);
	cv.h_d2_d3 = new TH2D("h_d2_d3_sec", "PID D2-D3 (sec);D3 Energy (MeV);D2 Energy (MeV)",
		1000, 0, 400, 1000, 0, 400);
	cv.h_d2_d3->SetDirectory(0);
	cv.h_d3_d4 = new TH2D("h_d3_d4_sec", "PID D3-D4 (sec);D4 Energy (MeV);D3 Energy (MeV)",
		1000, 0, 300, 1000, 0, 300);
	cv.h_d3_d4->SetDirectory(0);
	cv.h_d4_t0s = new TH2D("h_d4_t0s_sec", "PID D4-T0S (sec);T0S Energy (MeV);D4 Energy (MeV)",
		1000, 0, 300, 1000, 0, 300);
	cv.h_d4_t0s->SetDirectory(0);
}

static void DrawMainCanvas() {
	auto &cv = g_ctx.main_cv;
	if (!cv.canvas) return;
	cv.canvas->Clear();
	cv.canvas->Divide(2, 2);
	cv.canvas->cd(1); cv.h_d1_d2->Draw("colz");
	cv.canvas->cd(2); cv.h_d2_d3->Draw("colz");
	cv.canvas->cd(3); cv.h_d3_d4->Draw("colz");
	cv.canvas->cd(4); cv.h_d4_t0s->Draw("colz");
	cv.canvas->Modified();
	cv.canvas->Update();
}

static void DrawSecCanvas() {
	auto &cv = g_ctx.sec_cv;
	if (!cv.canvas) return;
	cv.canvas->Clear();
	cv.canvas->Divide(2, 2);

	int pad_idx = 0;
	TH2D *histos[4] = {cv.h_d1_d2, cv.h_d2_d3, cv.h_d3_d4, cv.h_d4_t0s};

	for (int i = 0; i < 4; ++i) {
		if (histos[i]->GetEntries() == 0) continue;
		pad_idx++;
		cv.canvas->cd(pad_idx);
		histos[i]->Draw("colz");
	}

	cv.canvas->Modified();
	cv.canvas->Update();
}

static void OnDraw() {
	printf("=== PID Draw ===\n");

	int run_start = g_ctx.run_start->GetIntNumber();
	int run_end = g_ctx.run_end->GetIntNumber();
	int hd1 = g_ctx.hit_d1->GetIntNumber();
	int hd2 = g_ctx.hit_d2->GetIntNumber();
	int hd3 = g_ctx.hit_d3->GetIntNumber();
	int hd4 = g_ctx.hit_d4->GetIntNumber();

	int trigger_id = g_ctx.trigger_combo->GetSelected();
	std::string trigger = (trigger_id == 0) ? "main" : "t1";

	printf("  trigger=%s, run=%d-%d, hit=(%d,%d,%d,%d)\n",
		trigger.c_str(), run_start, run_end, hd1, hd2, hd3, hd4);

	std::string trigger_infix = brill::TriggerInfix(trigger);

	RebuildMainHistograms();
	RebuildSecHistograms();

	g_ctx.main_cv.h_d1_d2->SetTitle("PID D1-D2 (hit0);D2 Energy (MeV);D1 Energy (MeV)");
	g_ctx.main_cv.h_d2_d3->SetTitle("PID D2-D3 (hit0);D3 Energy (MeV);D2 Energy (MeV)");
	g_ctx.main_cv.h_d3_d4->SetTitle("PID D3-D4 (hit0);D4 Energy (MeV);D3 Energy (MeV)");
	g_ctx.main_cv.h_d4_t0s->SetTitle("PID D4-T0S (hit0);T0S Energy (MeV);D4 Energy (MeV)");

	const char *d1s = (hd1 >= 2) ? "hit1" : "hit0";
	const char *d2s = (hd2 >= 2) ? "hit1" : "hit0";
	const char *d3s = (hd3 >= 2) ? "hit1" : "hit0";
	const char *d4s = (hd4 >= 2) ? "hit1" : "hit0";

	g_ctx.sec_cv.h_d1_d2->SetTitle(TString::Format(
		"PID D1(%s)-D2(%s);D2 Energy (MeV);D1 Energy (MeV)", d1s, d2s));
	g_ctx.sec_cv.h_d2_d3->SetTitle(TString::Format(
		"PID D2(%s)-D3(%s);D3 Energy (MeV);D2 Energy (MeV)", d2s, d3s));
	g_ctx.sec_cv.h_d3_d4->SetTitle(TString::Format(
		"PID D3(%s)-D4(%s);D4 Energy (MeV);D3 Energy (MeV)", d3s, d4s));
	g_ctx.sec_cv.h_d4_t0s->SetTitle(TString::Format(
		"PID D4(%s)-T0S;T0S Energy (MeV);D4 Energy (MeV)", d4s));

	int main_filled = 0;
	int main_skipped = 0;
	int sec_events = 0;

	for (int run = run_start; run <= run_end; ++run) {
		if (brill::IsJumpRun(g_ctx.config, run)) continue;

		TString d1_path = TString::Format("%s/t0d1_%s%04d.root",
			g_ctx.match_dir.c_str(), trigger_infix.c_str(), run);
		TString d2_path = TString::Format("%s/t0d2_%s%04d.root",
			g_ctx.match_dir.c_str(), trigger_infix.c_str(), run);
		TString d3_path = TString::Format("%s/t0d3_%s%04d.root",
			g_ctx.match_dir.c_str(), trigger_infix.c_str(), run);
		TString d4_path = TString::Format("%s/t0d4_%s%04d.root",
			g_ctx.match_dir.c_str(), trigger_infix.c_str(), run);
		TString t0s_path = TString::Format("%s/t0s_%s%04d.root",
			g_ctx.ingot_dir.c_str(), trigger_infix.c_str(), run);

		if (!std::filesystem::exists(d1_path.Data())) {
			printf("  Skip run %d: d1 file not found\n", run);
			continue;
		}

		TFile f1(d1_path, "read");
		TFile f2(d2_path, "read");
		TFile f3(d3_path, "read");
		TFile f4(d4_path, "read");
		TFile fs(t0s_path, "read");
		TTree *t1 = (TTree*)f1.Get("tree");
		TTree *t2 = (TTree*)f2.Get("tree");
		TTree *t3 = (TTree*)f3.Get("tree");
		TTree *t4 = (TTree*)f4.Get("tree");
		TTree *ts = (TTree*)fs.Get("tree");

		if (!t1) { printf("  Skip run %d: no tree in d1\n", run); continue; }

		Long64_t n = t1->GetEntries();
		Long64_t report_step = (n > 0) ? n / 10 : 1;
		if (report_step < 1) report_step = 1;
		printf("  Run %d: %lld entries\n", run, n);

		brill::DssdMatchEvent d1_ev, d2_ev, d3_ev, d4_ev;
		brill::SiliconEvent t0s_ev;
		brill::SetupInput(t1, d1_ev);
		if (t2) brill::SetupInput(t2, d2_ev);
		if (t3) brill::SetupInput(t3, d3_ev);
		if (t4) brill::SetupInput(t4, d4_ev);
		if (ts) brill::SetupInput(ts, t0s_ev);

		for (Long64_t entry = 0; entry < n; ++entry) {
			t1->GetEntry(entry);
			if (t2) t2->GetEntry(entry);
			if (t3) t3->GetEntry(entry);
			if (t4) t4->GetEntry(entry);
			if (ts) ts->GetEntry(entry);

			if (entry % report_step == 0) {
				printf("\r    Progress: %.0f%%", 100.0 * (double)entry / (double)n);
				fflush(stdout);
			}

			bool d1_ok = (d1_ev.num == hd1);
			bool d2_ok = (d2_ev.num == hd2);
			bool d3_ok = (d3_ev.num == hd3);
			bool d4_ok = (d4_ev.num == hd4);

			if (!d1_ok || !d2_ok || !d3_ok || !d4_ok) {
				main_skipped++;
				continue;
			}

			main_filled++;

			double e1 = brill::CalibrateD6LiEnergy(g_ctx.calib, 0, d1_ev.energy[0]);
			double e2 = brill::CalibrateD6LiEnergy(g_ctx.calib, 1, d2_ev.energy[0]);
			double e3 = brill::CalibrateD6LiEnergy(g_ctx.calib, 2, d3_ev.energy[0]);
			double e4 = brill::CalibrateD6LiEnergy(g_ctx.calib, 3, d4_ev.energy[0]);
			double es = t0s_ev.valid ? brill::CalibrateD6LiEnergy(g_ctx.calib, 4, (double)t0s_ev.energy) : 0.0;

			g_ctx.main_cv.h_d1_d2->Fill(e2, e1);
			g_ctx.main_cv.h_d2_d3->Fill(e3, e2);
			g_ctx.main_cv.h_d3_d4->Fill(e4, e3);
			g_ctx.main_cv.h_d4_t0s->Fill(es, e4);

			bool d1_sec = (d1_ev.num >= 2);
			bool d2_sec = (d2_ev.num >= 2);
			bool d3_sec = (d3_ev.num >= 2);
			bool d4_sec = (d4_ev.num >= 2);

			if (d1_sec || d2_sec || d3_sec || d4_sec) sec_events++;

			double e1s = d1_sec ? brill::CalibrateD6LiEnergy(g_ctx.calib, 0, d1_ev.energy[1]) : e1;
			double e2s = d2_sec ? brill::CalibrateD6LiEnergy(g_ctx.calib, 1, d2_ev.energy[1]) : e2;
			double e3s = d3_sec ? brill::CalibrateD6LiEnergy(g_ctx.calib, 2, d3_ev.energy[1]) : e3;
			double e4s = d4_sec ? brill::CalibrateD6LiEnergy(g_ctx.calib, 3, d4_ev.energy[1]) : e4;

			if (d1_sec) g_ctx.sec_cv.h_d1_d2->Fill(e2s, e1s);
			if (d2_sec) g_ctx.sec_cv.h_d2_d3->Fill(e3s, e2s);
			if (d3_sec) g_ctx.sec_cv.h_d3_d4->Fill(e4s, e3s);
			if (d4_sec) g_ctx.sec_cv.h_d4_t0s->Fill(es, e4s);
		}

		printf("\r    Run %d: 100%% complete\n", run);
		fflush(stdout);

		if (trigger == "main") {
			DrawMainCanvas();
			DrawSecCanvas();
			g_ctx.status_bar->SetText(TString::Format(
				"trigger=%s  run=%d-%d  hit=(%d,%d,%d,%d)  main=%d  sec=%d (reading...)",
				trigger.c_str(), run_start, run, hd1, hd2, hd3, hd4,
				main_filled, sec_events));
			gSystem->ProcessEvents();
		}
	}

	printf("  Main: %d filled, %d skipped (hit_num out of range)\n",
		main_filled, main_skipped);
	printf("  Secondary: %d events with extra hits\n", sec_events);

	DrawMainCanvas();
	DrawSecCanvas();

	g_ctx.status_bar->SetText(TString::Format(
		"trigger=%s  run=%d-%d  hit=(%d,%d,%d,%d)  main=%d  sec=%d",
		trigger.c_str(), run_start, run_end, hd1, hd2, hd3, hd4,
		main_filled, sec_events));

	printf("  Done.\n");
}

int main(int argc, char **argv) {
	cxxopts::Options options("GUI_pid", "GUI for DSSD PID visualization");
	options.add_options()
		("c,config", "Config file path.", cxxopts::value<std::string>()->default_value("config.toml"))
		("h,help", "Print help information.");
	auto result = options.parse(argc, argv);

	if (result.count("help")) {
		std::cout << options.help() << "\n";
		return 0;
	}

	g_ctx.config_path = result["config"].as<std::string>();
	if (brill::LoadConfig(g_ctx.config_path, g_ctx.config)) {
		std::cerr << "Error: Load config failed.\n";
		return 1;
	}
	g_ctx.match_dir = brill::JoinPath(g_ctx.config.workspace, g_ctx.config.paths.match);
	g_ctx.ingot_dir = brill::JoinPath(g_ctx.config.workspace, g_ctx.config.paths.ingot);

	std::string calib_path = TString::Format(
		"%s/t0.txt",
		brill::JoinPath(g_ctx.config.workspace, g_ctx.config.paths.calibration).c_str()
	).Data();
	if (brill::ReadD6LiCalibration(calib_path, g_ctx.calib)) {
		std::cerr << "Error: Read calibration from " << calib_path << " failed.\n";
		return 1;
	}
	printf("Calibration loaded: %s\n", calib_path.c_str());

	TApplication app("GUI_pid", &argc, argv);
	gStyle->SetPalette(kRainBow);

	TString decl_action = TString::Format(
		"volatile int &g_menu_action = *((volatile int*)%lu);",
		(unsigned long)&g_menu_action
	);
	gInterpreter->Declare(decl_action.Data());
	gInterpreter->Declare("void HandleMenuSlot(Int_t id) { g_menu_action = (int)id; }");

	TGMainFrame *main_frame = new TGMainFrame(gClient->GetRoot(), 1200, 900);
	main_frame->SetWindowName("GUI_pid");
	g_ctx.main_frame = main_frame;

	TGPopupMenu *menu_file = new TGPopupMenu(gClient->GetRoot());
	menu_file->AddEntry("&Quit", 2);
	menu_file->Connect("Activated(Int_t)", nullptr, nullptr, "HandleMenuSlot(Int_t)");

	TGMenuBar *menu_bar = new TGMenuBar(main_frame, 1, 1, kHorizontalFrame);
	menu_bar->AddPopup("&File", menu_file, new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 0, 0, 0));
	main_frame->AddFrame(menu_bar, new TGLayoutHints(kLHintsTop | kLHintsExpandX));

	TRootEmbeddedCanvas *embed = new TRootEmbeddedCanvas("embed_main", main_frame, 1200, 600);
	main_frame->AddFrame(embed, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));
	g_ctx.main_cv.embed = embed;
	g_ctx.main_cv.canvas = embed->GetCanvas();

	TGHorizontalFrame *ctrl_frame = new TGHorizontalFrame(main_frame, 1200, 40);
	main_frame->AddFrame(ctrl_frame, new TGLayoutHints(kLHintsBottom | kLHintsExpandX, 2, 2, 2, 2));

	ctrl_frame->AddFrame(new TGLabel(ctrl_frame, "Trigger:"),
		new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 2, 2, 2, 2));

	g_ctx.trigger_combo = new TGComboBox(ctrl_frame);
	g_ctx.trigger_combo->AddEntry("main", 0);
	g_ctx.trigger_combo->AddEntry("t1", 1);
	g_ctx.trigger_combo->Select(0);
	g_ctx.trigger_combo->Resize(80, 22);
	ctrl_frame->AddFrame(g_ctx.trigger_combo,
		new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 2, 10, 2, 2));

	ctrl_frame->AddFrame(new TGLabel(ctrl_frame, "Run:"),
		new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 2, 2, 2, 2));

	g_ctx.run_start = new TGNumberEntryField(ctrl_frame, -1, 57, TGNumberFormat::kNESInteger);
	g_ctx.run_start->Resize(60, 22);
	ctrl_frame->AddFrame(g_ctx.run_start,
		new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 2, 2, 2, 2));

	ctrl_frame->AddFrame(new TGLabel(ctrl_frame, "-"),
		new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 2, 2, 2, 2));

	g_ctx.run_end = new TGNumberEntryField(ctrl_frame, -1, 92, TGNumberFormat::kNESInteger);
	g_ctx.run_end->Resize(60, 22);
	ctrl_frame->AddFrame(g_ctx.run_end,
		new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 2, 10, 2, 2));

	ctrl_frame->AddFrame(new TGLabel(ctrl_frame, "d1_hit:"),
		new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 2, 2, 2, 2));

	g_ctx.hit_d1 = new TGNumberEntry(ctrl_frame, -1, 0, 3, TGNumberFormat::kNESInteger,
		TGNumberFormat::kNEAAnyNumber, TGNumberFormat::kNELLimitMinMax, 0, 7);
	g_ctx.hit_d1->Resize(60, 22);
	ctrl_frame->AddFrame(g_ctx.hit_d1,
		new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 2, 5, 2, 2));

	ctrl_frame->AddFrame(new TGLabel(ctrl_frame, "d2_hit:"),
		new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 2, 2, 2, 2));

	g_ctx.hit_d2 = new TGNumberEntry(ctrl_frame, -1, 0, 3, TGNumberFormat::kNESInteger,
		TGNumberFormat::kNEAAnyNumber, TGNumberFormat::kNELLimitMinMax, 0, 7);
	g_ctx.hit_d2->Resize(60, 22);
	ctrl_frame->AddFrame(g_ctx.hit_d2,
		new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 2, 5, 2, 2));

	ctrl_frame->AddFrame(new TGLabel(ctrl_frame, "d3_hit:"),
		new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 2, 2, 2, 2));

	g_ctx.hit_d3 = new TGNumberEntry(ctrl_frame, -1, 0, 3, TGNumberFormat::kNESInteger,
		TGNumberFormat::kNEAAnyNumber, TGNumberFormat::kNELLimitMinMax, 0, 7);
	g_ctx.hit_d3->Resize(60, 22);
	ctrl_frame->AddFrame(g_ctx.hit_d3,
		new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 2, 5, 2, 2));

	ctrl_frame->AddFrame(new TGLabel(ctrl_frame, "d4_hit:"),
		new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 2, 2, 2, 2));

	g_ctx.hit_d4 = new TGNumberEntry(ctrl_frame, -1, 0, 3, TGNumberFormat::kNESInteger,
		TGNumberFormat::kNEAAnyNumber, TGNumberFormat::kNELLimitMinMax, 0, 7);
	g_ctx.hit_d4->Resize(60, 22);
	ctrl_frame->AddFrame(g_ctx.hit_d4,
		new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 2, 10, 2, 2));

	g_ctx.draw_btn = new TGTextButton(ctrl_frame, "Draw");
	g_ctx.draw_btn->Resize(60, 22);
	g_ctx.draw_btn->SetCommand("g_menu_action = 3;");
	ctrl_frame->AddFrame(g_ctx.draw_btn,
		new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 2, 2, 2, 2));

	g_ctx.sec_cv.canvas = new TCanvas("canvas_sec", "PID Secondary Hits", 1200, 800);

	TGStatusBar *status_bar = new TGStatusBar(main_frame, 1, 1);
	main_frame->AddFrame(status_bar, new TGLayoutHints(kLHintsBottom | kLHintsExpandX));
	g_ctx.status_bar = status_bar;
	status_bar->SetText("Ready. Set trigger/run/hit and click Draw.");

	RebuildMainHistograms();
	RebuildSecHistograms();
	DrawMainCanvas();
	DrawSecCanvas();

	main_frame->MapSubwindows();
	main_frame->Resize(main_frame->GetDefaultSize());
	main_frame->MapWindow();

	while (true) {
		gSystem->DispatchOneEvent(kFALSE);
		if (g_menu_action == 3) {
			OnDraw();
		} else if (g_menu_action == 2) {
			break;
		}
		g_menu_action = 0;
	}

	return 0;
}