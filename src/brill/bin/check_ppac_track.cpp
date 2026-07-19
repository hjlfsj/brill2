#include "include/event/ppac/ppac_track.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TPaveText.h>
#include <TAxis.h>
#include <TText.h>
#include <TApplication.h>
#include <TSystem.h>
#include <TPad.h>
#include <TRootEmbeddedCanvas.h>
#include <TROOT.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TLine.h>
#include <TEllipse.h>
#include <TList.h>

#include <TGClient.h>
#include <TGFrame.h>
#include <TGLabel.h>
#include <TGNumberEntry.h>
#include <TGButton.h>
#include <TGLayout.h>
#include <TStyle.h>

#include <filesystem>
#include <iostream>
#include <string>
#include <vector>
#include <cstdio>
#include <sys/select.h>

#include "external/cxxopts.hpp"
#include "include/config.h"
#include "include/event/ppac/ppac_track_event.h"
#include "include/utils.h"

void PrintUsage(const cxxopts::Options &options) {
	std::cout << options.help() << "\n";
}

static constexpr double kZMin = -1000.0;
static constexpr double kZMax = 500.0;
static constexpr double kPosMin = -25.0;
static constexpr double kPosMax = 25.0;

struct CheckContext {
	TCanvas *canvas = nullptr;
	TPad *pad_zx = nullptr;
	TPad *pad_target = nullptr;
	TPad *pad_zy = nullptr;
	TPad *pad_info = nullptr;
	TMultiGraph *mg_zx = nullptr;
	TMultiGraph *mg_zy = nullptr;
	TGraph *g_target = nullptr;
	TEllipse *target_circle = nullptr;
	TPaveText *info_text = nullptr;
	TH2F *frame_zx = nullptr;
	TH2F *frame_zy = nullptr;
	TH2F *frame_target = nullptr;
	TTree *track_tree = nullptr;
	brill::PpacTrackEvent *track_event = nullptr;
	brill::PpacOffsetParams *offset = nullptr;
	brill::PpacConfig *ppac_config = nullptr;
	TGNumberEntry *entry_number = nullptr;
	TGNumberEntry *entry_x_used = nullptr;
	TGNumberEntry *entry_y_used = nullptr;
	long long *current_entry = nullptr;
	long long total = 0;
	std::vector<TGraph*> zx_graphs;
	std::vector<TGraph*> zy_graphs;
	std::vector<TGraph*> track_zx_graphs;
	std::vector<TGraph*> track_zy_graphs;
};

static CheckContext g_ctx;

volatile bool g_check_ppac_read = false;
volatile bool g_check_ppac_next = false;
volatile bool g_check_ppac_add = false;
volatile bool g_check_ppac_clear = false;

void DrawStaticElements() {
	if (!g_ctx.frame_zx) {
		g_ctx.pad_zx->cd();
		gStyle->SetOptStat(0);
		g_ctx.frame_zx = new TH2F("frame_zx", "Z-X Track;z (mm);x (mm)",
			1500, kZMin, kZMax, 100, kPosMin, kPosMax);
		g_ctx.frame_zx->SetDirectory(nullptr);
		g_ctx.frame_zx->Draw();
		for (int i = 0; i < 3; ++i) {
			TLine *line = new TLine(g_ctx.ppac_config->z_x_mm[i], kPosMin,
				g_ctx.ppac_config->z_x_mm[i], kPosMax);
			line->SetLineColor(kRed);
			line->SetLineStyle(2);
			line->Draw();
		}
		TLine *target_line_zx = new TLine(0, kPosMin, 0, kPosMax);
		target_line_zx->SetLineColor(kGreen + 2);
		target_line_zx->SetLineStyle(2);
		target_line_zx->Draw();
	}

	if (!g_ctx.frame_zy) {
		g_ctx.pad_zy->cd();
		gStyle->SetOptStat(0);
		g_ctx.frame_zy = new TH2F("frame_zy", "Z-Y Track;z (mm);y (mm)",
			1500, kZMin, kZMax, 100, kPosMin, kPosMax);
		g_ctx.frame_zy->SetDirectory(nullptr);
		g_ctx.frame_zy->Draw();
		for (int i = 0; i < 3; ++i) {
			TLine *line = new TLine(g_ctx.ppac_config->z_y_mm[i], kPosMin,
				g_ctx.ppac_config->z_y_mm[i], kPosMax);
			line->SetLineColor(kRed);
			line->SetLineStyle(2);
			line->Draw();
		}
		TLine *target_line_zy = new TLine(0, kPosMin, 0, kPosMax);
		target_line_zy->SetLineColor(kGreen + 2);
		target_line_zy->SetLineStyle(2);
		target_line_zy->Draw();
	}

	if (!g_ctx.frame_target) {
		g_ctx.pad_target->cd();
		gStyle->SetOptStat(0);
		g_ctx.frame_target = new TH2F("frame_target", "Target;target_x (mm);target_y (mm)",
			100, -25, 25, 100, -25, 25);
		g_ctx.frame_target->SetDirectory(nullptr);
		g_ctx.frame_target->Draw();
		TEllipse *circle = new TEllipse(0, 0, 15);
		circle->SetLineColor(kRed);
		circle->SetLineWidth(2);
		circle->SetFillStyle(0);
		circle->Draw();
	}
}

void AddEventToCanvas(const brill::PpacTrackEvent &track) {
	int n_points = 0;
	double zx_points[3], zx_vals[3];

	for (int i = 0; i < 3; ++i) {
		if (!std::isnan(track.ppac_x[i])) {
			zx_vals[n_points] = g_ctx.ppac_config->z_x_mm[i];
			zx_points[n_points] = track.ppac_x[i];
			++n_points;
		}
	}
	if (n_points > 0) {
		TGraph *g_pt_zx = new TGraph(n_points, zx_vals, zx_points);
		g_pt_zx->SetMarkerStyle(20);
		g_pt_zx->SetMarkerColor(kBlue);
		g_pt_zx->SetMarkerSize(1.2);
		g_ctx.mg_zx->Add(g_pt_zx);
		g_ctx.zx_graphs.push_back(g_pt_zx);
	}

	if (track.valid) {
		TGraph *g_track_zx = new TGraph(2);
		g_track_zx->SetPoint(0, kZMin, track.target_x + track.dir_x * kZMin);
		g_track_zx->SetPoint(1, kZMax, track.target_x + track.dir_x * kZMax);
		g_track_zx->SetLineColor(kBlue);
		g_track_zx->SetLineWidth(2);
		g_ctx.mg_zx->Add(g_track_zx);
		g_ctx.track_zx_graphs.push_back(g_track_zx);

		TGraph *g_target_zx = new TGraph(1);
		g_target_zx->SetPoint(0, 0, track.target_x);
		g_target_zx->SetMarkerStyle(29);
		g_target_zx->SetMarkerColor(kRed);
		g_target_zx->SetMarkerSize(2.0);
		g_ctx.mg_zx->Add(g_target_zx);
		g_ctx.zx_graphs.push_back(g_target_zx);
	}

	int n_points_y = 0;
	double zy_vals_arr[3], zy_points_arr[3];
	for (int i = 0; i < 3; ++i) {
		if (!std::isnan(track.ppac_y[i])) {
			zy_vals_arr[n_points_y] = g_ctx.ppac_config->z_y_mm[i];
			zy_points_arr[n_points_y] = track.ppac_y[i];
			++n_points_y;
		}
	}
	if (n_points_y > 0) {
		TGraph *g_pt_zy = new TGraph(n_points_y, zy_vals_arr, zy_points_arr);
		g_pt_zy->SetMarkerStyle(20);
		g_pt_zy->SetMarkerColor(kBlue);
		g_pt_zy->SetMarkerSize(1.2);
		g_ctx.mg_zy->Add(g_pt_zy);
		g_ctx.zy_graphs.push_back(g_pt_zy);
	}

	if (track.valid) {
		TGraph *g_track_zy = new TGraph(2);
		g_track_zy->SetPoint(0, kZMin, track.target_y + track.dir_y * kZMin);
		g_track_zy->SetPoint(1, kZMax, track.target_y + track.dir_y * kZMax);
		g_track_zy->SetLineColor(kBlue);
		g_track_zy->SetLineWidth(2);
		g_ctx.mg_zy->Add(g_track_zy);
		g_ctx.track_zy_graphs.push_back(g_track_zy);

		TGraph *g_target_zy = new TGraph(1);
		g_target_zy->SetPoint(0, 0, track.target_y);
		g_target_zy->SetMarkerStyle(29);
		g_target_zy->SetMarkerColor(kRed);
		g_target_zy->SetMarkerSize(2.0);
		g_ctx.mg_zy->Add(g_target_zy);
		g_ctx.zy_graphs.push_back(g_target_zy);
	}

	if (track.valid) {
		int n_target = g_ctx.g_target->GetN();
		g_ctx.g_target->SetPoint(n_target, track.target_x, track.target_y);
	}
}

void UpdateInfoDisplay(const brill::PpacTrackEvent &track, long long entry_num) {
	g_ctx.pad_info->cd();
	g_ctx.info_text->Clear();
	g_ctx.info_text->AddText(TString::Format("Entry: %lld / %lld", entry_num, g_ctx.total));
	g_ctx.info_text->AddText(" ");
	g_ctx.info_text->AddText(TString::Format("Track valid: %d", track.valid));
	g_ctx.info_text->AddText(TString::Format("x_used_ppac: %d  y_used_ppac: %d",
		track.x_used_ppac, track.y_used_ppac));
	g_ctx.info_text->AddText(TString::Format("target_x: %.3f  target_y: %.3f",
		track.target_x, track.target_y));
	g_ctx.info_text->AddText(TString::Format("dir_x: %.5f  dir_y: %.5f",
		track.dir_x, track.dir_y));
	g_ctx.info_text->AddText(" ");
	g_ctx.info_text->AddText("--- PPAC Positions ---");
	for (int i = 0; i < 3; ++i) {
		if (!std::isnan(track.ppac_x[i])) {
			g_ctx.info_text->AddText(TString::Format("PPAC%d: x=%.3f  y=%.3f",
				i + 1, track.ppac_x[i], track.ppac_y[i]));
		} else {
			g_ctx.info_text->AddText(TString::Format("PPAC%d: x=N/A  y=N/A", i + 1));
		}
	}
	g_ctx.info_text->Draw();
}

void PrintEvent(const brill::PpacTrackEvent &track, long long entry_num) {
	printf("\n=== Entry %lld/%lld ===\n", entry_num, g_ctx.total);
	printf("Track valid: %d  x_used_ppac: %d  y_used_ppac: %d\n",
		track.valid, track.x_used_ppac, track.y_used_ppac);
	printf("target_x: %.3f  target_y: %.3f  dir_x: %.5f  dir_y: %.5f\n",
		track.target_x, track.target_y, track.dir_x, track.dir_y);
	printf("PPAC positions:\n");
	for (int i = 0; i < 3; ++i) {
		printf("  PPAC%d: x=%.3f  y=%.3f\n", i + 1, track.ppac_x[i], track.ppac_y[i]);
	}
}

bool PassFilter() {
	int x_filter = g_ctx.entry_x_used->GetIntNumber();
	int y_filter = g_ctx.entry_y_used->GetIntNumber();
	bool x_ok = (x_filter == 0) || ((int)g_ctx.track_event->x_used_ppac == x_filter);
	bool y_ok = (y_filter == 0) || ((int)g_ctx.track_event->y_used_ppac == y_filter);
	return x_ok && y_ok;
}

long long FindNextMatching(long long start) {
	while (start < g_ctx.total) {
		g_ctx.track_tree->GetEntry(start);
		if (PassFilter()) return start;
		++start;
	}
	return -1;
}

void ClearCanvas() {
	for (auto *g : g_ctx.zx_graphs) {
		TList *list = g_ctx.mg_zx->GetListOfGraphs();
		if (list) list->Remove(g);
		delete g;
	}
	for (auto *g : g_ctx.zy_graphs) {
		TList *list = g_ctx.mg_zy->GetListOfGraphs();
		if (list) list->Remove(g);
		delete g;
	}
	for (auto *g : g_ctx.track_zx_graphs) {
		TList *list = g_ctx.mg_zx->GetListOfGraphs();
		if (list) list->Remove(g);
		delete g;
	}
	for (auto *g : g_ctx.track_zy_graphs) {
		TList *list = g_ctx.mg_zy->GetListOfGraphs();
		if (list) list->Remove(g);
		delete g;
	}
	g_ctx.zx_graphs.clear();
	g_ctx.zy_graphs.clear();
	g_ctx.track_zx_graphs.clear();
	g_ctx.track_zy_graphs.clear();
	g_ctx.g_target->Set(0);
}

void ProcessRead() {
	long long entry = g_ctx.entry_number->GetIntNumber();
	if (entry < 0) entry = 0;
	if (entry >= g_ctx.total) entry = g_ctx.total - 1;

	long long found = FindNextMatching(entry);
	if (found < 0) {
		printf("\nNo more matching events.\n");
		return;
	}

	ClearCanvas();
	g_ctx.track_tree->GetEntry(found);

	*g_ctx.current_entry = found;
	g_ctx.entry_number->SetIntNumber(found);
	AddEventToCanvas(*g_ctx.track_event);
	UpdateInfoDisplay(*g_ctx.track_event, found);
	PrintEvent(*g_ctx.track_event, found);

	g_ctx.pad_zx->cd();
	{
		TList *list = g_ctx.mg_zx->GetListOfGraphs();
		if (list && list->GetSize() > 0) g_ctx.mg_zx->Draw("LP");
	}
	g_ctx.pad_zy->cd();
	{
		TList *list = g_ctx.mg_zy->GetListOfGraphs();
		if (list && list->GetSize() > 0) g_ctx.mg_zy->Draw("LP");
	}
	g_ctx.pad_target->cd();
	if (g_ctx.g_target->GetN() > 0) g_ctx.g_target->Draw("P same");

	g_ctx.canvas->Modified();
	g_ctx.canvas->Update();
}

void ProcessNext() {
	long long start = *g_ctx.current_entry + 1;
	long long found = FindNextMatching(start);
	if (found < 0) {
		printf("\nNo more matching events.\n");
		return;
	}

	ClearCanvas();
	g_ctx.track_tree->GetEntry(found);

	*g_ctx.current_entry = found;
	g_ctx.entry_number->SetIntNumber(found);
	AddEventToCanvas(*g_ctx.track_event);
	UpdateInfoDisplay(*g_ctx.track_event, found);
	PrintEvent(*g_ctx.track_event, found);

	g_ctx.pad_zx->cd();
	{
		TList *list = g_ctx.mg_zx->GetListOfGraphs();
		if (list && list->GetSize() > 0) g_ctx.mg_zx->Draw("LP");
	}
	g_ctx.pad_zy->cd();
	{
		TList *list = g_ctx.mg_zy->GetListOfGraphs();
		if (list && list->GetSize() > 0) g_ctx.mg_zy->Draw("LP");
	}
	g_ctx.pad_target->cd();
	if (g_ctx.g_target->GetN() > 0) g_ctx.g_target->Draw("P same");

	g_ctx.canvas->Modified();
	g_ctx.canvas->Update();
}

void ProcessAdd() {
	long long start = *g_ctx.current_entry + 1;
	long long found = FindNextMatching(start);
	if (found < 0) {
		printf("\nNo more matching events.\n");
		return;
	}

	g_ctx.track_tree->GetEntry(found);

	*g_ctx.current_entry = found;
	g_ctx.entry_number->SetIntNumber(found);
	AddEventToCanvas(*g_ctx.track_event);
	UpdateInfoDisplay(*g_ctx.track_event, found);
	PrintEvent(*g_ctx.track_event, found);

	g_ctx.pad_zx->cd();
	{
		TList *list = g_ctx.mg_zx->GetListOfGraphs();
		if (list && list->GetSize() > 0) g_ctx.mg_zx->Draw("LP");
	}
	g_ctx.pad_zy->cd();
	{
		TList *list = g_ctx.mg_zy->GetListOfGraphs();
		if (list && list->GetSize() > 0) g_ctx.mg_zy->Draw("LP");
	}
	g_ctx.pad_target->cd();
	if (g_ctx.g_target->GetN() > 0) g_ctx.g_target->Draw("P same");

	g_ctx.canvas->Modified();
	g_ctx.canvas->Update();
}

int main(int argc, char **argv) {
	cxxopts::Options options("check_ppac_track", "Interactive PPAC track checker.");
	options.add_options()
		("h,help", "Print help information.")
		("r,run", "Run number.", cxxopts::value<int>(), "run")
		("t,trigger", "Trigger type.", cxxopts::value<std::string>()->default_value("main"), "trigger")
		("c,config", "Config file path.",
			cxxopts::value<std::string>()->default_value("config.toml"), "file");

	auto result = options.parse(argc, argv);
	if (result.count("help")) { PrintUsage(options); return 0; }
	if (!result.count("run")) {
		std::cerr << "Error: Missing required option --run.\n";
		PrintUsage(options); return 1;
	}

	brill::AppConfig config;
	if (brill::LoadConfig(result["config"].as<std::string>(), config)) return 1;
	config.trigger = result["trigger"].as<std::string>();
	const int run = result["run"].as<int>();

	std::string track_dir = brill::JoinPath(config.workspace, config.paths.track);
	std::string normalize_dir = brill::JoinPath(config.workspace, config.paths.normalize);

	std::string prefix = (config.trigger == "main") ? "" : "t1_";

	TString track_path = TString::Format("%s/ppac_%s%04d.root",
		track_dir.c_str(), prefix.c_str(), run);
	TString offset_path = TString::Format("%s/ppac_offset_%04d.txt",
		normalize_dir.c_str(), run);

	std::cout << "Track file:  " << track_path << "\n";
	std::cout << "Offset file: " << offset_path << "\n";

	if (!std::filesystem::exists(track_path.Data())) {
		std::cerr << "Error: Track file not found: " << track_path << "\n";
		return 1;
	}

	TFile *track_file = new TFile(track_path, "read");
	TTree *track_tree = (TTree*)track_file->Get("tree");
	if (!track_tree) { std::cerr << "Error: Get track tree failed.\n"; return 1; }

	brill::PpacTrackEvent track_event;
	brill::SetupInput(track_tree, track_event);

	brill::PpacOffsetParams offset[3];
	if (brill::ReadPpacOffsetParams(offset_path.Data(), offset)) return 1;

	long long total = track_tree->GetEntries();
	printf("Total events: %lld\n", total);

	TApplication app("check_ppac_track", nullptr, nullptr);

	gStyle->SetPalette(kRainBow);

	gInterpreter->Declare("#include <cstdint>");
	gInterpreter->Declare(
		TString::Format("volatile bool &g_check_ppac_read = *((volatile bool*)%lu);",
			(unsigned long)&g_check_ppac_read).Data()
	);
	gInterpreter->Declare(
		TString::Format("volatile bool &g_check_ppac_next = *((volatile bool*)%lu);",
			(unsigned long)&g_check_ppac_next).Data()
	);
	gInterpreter->Declare(
		TString::Format("volatile bool &g_check_ppac_add = *((volatile bool*)%lu);",
			(unsigned long)&g_check_ppac_add).Data()
	);
	gInterpreter->Declare(
		TString::Format("volatile bool &g_check_ppac_clear = *((volatile bool*)%lu);",
			(unsigned long)&g_check_ppac_clear).Data()
	);

	TGMainFrame *main_frame = new TGMainFrame(gClient->GetRoot(), 1400, 950);
	TString title = TString::Format("check_ppac_track - Run %d, Trigger %s",
		run, config.trigger.c_str());
	main_frame->SetWindowName(title.Data());

	TRootEmbeddedCanvas *embed = new TRootEmbeddedCanvas("embed", main_frame, 1400, 800);
	main_frame->AddFrame(embed, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));
	TCanvas *canvas = embed->GetCanvas();

	TPad *pad_zx = new TPad("pad_zx", "Z-X Track", 0.0, 0.5, 0.72, 1.0);
	TPad *pad_target = new TPad("pad_target", "Target", 0.72, 0.5, 1.0, 1.0);
	TPad *pad_zy = new TPad("pad_zy", "Z-Y Track", 0.0, 0.0, 0.72, 0.5);
	TPad *pad_info = new TPad("pad_info", "PPAC Info", 0.72, 0.0, 1.0, 0.5);
	pad_zx->Draw();
	pad_target->Draw();
	pad_zy->Draw();
	pad_info->Draw();

	pad_zx->SetGrid();
	pad_zy->SetGrid();

	TMultiGraph *mg_zx = new TMultiGraph();
	TMultiGraph *mg_zy = new TMultiGraph();
	TGraph *g_target = new TGraph();
	g_target->SetName("g_target");
	g_target->SetMarkerStyle(20);
	g_target->SetMarkerColor(kBlue);
	g_target->SetMarkerSize(0.8);

	TPaveText *info_text = new TPaveText(0.05, 0.05, 0.95, 0.95, "NB");
	info_text->SetTextSize(0.035);
	info_text->SetFillColor(kWhite);
	info_text->SetBorderSize(1);

	long long current_entry = 0;

	g_ctx.canvas = canvas;
	g_ctx.pad_zx = pad_zx;
	g_ctx.pad_target = pad_target;
	g_ctx.pad_zy = pad_zy;
	g_ctx.pad_info = pad_info;
	g_ctx.mg_zx = mg_zx;
	g_ctx.mg_zy = mg_zy;
	g_ctx.g_target = g_target;
	g_ctx.info_text = info_text;
	g_ctx.track_tree = track_tree;
	g_ctx.track_event = &track_event;
	g_ctx.offset = offset;
	g_ctx.ppac_config = &config.ppac;
	g_ctx.current_entry = &current_entry;
	g_ctx.total = total;

	DrawStaticElements();

	TGHorizontalFrame *ctrl_frame = new TGHorizontalFrame(main_frame, 1400, 40);

	ctrl_frame->AddFrame(new TGLabel(ctrl_frame, "Entry:"),
		new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 5, 2, 5, 5));
	TGNumberEntry *entry_number = new TGNumberEntry(ctrl_frame, 0, 8, -1,
		TGNumberFormat::kNESInteger,
		TGNumberFormat::kNEANonNegative,
		TGNumberFormat::kNELLimitMinMax, 0, total - 1);
	ctrl_frame->AddFrame(entry_number,
		new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 5, 5, 5));
	g_ctx.entry_number = entry_number;

	ctrl_frame->AddFrame(new TGLabel(ctrl_frame, "  x_used_ppac:"),
		new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 8, 2, 5, 5));
	TGNumberEntry *entry_x_used = new TGNumberEntry(ctrl_frame, 0, 3, -1,
		TGNumberFormat::kNESInteger,
		TGNumberFormat::kNEANonNegative,
		TGNumberFormat::kNELLimitMinMax, 0, 7);
	ctrl_frame->AddFrame(entry_x_used,
		new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 5, 5, 5));
	g_ctx.entry_x_used = entry_x_used;

	ctrl_frame->AddFrame(new TGLabel(ctrl_frame, "  y_used_ppac:"),
		new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 8, 2, 5, 5));
	TGNumberEntry *entry_y_used = new TGNumberEntry(ctrl_frame, 0, 3, -1,
		TGNumberFormat::kNESInteger,
		TGNumberFormat::kNEANonNegative,
		TGNumberFormat::kNELLimitMinMax, 0, 7);
	ctrl_frame->AddFrame(entry_y_used,
		new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 5, 5, 5));
	g_ctx.entry_y_used = entry_y_used;

	TGTextButton *btn_read = new TGTextButton(ctrl_frame, "Read");
	btn_read->SetCommand("g_check_ppac_read = true;");
	ctrl_frame->AddFrame(btn_read,
		new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 15, 5, 5, 5));

	TGTextButton *btn_next = new TGTextButton(ctrl_frame, "Next");
	btn_next->SetCommand("g_check_ppac_next = true;");
	ctrl_frame->AddFrame(btn_next,
		new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 5, 5, 5, 5));

	TGTextButton *btn_add = new TGTextButton(ctrl_frame, "Add");
	btn_add->SetCommand("g_check_ppac_add = true;");
	ctrl_frame->AddFrame(btn_add,
		new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 5, 5, 5, 5));

	TGTextButton *btn_clear = new TGTextButton(ctrl_frame, "Clear");
	btn_clear->SetCommand("g_check_ppac_clear = true;");
	ctrl_frame->AddFrame(btn_clear,
		new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 5, 5, 5, 5));

	main_frame->AddFrame(ctrl_frame,
		new TGLayoutHints(kLHintsBottom | kLHintsExpandX, 2, 2, 2, 2));

	main_frame->SetMWMHints(kMWMDecorAll, kMWMFuncAll, kMWMInputModeless);
	main_frame->MapSubwindows();
	main_frame->Resize(main_frame->GetDefaultSize());
	main_frame->MapWindow();

	printf("\n");
	printf("=== check_ppac_track ===\n");
	printf("Left-top:    Z-X track view\n");
	printf("Right-top:   target_x vs target_y scatter\n");
	printf("Left-bottom: Z-Y track view\n");
	printf("Right-bottom: PPAC event details\n");
	printf("Controls:\n");
	printf("  Entry:       set entry number\n");
	printf("  x_used_ppac: filter by x_used_ppac (0=all, 1-7=specific)\n");
	printf("  y_used_ppac: filter by y_used_ppac (0=all, 1-7=specific)\n");
	printf("  Read:  clear + read next matching event\n");
	printf("  Next:  clear + draw next matching event only\n");
	printf("  Add:   overlay next matching event\n");
	printf("  Clear: clear all drawings\n");
	printf("========================\n\n");

	while (true) {
		gSystem->ProcessEvents();
		if (g_check_ppac_read) {
			g_check_ppac_read = false;
			ProcessRead();
		}
		if (g_check_ppac_next) {
			g_check_ppac_next = false;
			ProcessNext();
		}
		if (g_check_ppac_add) {
			g_check_ppac_add = false;
			ProcessAdd();
		}
		if (g_check_ppac_clear) {
			g_check_ppac_clear = false;
			ClearCanvas();
		}
		fd_set rfds;
		FD_ZERO(&rfds);
		struct timeval tv = {0, 50000};
		select(0, &rfds, nullptr, nullptr, &tv);
	}

	return 0;
}