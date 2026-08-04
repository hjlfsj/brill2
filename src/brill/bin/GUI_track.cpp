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
#include <TH3.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>

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
#include "include/event/t0/dssd_match_event.h"
#include "include/utils.h"

void PrintUsage(const cxxopts::Options &options) {
	std::cout << options.help() << "\n";
}

static constexpr double kZMin = -650.0;
static constexpr double kZMax = 150.0;
static constexpr double kZMinZoom = 80.0;
static constexpr double kZMaxZoom = 150.0;
static constexpr double kPosMin = -35.0;
static constexpr double kPosMax = 35.0;

struct T0DetectorEntry {
	std::string name;
	TFile *file = nullptr;
	TTree *tree = nullptr;
	brill::DssdMatchEvent *event = nullptr;
	double z_mm = 0.0;
	int color = kGreen;
	int marker_style = 21;
};

struct CheckContext {
	TCanvas *canvas1 = nullptr;
	TCanvas *canvas2 = nullptr;
	TCanvas *canvas3 = nullptr;
	TCanvas *canvas4 = nullptr;
	TPad *pad_zx = nullptr;
	TPad *pad_zy = nullptr;
	TPad *pad_zx_zoom = nullptr;
	TPad *pad_zy_zoom = nullptr;
	TPad *pad_3d = nullptr;
	TPad *pad_t0_info = nullptr;
	TPad *pad_fitting = nullptr;
	TMultiGraph *mg_zx = nullptr;
	TMultiGraph *mg_zy = nullptr;
	TMultiGraph *mg_zx_zoom = nullptr;
	TMultiGraph *mg_zy_zoom = nullptr;
	TPaveText *t0_info_text = nullptr;
	TH2F *frame_zx = nullptr;
	TH2F *frame_zy = nullptr;
	TH2F *frame_zx_zoom = nullptr;
	TH2F *frame_zy_zoom = nullptr;
	TH3F *frame_3d = nullptr;
	TTree *input_tree = nullptr;
	brill::PpacTrackEvent *track_event = nullptr;
	brill::PpacConfig *ppac_config = nullptr;
	TGNumberEntry *entry_number = nullptr;
	TGNumberEntry *entry_x_used = nullptr;
	TGNumberEntry *entry_y_used = nullptr;
	TGNumberEntry *entry_d1_num = nullptr;
	TGNumberEntry *entry_d2_num = nullptr;
	TGNumberEntry *entry_d3_num = nullptr;
	TGNumberEntry *entry_d4_num = nullptr;
	long long *current_entry = nullptr;
	long long total = 0;
	double t0_z_min = 0.0;
	double t0_z_max = 500.0;
	std::vector<TGraph*> zx_graphs;
	std::vector<TGraph*> zy_graphs;
	std::vector<TGraph*> zx_zoom_graphs;
	std::vector<TGraph*> zy_zoom_graphs;
	std::vector<TObject*> objects_3d;
	std::vector<T0DetectorEntry> t0_dets;
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
			800, kZMin, kZMax, 100, kPosMin, kPosMax);
		g_ctx.frame_zx->SetDirectory(nullptr);
		g_ctx.frame_zx->Draw();
		for (int i = 0; i < 3; ++i) {
			TLine *line = new TLine(g_ctx.ppac_config->z_x_mm[i], kPosMin,
				g_ctx.ppac_config->z_x_mm[i], kPosMax);
			line->SetLineColor(kRed);
			line->SetLineStyle(2);
			line->Draw();
			TText *label = new TText(g_ctx.ppac_config->z_x_mm[i], kPosMax - 2.0,
				TString::Format("PPAC%d", i + 1));
			label->SetTextSize(0.025);
			label->SetTextAlign(22);
			label->SetTextColor(kRed);
			label->Draw();
		}
		TLine *target_line_zx = new TLine(0, kPosMin, 0, kPosMax);
		target_line_zx->SetLineColor(kGreen + 2);
		target_line_zx->SetLineStyle(2);
		target_line_zx->Draw();
		{
			TText *label = new TText(0, kPosMax - 2.0, "target");
			label->SetTextSize(0.025);
			label->SetTextAlign(22);
			label->SetTextColor(kGreen + 2);
			label->Draw();
		}
		for (size_t i = 0; i < g_ctx.t0_dets.size(); ++i) {
			auto &det = g_ctx.t0_dets[i];
			TLine *t0_line = new TLine(det.z_mm, kPosMin, det.z_mm, kPosMax);
			t0_line->SetLineColor(det.color);
			t0_line->SetLineStyle(3);
			t0_line->SetLineWidth(3);
			t0_line->Draw();
			TText *label = new TText(det.z_mm, kPosMax - 4.0,
				TString::Format("d%d", (int)i + 1));
			label->SetTextSize(0.025);
			label->SetTextAlign(22);
			label->SetTextColor(det.color);
			label->Draw();
		}
	}

	if (!g_ctx.frame_zy) {
		g_ctx.pad_zy->cd();
		gStyle->SetOptStat(0);
		g_ctx.frame_zy = new TH2F("frame_zy", "Z-Y Track;z (mm);y (mm)",
			800, kZMin, kZMax, 100, kPosMin, kPosMax);
		g_ctx.frame_zy->SetDirectory(nullptr);
		g_ctx.frame_zy->Draw();
		for (int i = 0; i < 3; ++i) {
			TLine *line = new TLine(g_ctx.ppac_config->z_y_mm[i], kPosMin,
				g_ctx.ppac_config->z_y_mm[i], kPosMax);
			line->SetLineColor(kRed);
			line->SetLineStyle(2);
			line->Draw();
			TText *label = new TText(g_ctx.ppac_config->z_y_mm[i], kPosMax - 2.0,
				TString::Format("PPAC%d", i + 1));
			label->SetTextSize(0.025);
			label->SetTextAlign(22);
			label->SetTextColor(kRed);
			label->Draw();
		}
		TLine *target_line_zy = new TLine(0, kPosMin, 0, kPosMax);
		target_line_zy->SetLineColor(kGreen + 2);
		target_line_zy->SetLineStyle(2);
		target_line_zy->Draw();
		{
			TText *label = new TText(0, kPosMax - 2.0, "target");
			label->SetTextSize(0.025);
			label->SetTextAlign(22);
			label->SetTextColor(kGreen + 2);
			label->Draw();
		}
		for (size_t i = 0; i < g_ctx.t0_dets.size(); ++i) {
			auto &det = g_ctx.t0_dets[i];
			TLine *t0_line = new TLine(det.z_mm, kPosMin, det.z_mm, kPosMax);
			t0_line->SetLineColor(det.color);
			t0_line->SetLineStyle(3);
			t0_line->SetLineWidth(3);
			t0_line->Draw();
			TText *label = new TText(det.z_mm, kPosMax - 4.0,
				TString::Format("d%d", (int)i + 1));
			label->SetTextSize(0.025);
			label->SetTextAlign(22);
			label->SetTextColor(det.color);
			label->Draw();
		}
	}
}

void DrawStaticElementsZoom() {
	if (!g_ctx.frame_zx_zoom) {
		g_ctx.pad_zx_zoom->cd();
		gStyle->SetOptStat(0);
		g_ctx.frame_zx_zoom = new TH2F("frame_zx_zoom", "Z-X Zoom;z (mm);x (mm)",
			200, kZMinZoom, kZMaxZoom, 100, kPosMin, kPosMax);
		g_ctx.frame_zx_zoom->SetDirectory(nullptr);
		g_ctx.frame_zx_zoom->Draw();
		TLine *target_line = new TLine(0, kPosMin, 0, kPosMax);
		target_line->SetLineColor(kGreen + 2);
		target_line->SetLineStyle(2);
		target_line->Draw();
		TText *target_label = new TText(0, kPosMax - 2.0, "target");
		target_label->SetTextSize(0.025);
		target_label->SetTextAlign(22);
		target_label->SetTextColor(kGreen + 2);
		target_label->Draw();
		for (size_t i = 0; i < g_ctx.t0_dets.size(); ++i) {
			auto &det = g_ctx.t0_dets[i];
			if (det.z_mm < kZMinZoom || det.z_mm > kZMaxZoom) continue;
			TLine *t0_line = new TLine(det.z_mm, kPosMin, det.z_mm, kPosMax);
			t0_line->SetLineColor(det.color);
			t0_line->SetLineStyle(3);
			t0_line->SetLineWidth(3);
			t0_line->Draw();
			TText *label = new TText(det.z_mm, kPosMax - 4.0,
				TString::Format("d%d", (int)i + 1));
			label->SetTextSize(0.025);
			label->SetTextAlign(22);
			label->SetTextColor(det.color);
			label->Draw();
		}
	}

	if (!g_ctx.frame_zy_zoom) {
		g_ctx.pad_zy_zoom->cd();
		gStyle->SetOptStat(0);
		g_ctx.frame_zy_zoom = new TH2F("frame_zy_zoom", "Z-Y Zoom;z (mm);y (mm)",
			200, kZMinZoom, kZMaxZoom, 100, kPosMin, kPosMax);
		g_ctx.frame_zy_zoom->SetDirectory(nullptr);
		g_ctx.frame_zy_zoom->Draw();
		TLine *target_line = new TLine(0, kPosMin, 0, kPosMax);
		target_line->SetLineColor(kGreen + 2);
		target_line->SetLineStyle(2);
		target_line->Draw();
		TText *target_label = new TText(0, kPosMax - 2.0, "target");
		target_label->SetTextSize(0.025);
		target_label->SetTextAlign(22);
		target_label->SetTextColor(kGreen + 2);
		target_label->Draw();
		for (size_t i = 0; i < g_ctx.t0_dets.size(); ++i) {
			auto &det = g_ctx.t0_dets[i];
			if (det.z_mm < kZMinZoom || det.z_mm > kZMaxZoom) continue;
			TLine *t0_line = new TLine(det.z_mm, kPosMin, det.z_mm, kPosMax);
			t0_line->SetLineColor(det.color);
			t0_line->SetLineStyle(3);
			t0_line->SetLineWidth(3);
			t0_line->Draw();
			TText *label = new TText(det.z_mm, kPosMax - 4.0,
				TString::Format("d%d", (int)i + 1));
			label->SetTextSize(0.025);
			label->SetTextAlign(22);
			label->SetTextColor(det.color);
			label->Draw();
		}
	}
}

void Draw3DStaticElements() {
	if (!g_ctx.frame_3d && g_ctx.pad_3d) {
		g_ctx.pad_3d->cd();
		gStyle->SetOptStat(0);
		g_ctx.frame_3d = new TH3F("frame_3d", "3D View;z (mm);x (mm);y (mm)",
			40, 80, 150, 40, kPosMin, kPosMax, 40, kPosMin, kPosMax);
		g_ctx.frame_3d->SetDirectory(nullptr);
		g_ctx.frame_3d->Draw();
	}
}

void AddEventToCanvas(const brill::PpacTrackEvent &track) {
	{
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
			g_ctx.zx_graphs.push_back(g_track_zx);
		}
	}

	{
		int n_points = 0;
		double zy_vals_arr[3], zy_points_arr[3];
		for (int i = 0; i < 3; ++i) {
			if (!std::isnan(track.ppac_y[i])) {
				zy_vals_arr[n_points] = g_ctx.ppac_config->z_y_mm[i];
				zy_points_arr[n_points] = track.ppac_y[i];
				++n_points;
			}
		}
		if (n_points > 0) {
			TGraph *g_pt_zy = new TGraph(n_points, zy_vals_arr, zy_points_arr);
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
			g_ctx.zy_graphs.push_back(g_track_zy);
		}
	}

	for (auto &det : g_ctx.t0_dets) {
		if (det.event == nullptr || det.tree == nullptr) continue;
		det.tree->GetEntry(*g_ctx.current_entry);
		for (int i = 0; i < det.event->num; ++i) {
			double x = det.event->x[i];
			double y = det.event->y[i];
			double z = det.z_mm;

			TGraph *g_t0_zx = new TGraph(1);
			g_t0_zx->SetPoint(0, z, x);
			g_t0_zx->SetMarkerStyle(det.marker_style);
			g_t0_zx->SetMarkerColor(det.color);
			g_t0_zx->SetMarkerSize(2.5);
			g_ctx.mg_zx->Add(g_t0_zx);
			g_ctx.zx_graphs.push_back(g_t0_zx);

			TGraph *g_t0_zy = new TGraph(1);
			g_t0_zy->SetPoint(0, z, y);
			g_t0_zy->SetMarkerStyle(det.marker_style);
			g_t0_zy->SetMarkerColor(det.color);
			g_t0_zy->SetMarkerSize(2.5);
			g_ctx.mg_zy->Add(g_t0_zy);
			g_ctx.zy_graphs.push_back(g_t0_zy);
		}
	}

	if (track.valid) {
		{
			TGraph *g_target_zx = new TGraph(1);
			g_target_zx->SetPoint(0, 0, track.target_x);
			g_target_zx->SetMarkerStyle(29);
			g_target_zx->SetMarkerColor(kBlack);
			g_target_zx->SetMarkerSize(2.0);
			g_ctx.mg_zx->Add(g_target_zx);
			g_ctx.zx_graphs.push_back(g_target_zx);
		}
		{
			TGraph *g_target_zy = new TGraph(1);
			g_target_zy->SetPoint(0, 0, track.target_y);
			g_target_zy->SetMarkerStyle(29);
			g_target_zy->SetMarkerColor(kBlack);
			g_target_zy->SetMarkerSize(2.0);
			g_ctx.mg_zy->Add(g_target_zy);
			g_ctx.zy_graphs.push_back(g_target_zy);
		}

		if (g_ctx.t0_dets.size() > 1 && g_ctx.t0_dets[1].event != nullptr &&
			g_ctx.t0_dets[1].tree != nullptr) {
			auto &d2 = g_ctx.t0_dets[1];
			d2.tree->GetEntry(*g_ctx.current_entry);
			for (int i = 0; i < d2.event->num; ++i) {
				double ex_x = track.target_x + (d2.event->x[i] - track.target_x) / d2.z_mm * 150.0;
				double ex_y = track.target_y + (d2.event->y[i] - track.target_y) / d2.z_mm * 150.0;
				{
					TGraph *g_line_zx = new TGraph(2);
					g_line_zx->SetPoint(0, 0, track.target_x);
					g_line_zx->SetPoint(1, 150, ex_x);
					g_line_zx->SetLineColor(i + 2);
					g_line_zx->SetLineWidth(1);
					g_line_zx->SetLineStyle(7);
					g_ctx.mg_zx->Add(g_line_zx);
					g_ctx.zx_graphs.push_back(g_line_zx);
				}
				{
					TGraph *g_line_zy = new TGraph(2);
					g_line_zy->SetPoint(0, 0, track.target_y);
					g_line_zy->SetPoint(1, 150, ex_y);
					g_line_zy->SetLineColor(i + 2);
					g_line_zy->SetLineWidth(1);
					g_line_zy->SetLineStyle(7);
					g_ctx.mg_zy->Add(g_line_zy);
					g_ctx.zy_graphs.push_back(g_line_zy);
				}
			}
		}
	}

	if (g_ctx.mg_zx_zoom) {
		if (track.valid) {
			TGraph *g_track_zx = new TGraph(2);
			g_track_zx->SetPoint(0, kZMinZoom, track.target_x + track.dir_x * kZMinZoom);
			g_track_zx->SetPoint(1, kZMaxZoom, track.target_x + track.dir_x * kZMaxZoom);
			g_track_zx->SetLineColor(kBlue);
			g_track_zx->SetLineWidth(2);
			g_ctx.mg_zx_zoom->Add(g_track_zx);
			g_ctx.zx_zoom_graphs.push_back(g_track_zx);
		}
		for (auto &det : g_ctx.t0_dets) {
			if (det.event == nullptr || det.tree == nullptr) continue;
			if (det.z_mm < kZMinZoom || det.z_mm > kZMaxZoom) continue;
			det.tree->GetEntry(*g_ctx.current_entry);
			for (int i = 0; i < det.event->num; ++i) {
				TGraph *g_t0 = new TGraph(1);
				g_t0->SetPoint(0, det.z_mm, det.event->x[i]);
				g_t0->SetMarkerStyle(det.marker_style);
				g_t0->SetMarkerColor(det.color);
				g_t0->SetMarkerSize(2.5);
				g_ctx.mg_zx_zoom->Add(g_t0);
				g_ctx.zx_zoom_graphs.push_back(g_t0);
			}
		}
		if (track.valid) {
			TGraph *g_target = new TGraph(1);
			g_target->SetPoint(0, 0, track.target_x);
			g_target->SetMarkerStyle(29);
			g_target->SetMarkerColor(kBlack);
			g_target->SetMarkerSize(2.0);
			g_ctx.mg_zx_zoom->Add(g_target);
			g_ctx.zx_zoom_graphs.push_back(g_target);

			if (g_ctx.t0_dets.size() > 1 && g_ctx.t0_dets[1].event != nullptr &&
				g_ctx.t0_dets[1].tree != nullptr) {
				auto &d2 = g_ctx.t0_dets[1];
				d2.tree->GetEntry(*g_ctx.current_entry);
				for (int i = 0; i < d2.event->num; ++i) {
					double ex = track.target_x + (d2.event->x[i] - track.target_x) / d2.z_mm * 150.0;
					TGraph *g_line = new TGraph(2);
					g_line->SetPoint(0, 0, track.target_x);
					g_line->SetPoint(1, 150, ex);
					g_line->SetLineColor(i + 2);
					g_line->SetLineWidth(1);
					g_line->SetLineStyle(7);
					g_ctx.mg_zx_zoom->Add(g_line);
					g_ctx.zx_zoom_graphs.push_back(g_line);
				}
			}
		}
	}

	if (g_ctx.mg_zy_zoom) {
		if (track.valid) {
			TGraph *g_track_zy = new TGraph(2);
			g_track_zy->SetPoint(0, kZMinZoom, track.target_y + track.dir_y * kZMinZoom);
			g_track_zy->SetPoint(1, kZMaxZoom, track.target_y + track.dir_y * kZMaxZoom);
			g_track_zy->SetLineColor(kBlue);
			g_track_zy->SetLineWidth(2);
			g_ctx.mg_zy_zoom->Add(g_track_zy);
			g_ctx.zy_zoom_graphs.push_back(g_track_zy);
		}
		for (auto &det : g_ctx.t0_dets) {
			if (det.event == nullptr || det.tree == nullptr) continue;
			if (det.z_mm < kZMinZoom || det.z_mm > kZMaxZoom) continue;
			det.tree->GetEntry(*g_ctx.current_entry);
			for (int i = 0; i < det.event->num; ++i) {
				TGraph *g_t0 = new TGraph(1);
				g_t0->SetPoint(0, det.z_mm, det.event->y[i]);
				g_t0->SetMarkerStyle(det.marker_style);
				g_t0->SetMarkerColor(det.color);
				g_t0->SetMarkerSize(2.5);
				g_ctx.mg_zy_zoom->Add(g_t0);
				g_ctx.zy_zoom_graphs.push_back(g_t0);
			}
		}
		if (track.valid) {
			TGraph *g_target = new TGraph(1);
			g_target->SetPoint(0, 0, track.target_y);
			g_target->SetMarkerStyle(29);
			g_target->SetMarkerColor(kBlack);
			g_target->SetMarkerSize(2.0);
			g_ctx.mg_zy_zoom->Add(g_target);
			g_ctx.zy_zoom_graphs.push_back(g_target);

			if (g_ctx.t0_dets.size() > 1 && g_ctx.t0_dets[1].event != nullptr &&
				g_ctx.t0_dets[1].tree != nullptr) {
				auto &d2 = g_ctx.t0_dets[1];
				d2.tree->GetEntry(*g_ctx.current_entry);
				for (int i = 0; i < d2.event->num; ++i) {
					double ey = track.target_y + (d2.event->y[i] - track.target_y) / d2.z_mm * 150.0;
					TGraph *g_line = new TGraph(2);
					g_line->SetPoint(0, 0, track.target_y);
					g_line->SetPoint(1, 150, ey);
					g_line->SetLineColor(i + 2);
					g_line->SetLineWidth(1);
					g_line->SetLineStyle(7);
					g_ctx.mg_zy_zoom->Add(g_line);
					g_ctx.zy_zoom_graphs.push_back(g_line);
				}
			}
		}
	}

	if (g_ctx.pad_3d && track.valid) {
		g_ctx.pad_3d->cd();

		TPolyMarker3D *pm_target = new TPolyMarker3D(1);
		pm_target->SetPoint(0, 0, track.target_x, track.target_y);
		pm_target->SetMarkerStyle(29);
		pm_target->SetMarkerColor(kBlack);
		pm_target->SetMarkerSize(2);
		pm_target->Draw();
		g_ctx.objects_3d.push_back(pm_target);

		if (g_ctx.t0_dets.size() > 1 && g_ctx.t0_dets[1].event != nullptr &&
			g_ctx.t0_dets[1].tree != nullptr) {
			auto &d2 = g_ctx.t0_dets[1];
			d2.tree->GetEntry(*g_ctx.current_entry);
			for (int i = 0; i < d2.event->num; ++i) {
				double slope_x = (d2.event->x[i] - track.target_x) / d2.z_mm;
				double slope_y = (d2.event->y[i] - track.target_y) / d2.z_mm;
				double x_at_140 = track.target_x + slope_x * 140;
				double y_at_140 = track.target_y + slope_y * 140;

				TPolyLine3D *pl = new TPolyLine3D(2);
				pl->SetPoint(0, 0, track.target_x, track.target_y);
				pl->SetPoint(1, 140, x_at_140, y_at_140);
				pl->SetLineColor(i + 2);
				pl->SetLineWidth(2);
				pl->Draw();
				g_ctx.objects_3d.push_back(pl);
			}
		}

		for (auto &det : g_ctx.t0_dets) {
			if (det.z_mm < 80 || det.z_mm > 150) continue;
			if (det.tree == nullptr || det.event == nullptr) continue;
			det.tree->GetEntry(*g_ctx.current_entry);
			if (det.event->num > 0) {
				TPolyMarker3D *pm_hits = new TPolyMarker3D(det.event->num);
				for (int i = 0; i < det.event->num; ++i) {
					pm_hits->SetPoint(i, det.z_mm, det.event->x[i], det.event->y[i]);
				}
				pm_hits->SetMarkerStyle(det.marker_style);
				pm_hits->SetMarkerColor(det.color);
				pm_hits->SetMarkerSize(2.5);
				pm_hits->Draw();
				g_ctx.objects_3d.push_back(pm_hits);
			}
		}
	}
}

void UpdateInfoDisplay(const brill::PpacTrackEvent &track, long long entry_num) {
	g_ctx.pad_t0_info->cd();
	g_ctx.t0_info_text->Clear();
	g_ctx.t0_info_text->AddText(TString::Format("Entry: %lld / %lld", entry_num, g_ctx.total));
	g_ctx.t0_info_text->AddText(" ");
	g_ctx.t0_info_text->AddText(TString::Format("Track valid: %d", track.valid));
	g_ctx.t0_info_text->AddText(TString::Format("x_used: %d  y_used: %d",
		track.x_used_ppac, track.y_used_ppac));
	g_ctx.t0_info_text->AddText(TString::Format("target_x: %.3f  target_y: %.3f",
		track.target_x, track.target_y));
	g_ctx.t0_info_text->AddText(TString::Format("dir_x: %.5f  dir_y: %.5f",
		track.dir_x, track.dir_y));
	g_ctx.t0_info_text->AddText(" ");
	g_ctx.t0_info_text->AddText("--- T0 Detectors ---");
	for (auto &det : g_ctx.t0_dets) {
		if (det.event == nullptr || det.tree == nullptr) {
			g_ctx.t0_info_text->AddText(TString::Format("%s: N/A", det.name.c_str()));
			continue;
		}
		det.tree->GetEntry(*g_ctx.current_entry);
		g_ctx.t0_info_text->AddText(TString::Format("%s: num=%d", det.name.c_str(), det.event->num));
		for (int i = 0; i < det.event->num; ++i) {
			g_ctx.t0_info_text->AddText(TString::Format("  hit%d: x=%.3f y=%.3f E=%.3f",
				i + 1, det.event->x[i], det.event->y[i], det.event->energy[i]));
		}
	}
	g_ctx.t0_info_text->Draw();
	g_ctx.canvas2->Modified();
	g_ctx.canvas2->Update();
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
	printf("T0 detectors:\n");
	for (auto &det : g_ctx.t0_dets) {
		if (det.event == nullptr || det.tree == nullptr) {
			printf("  %s: N/A\n", det.name.c_str());
			continue;
		}
		det.tree->GetEntry(*g_ctx.current_entry);
		printf("  %s: num=%d\n", det.name.c_str(), det.event->num);
		for (int i = 0; i < det.event->num; ++i) {
			printf("    hit%d: x=%.3f y=%.3f E=%.3f\n",
				i + 1, det.event->x[i], det.event->y[i], det.event->energy[i]);
		}
	}
}

bool PassFilter() {
	int x_filter = g_ctx.entry_x_used->GetIntNumber();
	int y_filter = g_ctx.entry_y_used->GetIntNumber();
	bool x_ok = (x_filter == 0) || ((int)g_ctx.track_event->x_used_ppac == x_filter);
	bool y_ok = (y_filter == 0) || ((int)g_ctx.track_event->y_used_ppac == y_filter);
	if (!x_ok || !y_ok) return false;

	int d1_filter = g_ctx.entry_d1_num->GetIntNumber();
	int d2_filter = g_ctx.entry_d2_num->GetIntNumber();
	int d3_filter = g_ctx.entry_d3_num->GetIntNumber();
	int d4_filter = g_ctx.entry_d4_num->GetIntNumber();

	if (d1_filter != -1) {
		if (g_ctx.t0_dets.size() < 1 || g_ctx.t0_dets[0].event == nullptr ||
			g_ctx.t0_dets[0].tree == nullptr) return false;
		g_ctx.t0_dets[0].tree->GetEntry(*g_ctx.current_entry);
		if (g_ctx.t0_dets[0].event->num != d1_filter) return false;
	}
	if (d2_filter != -1) {
		if (g_ctx.t0_dets.size() < 2 || g_ctx.t0_dets[1].event == nullptr ||
			g_ctx.t0_dets[1].tree == nullptr) return false;
		g_ctx.t0_dets[1].tree->GetEntry(*g_ctx.current_entry);
		if (g_ctx.t0_dets[1].event->num != d2_filter) return false;
	}
	if (d3_filter != -1) {
		if (g_ctx.t0_dets.size() < 3 || g_ctx.t0_dets[2].event == nullptr ||
			g_ctx.t0_dets[2].tree == nullptr) return false;
		g_ctx.t0_dets[2].tree->GetEntry(*g_ctx.current_entry);
		if (g_ctx.t0_dets[2].event->num != d3_filter) return false;
	}
	if (d4_filter != -1) {
		if (g_ctx.t0_dets.size() < 4 || g_ctx.t0_dets[3].event == nullptr ||
			g_ctx.t0_dets[3].tree == nullptr) return false;
		g_ctx.t0_dets[3].tree->GetEntry(*g_ctx.current_entry);
		if (g_ctx.t0_dets[3].event->num != d4_filter) return false;
	}

	return true;
}

long long FindNextMatching(long long start) {
	while (start < g_ctx.total) {
		g_ctx.input_tree->GetEntry(start);
		*g_ctx.current_entry = start;
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
	g_ctx.zx_graphs.clear();
	g_ctx.zy_graphs.clear();

	if (g_ctx.mg_zx_zoom) {
		for (auto *g : g_ctx.zx_zoom_graphs) {
			TList *list = g_ctx.mg_zx_zoom->GetListOfGraphs();
			if (list) list->Remove(g);
			delete g;
		}
		g_ctx.zx_zoom_graphs.clear();
	}
	if (g_ctx.mg_zy_zoom) {
		for (auto *g : g_ctx.zy_zoom_graphs) {
			TList *list = g_ctx.mg_zy_zoom->GetListOfGraphs();
			if (list) list->Remove(g);
			delete g;
		}
		g_ctx.zy_zoom_graphs.clear();
	}

	g_ctx.objects_3d.clear();
	if (g_ctx.pad_3d) {
		g_ctx.pad_3d->cd();
		g_ctx.pad_3d->Clear();
		g_ctx.frame_3d = nullptr;
		Draw3DStaticElements();
	}
}

void DrawZoomCanvas() {
	if (!g_ctx.canvas3) return;
	if (g_ctx.mg_zx_zoom) {
		g_ctx.pad_zx_zoom->cd();
		TList *list = g_ctx.mg_zx_zoom->GetListOfGraphs();
		if (list && list->GetSize() > 0) {
			for (int i = 0; i < list->GetSize(); ++i) {
				TGraph *g = (TGraph*)list->At(i);
				if (g) g->Draw("LP same");
			}
		}
	}
	if (g_ctx.mg_zy_zoom) {
		g_ctx.pad_zy_zoom->cd();
		TList *list = g_ctx.mg_zy_zoom->GetListOfGraphs();
		if (list && list->GetSize() > 0) {
			for (int i = 0; i < list->GetSize(); ++i) {
				TGraph *g = (TGraph*)list->At(i);
				if (g) g->Draw("LP same");
			}
		}
	}
	g_ctx.canvas3->Modified();
	g_ctx.canvas3->Update();
}

void Draw3DCanvas() {
	if (!g_ctx.canvas4) return;
	g_ctx.canvas4->Modified();
	g_ctx.canvas4->Update();
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
	*g_ctx.current_entry = found;
	g_ctx.entry_number->SetIntNumber(found);
	AddEventToCanvas(*g_ctx.track_event);
	UpdateInfoDisplay(*g_ctx.track_event, found);
	PrintEvent(*g_ctx.track_event, found);

	g_ctx.pad_zx->cd();
	{
		TList *list = g_ctx.mg_zx->GetListOfGraphs();
		if (list && list->GetSize() > 0) {
			for (int i = 0; i < list->GetSize(); ++i) {
				TGraph *g = (TGraph*)list->At(i);
				if (g) g->Draw("LP same");
			}
		}
	}
	g_ctx.pad_zy->cd();
	{
		TList *list = g_ctx.mg_zy->GetListOfGraphs();
		if (list && list->GetSize() > 0) {
			for (int i = 0; i < list->GetSize(); ++i) {
				TGraph *g = (TGraph*)list->At(i);
				if (g) g->Draw("LP same");
			}
		}
	}

	g_ctx.canvas1->Modified();
	g_ctx.canvas1->Update();
	DrawZoomCanvas();
	Draw3DCanvas();
}

void ProcessNext() {
	long long start = *g_ctx.current_entry + 1;
	long long found = FindNextMatching(start);
	if (found < 0) {
		printf("\nNo more matching events.\n");
		return;
	}

	ClearCanvas();
	*g_ctx.current_entry = found;
	g_ctx.entry_number->SetIntNumber(found);
	AddEventToCanvas(*g_ctx.track_event);
	UpdateInfoDisplay(*g_ctx.track_event, found);
	PrintEvent(*g_ctx.track_event, found);

	g_ctx.pad_zx->cd();
	{
		TList *list = g_ctx.mg_zx->GetListOfGraphs();
		if (list && list->GetSize() > 0) {
			for (int i = 0; i < list->GetSize(); ++i) {
				TGraph *g = (TGraph*)list->At(i);
				if (g) g->Draw("LP same");
			}
		}
	}
	g_ctx.pad_zy->cd();
	{
		TList *list = g_ctx.mg_zy->GetListOfGraphs();
		if (list && list->GetSize() > 0) {
			for (int i = 0; i < list->GetSize(); ++i) {
				TGraph *g = (TGraph*)list->At(i);
				if (g) g->Draw("LP same");
			}
		}
	}

	g_ctx.canvas1->Modified();
	g_ctx.canvas1->Update();
	DrawZoomCanvas();
	Draw3DCanvas();
}

void ProcessAdd() {
	long long start = *g_ctx.current_entry + 1;
	long long found = FindNextMatching(start);
	if (found < 0) {
		printf("\nNo more matching events.\n");
		return;
	}

	*g_ctx.current_entry = found;
	g_ctx.entry_number->SetIntNumber(found);
	AddEventToCanvas(*g_ctx.track_event);
	UpdateInfoDisplay(*g_ctx.track_event, found);
	PrintEvent(*g_ctx.track_event, found);

	g_ctx.pad_zx->cd();
	{
		TList *list = g_ctx.mg_zx->GetListOfGraphs();
		if (list && list->GetSize() > 0) {
			for (int i = 0; i < list->GetSize(); ++i) {
				TGraph *g = (TGraph*)list->At(i);
				if (g) g->Draw("LP same");
			}
		}
	}
	g_ctx.pad_zy->cd();
	{
		TList *list = g_ctx.mg_zy->GetListOfGraphs();
		if (list && list->GetSize() > 0) {
			for (int i = 0; i < list->GetSize(); ++i) {
				TGraph *g = (TGraph*)list->At(i);
				if (g) g->Draw("LP same");
			}
		}
	}

	g_ctx.canvas1->Modified();
	g_ctx.canvas1->Update();
	DrawZoomCanvas();
	Draw3DCanvas();
}

int main(int argc, char **argv) {
	cxxopts::Options options("GUI_track", "Interactive PPAC+T0 track checker.");
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
	std::string match_dir = brill::JoinPath(config.workspace, config.paths.match);

	TString track_path = TString::Format("%s/ppac_%s%04d.root",
		track_dir.c_str(), brill::TriggerInfix(config.trigger).c_str(), run);

	std::cout << "PPAC track:  " << track_path << "\n";
	std::cout << "Match dir:   " << match_dir << "\n";

	if (!std::filesystem::exists(track_path.Data())) {
		std::cerr << "Error: PPAC track file not found: " << track_path << "\n";
		return 1;
	}

	TFile *input_file = new TFile(track_path, "read");
	TTree *input_tree = (TTree*)input_file->Get("tree");
	if (!input_tree) { std::cerr << "Error: Get tree failed.\n"; return 1; }

	brill::PpacTrackEvent track_event;
	brill::SetupInput(input_tree, track_event, "");

	long long total = input_tree->GetEntries();
	printf("PPAC total events: %lld\n", total);

	double t0_z_min = 0.0, t0_z_max = 500.0;
	{
		bool first = true;
		for (const char *name : {"t0d1", "t0d2", "t0d3", "t0d4"}) {
			const auto *det = brill::FindDetectorConfig(config, name);
			if (det) {
				if (first) {
					t0_z_min = det->z_mm;
					t0_z_max = det->z_mm;
					first = false;
				} else {
					if (det->z_mm < t0_z_min) t0_z_min = det->z_mm;
					if (det->z_mm > t0_z_max) t0_z_max = det->z_mm;
				}
			}
		}
		if (!first) {
			t0_z_min -= 50.0;
			t0_z_max += 50.0;
		}
	}
	g_ctx.t0_z_min = t0_z_min;
	g_ctx.t0_z_max = t0_z_max;

	TApplication app("GUI_track", nullptr, nullptr);

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
	TString title = TString::Format("GUI_track - Run %d, Trigger %s",
		run, config.trigger.c_str());
	main_frame->SetWindowName(title.Data());

	TRootEmbeddedCanvas *embed1 = new TRootEmbeddedCanvas("embed1", main_frame, 1400, 700);
	main_frame->AddFrame(embed1, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));
	TCanvas *canvas1 = embed1->GetCanvas();

	TPad *pad_zx = new TPad("pad_zx", "Z-X Track", 0.00, 0.50, 1.00, 1.00);
	TPad *pad_zy = new TPad("pad_zy", "Z-Y Track", 0.00, 0.00, 1.00, 0.50);

	canvas1->cd();
	pad_zx->Draw();
	pad_zy->Draw();

	pad_zx->SetGrid();
	pad_zy->SetGrid();

	TMultiGraph *mg_zx = new TMultiGraph();
	TMultiGraph *mg_zy = new TMultiGraph();

	TCanvas *canvas2 = new TCanvas("canvas2", "T0 Info & Fitting", 800, 600);
	canvas2->SetWindowPosition(1420, 0);

	TPad *pad_t0_info = new TPad("pad_t0_info", "T0 Info", 0.00, 0.00, 0.50, 1.00);
	TPad *pad_fitting = new TPad("pad_fitting", "Fitting", 0.50, 0.00, 1.00, 1.00);

	canvas2->cd();
	pad_t0_info->Draw();
	pad_fitting->Draw();

	pad_fitting->cd();
	TPaveText *fitting_placeholder = new TPaveText(0.1, 0.4, 0.9, 0.6, "NB");
	fitting_placeholder->AddText("(fitting - reserved)");
	fitting_placeholder->SetTextSize(0.06);
	fitting_placeholder->SetFillColor(kWhite);
	fitting_placeholder->SetBorderSize(1);
	fitting_placeholder->Draw();

	TPaveText *t0_info_text = new TPaveText(0.05, 0.05, 0.95, 0.95, "NB");
	t0_info_text->SetTextSize(0.03);
	t0_info_text->SetFillColor(kWhite);
	t0_info_text->SetBorderSize(1);

	long long current_entry = 0;

	TCanvas *canvas3 = new TCanvas("canvas3", "Zoom View (80-150)", 800, 400);
	canvas3->SetWindowPosition(1420, 620);

	TPad *pad_zx_zoom = new TPad("pad_zx_zoom", "Z-X Zoom", 0.00, 0.00, 0.50, 1.00);
	TPad *pad_zy_zoom = new TPad("pad_zy_zoom", "Z-Y Zoom", 0.50, 0.00, 1.00, 1.00);

	canvas3->cd();
	pad_zx_zoom->Draw();
	pad_zy_zoom->Draw();

	pad_zx_zoom->SetGrid();
	pad_zy_zoom->SetGrid();

	TMultiGraph *mg_zx_zoom = new TMultiGraph();
	TMultiGraph *mg_zy_zoom = new TMultiGraph();

	TCanvas *canvas4 = new TCanvas("canvas4", "3D View", 500, 500);
	canvas4->SetWindowPosition(2230, 0);

	TPad *pad_3d = new TPad("pad_3d", "3D", 0.00, 0.00, 1.00, 1.00);
	canvas4->cd();
	pad_3d->Draw();

	g_ctx.canvas1 = canvas1;
	g_ctx.canvas2 = canvas2;
	g_ctx.canvas3 = canvas3;
	g_ctx.canvas4 = canvas4;
	g_ctx.pad_zx = pad_zx;
	g_ctx.pad_zy = pad_zy;
	g_ctx.pad_zx_zoom = pad_zx_zoom;
	g_ctx.pad_zy_zoom = pad_zy_zoom;
	g_ctx.pad_3d = pad_3d;
	g_ctx.pad_t0_info = pad_t0_info;
	g_ctx.pad_fitting = pad_fitting;
	g_ctx.mg_zx = mg_zx;
	g_ctx.mg_zy = mg_zy;
	g_ctx.mg_zx_zoom = mg_zx_zoom;
	g_ctx.mg_zy_zoom = mg_zy_zoom;
	g_ctx.t0_info_text = t0_info_text;
	g_ctx.input_tree = input_tree;
	g_ctx.track_event = &track_event;
	g_ctx.ppac_config = &config.ppac;
	g_ctx.current_entry = &current_entry;
	g_ctx.total = total;

	T0DetectorEntry t0d1;
	t0d1.name = "t0d1";
	{
		const auto *det = brill::FindDetectorConfig(config, "t0d1");
		if (det) {
			t0d1.z_mm = det->z_mm;
		} else {
			std::cerr << "Warning: t0d1 config not found, using z=0\n";
			t0d1.z_mm = 0.0;
		}
	}
	t0d1.color = kGreen + 2;
	t0d1.marker_style = 21;
	TString t0d1_path = TString::Format("%s/t0d1_%s%04d.root",
		match_dir.c_str(), brill::TriggerInfix(config.trigger).c_str(), run);
	if (std::filesystem::exists(t0d1_path.Data())) {
		t0d1.file = new TFile(t0d1_path, "read");
		t0d1.tree = (TTree*)t0d1.file->Get("tree");
		if (t0d1.tree) {
			t0d1.event = new brill::DssdMatchEvent();
			brill::SetupInput(t0d1.tree, *t0d1.event, "");
			printf("T0 loaded: %s -> %s (%lld entries)\n",
				t0d1.name.c_str(), t0d1_path.Data(), t0d1.tree->GetEntries());
		} else {
			std::cerr << "Warning: t0d1 tree not found in " << t0d1_path << "\n";
		}
	} else {
		std::cerr << "Warning: t0d1 file not found: " << t0d1_path << "\n";
	}
	g_ctx.t0_dets.push_back(t0d1);

	T0DetectorEntry t0d2;
	t0d2.name = "t0d2";
	{
		const auto *det = brill::FindDetectorConfig(config, "t0d2");
		if (det) {
			t0d2.z_mm = det->z_mm;
		} else {
			std::cerr << "Warning: t0d2 config not found, using z=0\n";
			t0d2.z_mm = 0.0;
		}
	}
	t0d2.color = kRed;
	t0d2.marker_style = 22;
	TString t0d2_path = TString::Format("%s/t0d2_%s%04d.root",
		match_dir.c_str(), brill::TriggerInfix(config.trigger).c_str(), run);
	if (std::filesystem::exists(t0d2_path.Data())) {
		t0d2.file = new TFile(t0d2_path, "read");
		t0d2.tree = (TTree*)t0d2.file->Get("tree");
		if (t0d2.tree) {
			t0d2.event = new brill::DssdMatchEvent();
			brill::SetupInput(t0d2.tree, *t0d2.event, "");
			printf("T0 loaded: %s -> %s (%lld entries)\n",
				t0d2.name.c_str(), t0d2_path.Data(), t0d2.tree->GetEntries());
		} else {
			std::cerr << "Warning: t0d2 tree not found in " << t0d2_path << "\n";
		}
	} else {
		std::cerr << "Warning: t0d2 file not found: " << t0d2_path << "\n";
	}
	g_ctx.t0_dets.push_back(t0d2);

	T0DetectorEntry t0d3;
	t0d3.name = "t0d3";
	{
		const auto *det = brill::FindDetectorConfig(config, "t0d3");
		if (det) {
			t0d3.z_mm = det->z_mm;
		} else {
			std::cerr << "Warning: t0d3 config not found, using z=0\n";
			t0d3.z_mm = 0.0;
		}
	}
	t0d3.color = kBlue;
	t0d3.marker_style = 23;
	TString t0d3_path = TString::Format("%s/t0d3_%s%04d.root",
		match_dir.c_str(), brill::TriggerInfix(config.trigger).c_str(), run);
	if (std::filesystem::exists(t0d3_path.Data())) {
		t0d3.file = new TFile(t0d3_path, "read");
		t0d3.tree = (TTree*)t0d3.file->Get("tree");
		if (t0d3.tree) {
			t0d3.event = new brill::DssdMatchEvent();
			brill::SetupInput(t0d3.tree, *t0d3.event, "");
			printf("T0 loaded: %s -> %s (%lld entries)\n",
				t0d3.name.c_str(), t0d3_path.Data(), t0d3.tree->GetEntries());
		} else {
			std::cerr << "Warning: t0d3 tree not found in " << t0d3_path << "\n";
		}
	} else {
		std::cerr << "Warning: t0d3 file not found: " << t0d3_path << "\n";
	}
	g_ctx.t0_dets.push_back(t0d3);

	T0DetectorEntry t0d4;
	t0d4.name = "t0d4";
	{
		const auto *det = brill::FindDetectorConfig(config, "t0d4");
		if (det) {
			t0d4.z_mm = det->z_mm;
		} else {
			std::cerr << "Warning: t0d4 config not found, using z=0\n";
			t0d4.z_mm = 0.0;
		}
	}
	t0d4.color = kMagenta;
	t0d4.marker_style = 24;
	TString t0d4_path = TString::Format("%s/t0d4_%s%04d.root",
		match_dir.c_str(), brill::TriggerInfix(config.trigger).c_str(), run);
	if (std::filesystem::exists(t0d4_path.Data())) {
		t0d4.file = new TFile(t0d4_path, "read");
		t0d4.tree = (TTree*)t0d4.file->Get("tree");
		if (t0d4.tree) {
			t0d4.event = new brill::DssdMatchEvent();
			brill::SetupInput(t0d4.tree, *t0d4.event, "");
			printf("T0 loaded: %s -> %s (%lld entries)\n",
				t0d4.name.c_str(), t0d4_path.Data(), t0d4.tree->GetEntries());
		} else {
			std::cerr << "Warning: t0d4 tree not found in " << t0d4_path << "\n";
		}
	} else {
		std::cerr << "Warning: t0d4 file not found: " << t0d4_path << "\n";
	}
	g_ctx.t0_dets.push_back(t0d4);

	DrawStaticElements();
	DrawStaticElementsZoom();
	Draw3DStaticElements();

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

	ctrl_frame->AddFrame(new TGLabel(ctrl_frame, "  d1_num:"),
		new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 8, 2, 5, 5));
	TGNumberEntry *entry_d1_num = new TGNumberEntry(ctrl_frame, -1, 3, -1,
		TGNumberFormat::kNESInteger,
		TGNumberFormat::kNEAAnyNumber,
		TGNumberFormat::kNELLimitMinMax, -1, 99);
	ctrl_frame->AddFrame(entry_d1_num,
		new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 5, 5, 5));
	g_ctx.entry_d1_num = entry_d1_num;

	ctrl_frame->AddFrame(new TGLabel(ctrl_frame, "  d2_num:"),
		new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 8, 2, 5, 5));
	TGNumberEntry *entry_d2_num = new TGNumberEntry(ctrl_frame, -1, 3, -1,
		TGNumberFormat::kNESInteger,
		TGNumberFormat::kNEAAnyNumber,
		TGNumberFormat::kNELLimitMinMax, -1, 99);
	ctrl_frame->AddFrame(entry_d2_num,
		new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 5, 5, 5));
	g_ctx.entry_d2_num = entry_d2_num;

	ctrl_frame->AddFrame(new TGLabel(ctrl_frame, "  d3_num:"),
		new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 8, 2, 5, 5));
	TGNumberEntry *entry_d3_num = new TGNumberEntry(ctrl_frame, -1, 3, -1,
		TGNumberFormat::kNESInteger,
		TGNumberFormat::kNEAAnyNumber,
		TGNumberFormat::kNELLimitMinMax, -1, 99);
	ctrl_frame->AddFrame(entry_d3_num,
		new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 5, 5, 5));
	g_ctx.entry_d3_num = entry_d3_num;

	ctrl_frame->AddFrame(new TGLabel(ctrl_frame, "  d4_num:"),
		new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 8, 2, 5, 5));
	TGNumberEntry *entry_d4_num = new TGNumberEntry(ctrl_frame, -1, 3, -1,
		TGNumberFormat::kNESInteger,
		TGNumberFormat::kNEAAnyNumber,
		TGNumberFormat::kNELLimitMinMax, -1, 99);
	ctrl_frame->AddFrame(entry_d4_num,
		new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 5, 5, 5));
	g_ctx.entry_d4_num = entry_d4_num;

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
	printf("=== GUI_track ===\n");
	printf("Canvas 1 - Track views:\n");
	printf("  Top:    Z-X track (PPAC + T0)\n");
	printf("  Bottom: Z-Y track (PPAC + T0)\n");
	printf("Canvas 2 - Info:\n");
	printf("  Left:  T0 detector raw info\n");
	printf("  Right: Fitting (reserved)\n");
	printf("Canvas 3 - Zoom (80-150):\n");
	printf("  Left:  Z-X track zoom\n");
	printf("  Right: Z-Y track zoom\n");
	printf("Canvas 4 - 3D View:\n");
	printf("  target -> t0d1 hit lines in 3D\n");
	printf("Controls:\n");
	printf("  Entry:       set entry number\n");
	printf("  x_used_ppac: filter by x_used_ppac (0=all, 1-7=specific)\n");
	printf("  y_used_ppac: filter by y_used_ppac (0=all, 1-7=specific)\n");
	printf("  d1-d4_num:   T0 hit multiplicity (-1=any, 0=strictly 0, N=exactly N)\n");
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