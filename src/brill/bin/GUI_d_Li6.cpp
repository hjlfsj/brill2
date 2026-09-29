#include "include/config.h"
#include "include/d_6Li/d_6Li_event.h"
#include "include/Lise++/e_theta.h"
#include "include/Lise++/target_energy_loss.h"
#include "include/physics/kinematics.h"
#include "include/rebuild/nuclear_data.h"
#include "include/rebuild/rebuild_d_6Li.h"
#include "include/utils.h"
#include "external/cxxopts.hpp"

#include <TApplication.h>
#include <TCanvas.h>
#include <TCutG.h>
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
#include <TGraph.h>
#include <TH2D.h>
#include <TPad.h>
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

struct ExcitationCanvas {
	TCanvas *canvas = nullptr;
	TH1D *h_excitation = nullptr;
	TH1D *h_excitation_10C_mm = nullptr;
	TH1D *h_excitation_6Li_mm = nullptr;
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
	TGNumberEntry *entry_excitation_bins = nullptr;
	TGNumberEntry *entry_total_E = nullptr;
	TGCheckButton *chk_energy_loss = nullptr;

	BeamCanvases bcs;
	AnalysisCanvas ac;
	ExcitationCanvas ec;

	std::string current_file;
	std::string config_path;
	std::string d_Li6_dir;
	std::string assets_dir;
	brill::AppConfig config;
	std::vector<brill::D6LiEvent> all_events;
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

static bool PassBeamFilter(const brill::D6LiEvent &ev) {
	switch (g_beam_filter) {
		case 1: return ev.is_14O;
		case 2: return ev.is_13N;
		case 3: return ev.is_12C;
		default: return true;
	}
}

static void RebuildHistograms() {
	auto &bc = g_ctx.bcs;
	if (bc.h_e1_10C_e2_10C) delete bc.h_e1_10C_e2_10C;
	if (bc.h_e1_6Li_e2_6Li) delete bc.h_e1_6Li_e2_6Li;
	if (bc.h_e2_10C_e3_10C) delete bc.h_e2_10C_e3_10C;
	if (bc.h_e3_10C_e4_10C) delete bc.h_e3_10C_e4_10C;

	TString label = TString::Format(" (%s)", BeamLabel());

	bc.h_e1_10C_e2_10C = new TH2D(
		"h_e1_10C_e2_10C",
		TString::Format("E1_10C vs E2_10C%s;E2_10C (MeV);E1_10C (MeV)", label.Data()),
		1000, 0, 400, 1000, 0, 25);
	bc.h_e1_10C_e2_10C->SetDirectory(0);
	bc.h_e1_6Li_e2_6Li = new TH2D(
		"h_e1_6Li_e2_6Li",
		TString::Format("E1_6Li vs E2_6Li%s;E2_6Li (MeV);E1_6Li (MeV)", label.Data()),
		1000, 0, 400, 1000, 0, 25);
	bc.h_e1_6Li_e2_6Li->SetDirectory(0);
	bc.h_e2_10C_e3_10C = new TH2D(
		"h_e2_10C_e3_10C",
		TString::Format("E2_10C vs E3_10C%s;E3_10C (MeV);E2_10C (MeV)", label.Data()),
		1000, 0, 350, 1000, 0, 400);
	bc.h_e2_10C_e3_10C->SetDirectory(0);
	bc.h_e3_10C_e4_10C = new TH2D(
		"h_e3_10C_e4_10C",
		TString::Format("E3_10C vs E4_10C%s;E4_10C (MeV);E3_10C (MeV)", label.Data()),
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
		120, 350, 470, 50, 30, 80);
	ac.h_E_6Li_E_10C->SetDirectory(0);
	ac.h_10C_e_theta = new TH2D(
		"h_10C_e_theta",
		"^{10}C E-#theta;#theta_{10C} (deg);E_{10C} (MeV)",
		20, 0, 20, 170, 300, 470);
	ac.h_10C_e_theta->SetDirectory(0);
	ac.h_6Li_e_theta = new TH2D(
		"h_6Li_e_theta",
		"^{6}Li E-#theta;#theta_{6Li} (deg);E_{6Li} (MeV)",
		25, 0, 25, 100, 0, 100);
	ac.h_6Li_e_theta->SetDirectory(0);
	ac.h_theta_theta = new TH2D(
		"h_theta_theta",
		"#theta_{10C} vs #theta_{6Li};#theta_{10C} (deg);#theta_{6Li} (deg)",
		20, 0, 20, 30, 0, 30);
	ac.h_theta_theta->SetDirectory(0);
}

static void RebuildExcitationHistograms() {
	auto &ec = g_ctx.ec;
	if (ec.h_excitation) delete ec.h_excitation;
	if (ec.h_excitation_10C_mm) delete ec.h_excitation_10C_mm;
	if (ec.h_excitation_6Li_mm) delete ec.h_excitation_6Li_mm;
	int nbins = g_ctx.entry_excitation_bins
		? g_ctx.entry_excitation_bins->GetIntNumber() : 30;
	if (nbins < 5) nbins = 5;
	double bin_keV = 30000.0 / nbins;
	TString bin_label = TString::Format("Counts / %.0f keV", bin_keV);

	ec.h_excitation = new TH1D(
		"h_excitation",
		TString::Format("^{10}C Excitation (p-consv, %s);E_{x} (MeV);%s",
			BeamLabel(), bin_label.Data()),
		nbins, -10, 20);
	ec.h_excitation->SetDirectory(0);

	ec.h_excitation_10C_mm = new TH1D(
		"h_excitation_10C_mm",
		TString::Format("^{10}C Excitation (^{10}C miss, %s);E_{x} (MeV);%s",
			BeamLabel(), bin_label.Data()),
		nbins, -10, 20);
	ec.h_excitation_10C_mm->SetDirectory(0);

	ec.h_excitation_6Li_mm = new TH1D(
		"h_excitation_6Li_mm",
		TString::Format("^{10}C Excitation (^{6}Li miss, %s);E_{x} (MeV);%s",
			BeamLabel(), bin_label.Data()),
		nbins, -10, 20);
	ec.h_excitation_6Li_mm->SetDirectory(0);
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

		if (!PassBeamFilter(ev)) continue;

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
	RebuildExcitationHistograms();

	if (g_ctx.all_events.empty()) return;

	const double m_10C = brill::GetMass(6, 10);
	const double m_6Li = brill::GetMass(3, 6);
	const double m_14O = brill::GetMass(8, 14);
	const double m_d = brill::GetMass(1, 2);
	const double Q = m_14O + m_d - m_10C - m_6Li;

	brill::TargetEnergyLoss tel_10C(g_ctx.config, 6, 10);
	brill::TargetEnergyLoss tel_6Li(g_ctx.config, 3, 6);

	int total_events = (int)g_ctx.all_events.size();
	int passed = 0;

	printf("  --- Selected events (run, entry) ---\n");
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

		if (!PassBeamFilter(ev)) continue;
		if (!ev.ppac_valid) continue;
		if (!cut->IsInside(ev.e2_6Li, ev.e1_6Li)) continue;

		printf("    run=%d  entry=%lld\n", ev.run_number, ev.entry);
		passed++;

		double E_10C_raw = ev.e1_10C + ev.e2_10C + ev.e3_10C + ev.e4_10C;
		double E_6Li_raw = ev.e1_6Li + ev.e2_6Li;

		bool apply_loss = g_ctx.chk_energy_loss
			? g_ctx.chk_energy_loss->IsOn() : true;
		double E_10C_corr, E_6Li_corr;
		if (apply_loss) {
			double dummy;
			E_10C_corr = tel_10C.IncidentEnergy(E_10C_raw, ev.theta_10C,
				brill::TargetType::kCD2, dummy);
			E_6Li_corr = tel_6Li.IncidentEnergy(E_6Li_raw, ev.theta_6Li,
				brill::TargetType::kCD2, dummy);
		} else {
			E_10C_corr = E_10C_raw;
			E_6Li_corr = E_6Li_raw;
		}

		g_ctx.ac.h_E_6Li_E_10C->Fill(E_10C_corr, E_6Li_corr);
		g_ctx.ac.h_10C_e_theta->Fill(ev.theta_10C, E_10C_corr);
		g_ctx.ac.h_6Li_e_theta->Fill(ev.theta_6Li, E_6Li_corr);
		g_ctx.ac.h_theta_theta->Fill(ev.theta_10C, ev.theta_6Li);

		double p_10C = std::sqrt(2.0 * m_10C * E_10C_corr);
		double p_6Li = std::sqrt(2.0 * m_6Li * E_6Li_corr);

		double sin_10C = std::sin(ev.theta_10C * M_PI / 180.0);
		double cos_10C = std::cos(ev.theta_10C * M_PI / 180.0);
		double sin_6Li = std::sin(ev.theta_6Li * M_PI / 180.0);
		double cos_6Li = std::cos(ev.theta_6Li * M_PI / 180.0);

		double pT_10C = p_10C * sin_10C;
		double pT_6Li = p_6Li * sin_6Li;
		// coplanarity: φ_6Li = φ_10C + π → transverse momenta opposite
		double pT_diff = pT_10C - pT_6Li;
		double pz_sum  = p_10C * cos_10C + p_6Li * cos_6Li;
		double p2_14O = pT_diff * pT_diff + pz_sum * pz_sum;

		double E_14O = p2_14O / (2.0 * m_14O);
		double E_x = E_14O - E_10C_corr - E_6Li_corr + Q;

		double total_E_cut = g_ctx.entry_total_E
			? g_ctx.entry_total_E->GetNumber() : 0.0;
		if (E_10C_corr + E_6Li_corr > total_E_cut) {
			g_ctx.ec.h_excitation->Fill(E_x);

			// use 10C position for φ reference, then stored θ + coplanarity
			double rho_10C = std::sqrt(
				(ev.t0d2_10C_x - ev.target_x) * (ev.t0d2_10C_x - ev.target_x) +
				(ev.t0d2_10C_y - ev.target_y) * (ev.t0d2_10C_y - ev.target_y));
			double cos_phi = rho_10C > 0
				? (ev.t0d2_10C_x - ev.target_x) / rho_10C : 1.0;
			double sin_phi = rho_10C > 0
				? (ev.t0d2_10C_y - ev.target_y) / rho_10C : 0.0;

			// relativistic momenta from stored θ + coplanarity
			double E_10C_tot_kin = m_10C + E_10C_corr;
			double p_rel_10C = std::sqrt(E_10C_tot_kin * E_10C_tot_kin - m_10C * m_10C);
			double px_10C = p_rel_10C * sin_10C * cos_phi;
			double py_10C = p_rel_10C * sin_10C * sin_phi;
			double pz_10C = p_rel_10C * cos_10C;

			double E_6Li_tot_kin = m_6Li + E_6Li_corr;
			double p_rel_6Li = std::sqrt(E_6Li_tot_kin * E_6Li_tot_kin - m_6Li * m_6Li);
			double px_6Li = p_rel_6Li * sin_6Li * cos_phi;
			double py_6Li = p_rel_6Li * sin_6Li * sin_phi;
			double pz_6Li = p_rel_6Li * cos_6Li;

			double E_14O_tot = m_14O + 35.3 * 14.0;
			double p_14O_sq = E_14O_tot * E_14O_tot - m_14O * m_14O;
			double p_beam_mag = std::sqrt(p_14O_sq);

			double bx = ev.dir_x;
			double by = ev.dir_y;
			double bz = 1.0;
			double br = std::sqrt(bx*bx + by*by + bz*bz);
			double pbx = p_beam_mag * bx / br;
			double pby = p_beam_mag * by / br;
			double pbz = p_beam_mag * bz / br;

			double E_6Li_tot = m_6Li + E_6Li_corr;
			double dE_10C_mm = E_14O_tot + m_d - E_6Li_tot;
			double dpx_10C = pbx - px_6Li;
			double dpy_10C = pby - py_6Li;
			double dpz_10C = pbz - pz_6Li;
			double dp2_10C = dpx_10C*dpx_10C + dpy_10C*dpy_10C + dpz_10C*dpz_10C;
			double M_10C_star = std::sqrt(dE_10C_mm * dE_10C_mm - dp2_10C);
			double E_x_10C_mm = M_10C_star - m_10C;

			double p10_sq_m = px_10C*px_10C + py_10C*py_10C + pz_10C*pz_10C;
			double dpx_6Li = pbx - px_10C;
			double dpy_6Li = pby - py_10C;
			double dpz_6Li = pbz - pz_10C;
			double dp2_6Li = dpx_6Li*dpx_6Li + dpy_6Li*dpy_6Li + dpz_6Li*dpz_6Li;
			double E_6Li_recon = std::sqrt(m_6Li * m_6Li + dp2_6Li);
			double dE_6Li_mm = E_14O_tot + m_d - E_6Li_recon;
			double M_10C_star_6Li = std::sqrt(dE_6Li_mm * dE_6Li_mm - p10_sq_m);
			double E_x_6Li_mm = M_10C_star_6Li - m_10C;

			g_ctx.ec.h_excitation_10C_mm->Fill(E_x_10C_mm);
			g_ctx.ec.h_excitation_6Li_mm->Fill(E_x_6Li_mm);
		}
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
	ref_10C->SetLineColor(kGreen + 2);
	ref_10C->SetLineWidth(2);
	TGraph *ref_6Li = brill::LoadEThetaCurve(
		brill::JoinPath(g_ctx.assets_dir, "14O_d_6Li_0+_e_theta.txt"), "6Li");
	ref_6Li->SetLineColor(kGreen + 2);
	ref_6Li->SetLineWidth(2);
	TGraph *ref_theta = brill::LoadThetaThetaCurve(
		brill::JoinPath(g_ctx.assets_dir, "14O_d_6Li_0+_theta_theta.txt"));
	ref_theta->SetLineColor(kGreen + 2);
	ref_theta->SetLineWidth(2);

	TGraph *ref_10C_2p = brill::LoadEThetaCurve(
		brill::JoinPath(g_ctx.assets_dir, "14O_d_6Li_2+_e_theta.txt"), "10C");
	ref_10C_2p->SetLineColor(kRed);
	ref_10C_2p->SetLineStyle(7);
	ref_10C_2p->SetLineWidth(2);
	TGraph *ref_6Li_2p = brill::LoadEThetaCurve(
		brill::JoinPath(g_ctx.assets_dir, "14O_d_6Li_2+_e_theta.txt"), "6Li");
	ref_6Li_2p->SetLineColor(kRed);
	ref_6Li_2p->SetLineStyle(7);
	ref_6Li_2p->SetLineWidth(2);
	TGraph *ref_theta_2p = brill::LoadThetaThetaCurve(
		brill::JoinPath(g_ctx.assets_dir, "14O_d_6Li_2+_theta_theta.txt"));
	ref_theta_2p->SetLineColor(kRed);
	ref_theta_2p->SetLineStyle(7);
	ref_theta_2p->SetLineWidth(2);

	ac.canvas->cd(1);
	ac.h_E_6Li_E_10C->Draw("colz");

	double total_E_line = g_ctx.entry_total_E
		? g_ctx.entry_total_E->GetNumber() : 0.0;
	if (total_E_line > 0) {
		TGraph *g_total = new TGraph(2);
		g_total->SetPoint(0, total_E_line - 80, 80);
		g_total->SetPoint(1, total_E_line - 30, 30);
		g_total->SetLineColor(kRed);
		g_total->SetLineWidth(2);
		g_total->Draw("L same");
	}

	ac.canvas->cd(2);
	ac.h_10C_e_theta->Draw("colz");
	ref_10C->Draw("L same");
	ref_10C_2p->Draw("L same");

	ac.canvas->cd(3);
	ac.h_6Li_e_theta->Draw("colz");
	ref_6Li->Draw("L same");
	ref_6Li_2p->Draw("L same");

	ac.canvas->cd(4);
	ac.h_theta_theta->Draw("colz");
	ref_theta->Draw("L same");
	ref_theta_2p->Draw("L same");

	ac.canvas->Modified();
	ac.canvas->Update();
}

static void DrawExcitationHistograms() {
	auto &ec = g_ctx.ec;
	if (!ec.canvas) return;
	ec.canvas->Clear();
	ec.canvas->Divide(2, 2);
	ec.canvas->cd(1);
	ec.h_excitation->Draw();
	ec.canvas->cd(2);
	ec.h_excitation_10C_mm->Draw();
	ec.canvas->cd(3);
	ec.h_excitation_6Li_mm->Draw();
	ec.canvas->Modified();
	ec.canvas->Update();
}

static void OnBeamChanged() {
	if (g_ctx.current_file.empty()) return;
	FillHistograms();
	DrawHistograms();

	std::string cut_path = "src/brill/Cut/cal_d1_d2_6Li_cut.C";
	TCutG *cut = brill::LoadCutGFromFile(cut_path);
	FillAnalysisHistograms(cut);
	DrawAnalysisHistograms();
	DrawExcitationHistograms();
}

static void OnBeamChanged();

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

	int run_min = g_ctx.entry_run_min ? g_ctx.entry_run_min->GetIntNumber() : 0;
	int run_max = g_ctx.entry_run_max ? g_ctx.entry_run_max->GetIntNumber() : 0;
	bool apply_run_filter = (run_min > 0 || run_max > 0);

	printf("  Loading %lld events from d_Li6 file...\n", n_entries);
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

	g_ctx.status_bar->SetText("Filling histograms...");
	gSystem->ProcessEvents();

	FillHistograms();
	DrawHistograms();

	std::string cut_path = "src/brill/Cut/cal_d1_d2_6Li_cut.C";
	printf("  Loading cut: %s\n", cut_path.c_str());
	TCutG *cut = brill::LoadCutGFromFile(cut_path);

	FillAnalysisHistograms(cut);
	DrawAnalysisHistograms();
	DrawExcitationHistograms();

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
	g_ctx.config = config;
	g_ctx.d_Li6_dir = brill::JoinPath(config.workspace, config.paths.d_Li6);
	g_ctx.assets_dir = "assets";
	brill::SetAssetsPath(config.assets);

	TApplication app("GUI_d_Li6", &argc, argv);
	gStyle->SetPalette(kRainBow);

	gInterpreter->Declare(
		TString::Format("volatile int &g_menu_action = *((volatile int*)%lu);",
			(unsigned long)&g_menu_action).Data()
	);
	gInterpreter->Declare("void HandleMenuSlot(Int_t id) { g_menu_action = (int)id; }");
	gInterpreter->Declare(
		TString::Format("volatile int &g_beam_filter = *((volatile int*)%lu);",
			(unsigned long)&g_beam_filter).Data());
	gInterpreter->Declare(
		TString::Format("volatile bool &g_beam_changed = *((volatile bool*)%lu);",
			(unsigned long)&g_beam_changed).Data());
	gInterpreter->Declare(
		TString::Format("volatile bool &g_redraw = *((volatile bool*)%lu);",
			(unsigned long)&g_redraw).Data());

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

	TGLabel *bins_label = new TGLabel(beam_frame, "  Bins: ");
	beam_frame->AddFrame(bins_label, new TGLayoutHints(kLHintsCenterY, 10, 2, 2, 2));

	TGNumberEntry *entry_excitation_bins = new TGNumberEntry(beam_frame, 30, 5, -1,
		TGNumberFormat::kNESInteger,
		TGNumberFormat::kNEANonNegative,
		TGNumberFormat::kNELLimitMinMax, 5, 500);
	TGLabel *bins_label2 = new TGLabel(beam_frame, "/bin");
	beam_frame->AddFrame(entry_excitation_bins, new TGLayoutHints(kLHintsCenterY, 2, 2, 2, 2));
	beam_frame->AddFrame(bins_label2, new TGLayoutHints(kLHintsCenterY, 2, 2, 2, 2));
	g_ctx.entry_excitation_bins = entry_excitation_bins;

	TGLabel *total_e_label = new TGLabel(beam_frame, "  E_{10C}+E_{6Li} > ");
	beam_frame->AddFrame(total_e_label, new TGLayoutHints(kLHintsCenterY, 10, 2, 2, 2));

	TGNumberEntry *entry_total_E = new TGNumberEntry(beam_frame, 0, 5, -1,
		TGNumberFormat::kNESReal,
		TGNumberFormat::kNEANonNegative,
		TGNumberFormat::kNELLimitMinMax, 0, 1000);
	beam_frame->AddFrame(entry_total_E, new TGLayoutHints(kLHintsCenterY, 2, 2, 2, 2));
	TGLabel *total_e_unit = new TGLabel(beam_frame, " MeV");
	beam_frame->AddFrame(total_e_unit, new TGLayoutHints(kLHintsCenterY, 2, 2, 2, 2));
	g_ctx.entry_total_E = entry_total_E;

	TGCheckButton *chk_loss = new TGCheckButton(beam_frame, "Eloss corr.");
	chk_loss->SetState(kButtonDown);
	beam_frame->AddFrame(chk_loss, new TGLayoutHints(kLHintsCenterY, 10, 2, 2, 2));
	g_ctx.chk_energy_loss = chk_loss;

	TGTextButton *btn_draw = new TGTextButton(beam_frame, "Draw");
	btn_draw->SetCommand("g_redraw = true;");
	beam_frame->AddFrame(btn_draw, new TGLayoutHints(kLHintsCenterY, 10, 2, 2, 2));
	g_ctx.btn_draw = btn_draw;

	main_frame->AddFrame(beam_frame, new TGLayoutHints(kLHintsTop | kLHintsLeft, 4, 4, 2, 2));

	TRootEmbeddedCanvas *embed = new TRootEmbeddedCanvas("embed_main", main_frame, 1200, 700);
	main_frame->AddFrame(embed, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));
	g_ctx.bcs.embed = embed;
	g_ctx.bcs.canvas = embed->GetCanvas();

	g_ctx.ac.canvas = new TCanvas("canvas_analysis", "d+6Li Analysis", 1200, 800);
	g_ctx.ec.canvas = new TCanvas("canvas_excitation", "10C Excitation Energy", 800, 600);

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
		if (g_beam_changed) {
			g_beam_changed = false;
			OnBeamChanged();
		}
		if (g_redraw) {
			g_redraw = false;
			if (!g_ctx.current_file.empty()) {
				OnFileOpen();
			}
		}
		g_menu_action = 0;
	}

	return 0;
}