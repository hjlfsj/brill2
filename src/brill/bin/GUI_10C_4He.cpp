#include "include/config.h"
#include "include/10C+4He/10C_4He_event.h"
#include "include/energy_calculator/range_energy_calculator.h"
#include "include/physics/kinematics.h"
#include "include/rebuild/nuclear_data.h"
#include "include/utils.h"
#include "external/cxxopts.hpp"

#include <catima/catima.h>
#include <TLorentzVector.h>
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
#include <TF1.h>
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
#include <fstream>
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

struct ExcitationCanvas {
	TCanvas *canvas = nullptr;
	TH1D *h_excitation_10C = nullptr;
	TH1D *h_14O_gate1 = nullptr;
	TH1D *h_14O_gate2 = nullptr;
	TH2D *h_recoil_ep = nullptr;
};

struct ExcitationCanvasFixed {
	TCanvas *canvas = nullptr;
	TH1D *h_excitation_10C = nullptr;
	TH1D *h_14O_gate1 = nullptr;
	TH1D *h_14O_gate2 = nullptr;
	TH2D *h_recoil_ep = nullptr;
};

struct Window6Canvas {
	TCanvas *canvas = nullptr;
	TH1D *h_missing_mass = nullptr;
	TH1D *h_inversion = nullptr;
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
	TGNumberEntry *entry_x1 = nullptr;
	TGNumberEntry *entry_x2 = nullptr;
	TGNumberEntry *entry_x3 = nullptr;
	TGNumberEntry *entry_z = nullptr;

	BeamCanvases bcs;
	FilteredCanvases fcs;
	AnalysisCanvas ac;
	ExcitationCanvas ec;
	ExcitationCanvasFixed ec_fixed;

	Window6Canvas w6;

	std::string current_file;
	std::string config_path;
	std::string c10_he4_dir;
	std::string beam_dir;
	std::string trigger;
	bool is_t1 = false;
	brill::AppConfig config;
	std::vector<brill::C10He4Event> all_events;

	double z_d1 = 0.0;
	double z_d2 = 0.0;
	double z_d3 = 0.0;
	double z_d4 = 0.0;
};

struct BeamFileCache {
	TFile *file = nullptr;
	TTree *tree = nullptr;
	double E = -1.0;
};

static std::map<std::pair<int, bool>, BeamFileCache> g_beam_cache;

static double ReadBeamEnergy(int run, Long64_t entry, const std::string &beam_dir, bool use_t1) {
	auto &cache = g_beam_cache[{run, use_t1}];

	if (!cache.tree) {
		TString beam_path = TString::Format("%s/beam_%s%04d.root",
			beam_dir.c_str(),
			use_t1 ? "t1_" : "",
			run);
		if (!std::filesystem::exists(beam_path.Data())) {
			fprintf(stderr, "Warning: beam file not found: %s\n", beam_path.Data());
			return -1.0;
		}
		cache.file = new TFile(beam_path, "read");
		cache.tree = (TTree*)cache.file->Get("tree");
		if (!cache.tree) {
			fprintf(stderr, "Warning: no tree in %s\n", beam_path.Data());
			delete cache.file;
			cache.file = nullptr;
			return -1.0;
		}
		if (!cache.tree->GetBranch("E")) {
			fprintf(stderr, "Warning: no E branch in %s\n", beam_path.Data());
			delete cache.file;
			cache.file = nullptr;
			cache.tree = nullptr;
			return -1.0;
		}
		cache.tree->SetBranchAddress("E", &cache.E);
	}

	Long64_t n_entries = cache.tree->GetEntries();
	if (entry < 0 || entry >= n_entries) {
		return -1.0;
	}

	cache.tree->GetEntry(entry);
	return cache.E;
}

static std::map<int, double> g_fixed_beam_cache;

static double ReadFixedBeamEnergy(int run, const std::string &beam_dir) {
	auto it = g_fixed_beam_cache.find(run);
	if (it != g_fixed_beam_cache.end()) return it->second;

	TString path = TString::Format("%s/check_tof_and_beam_%04d.txt",
		beam_dir.c_str(), run);
	std::ifstream fin(path.Data());
	if (!fin.is_open()) {
		fprintf(stderr,
			"Error: per-run beam energy file not found: %s\n"
			"  Please run 'check_tof_and_beam -r %d' first to generate it.\n",
			path.Data(), run);
		return -1.0;
	}

	std::string header;
	std::getline(fin, header);
	double E0 = -1.0, sigma_E = 0.0, t0 = 0.0, sigma_t0 = 0.0, dt = 0.0;
	fin >> E0 >> sigma_E >> t0 >> sigma_t0 >> dt;
	fin.close();

	if (E0 <= 0) {
		fprintf(stderr, "Error: invalid beam energy E0=%.4f in %s\n",
			E0, path.Data());
		return -1.0;
	}

	g_fixed_beam_cache[run] = E0;
	return E0;
}

static GUIContext g_ctx;
static volatile int g_menu_action = 0;
static volatile int g_beam_filter = 0;
static volatile bool g_beam_changed = false;
static volatile bool g_redraw = false;
static volatile bool g_is_track = false;

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
	int nbins = g_ctx.entry_excitation_bins
		? g_ctx.entry_excitation_bins->GetIntNumber() : 100;
	if (nbins < 10) nbins = 10;
	double bin_keV = 20000.0 / nbins;
	ac.h_excitation = new TH1D(
		"h_excitation",
		TString::Format("^{14}O Excitation Energy%s;E_{x} (MeV);Counts / %.0f keV",
			label.Data(), bin_keV),
		nbins, 5, 25);
	ac.h_excitation->SetDirectory(0);
}

static void RebuildExcitationHistogram() {
	auto &ec = g_ctx.ec;
	if (ec.h_excitation_10C) delete ec.h_excitation_10C;
	if (ec.h_14O_gate1) delete ec.h_14O_gate1;
	if (ec.h_14O_gate2) delete ec.h_14O_gate2;
	if (ec.h_recoil_ep) delete ec.h_recoil_ep;

	ec.h_excitation_10C = new TH1D(
		"h_excitation_10C",
		"^{10}C Excitation Energy (^{14}O beam);E_{x}(^{10}C) (MeV);Counts",
		100, -20, 80);
	ec.h_excitation_10C->SetDirectory(0);

	ec.h_recoil_ep = new TH2D(
		"h_recoil_ep",
		"Recoil E-P (^{14}O beam);p^{2}_{D} / 2u (MeV);T_{beam} - T_{10C} - T_{4He} (MeV)",
		200, 0, 200, 200, 0, 200);
	ec.h_recoil_ep->SetDirectory(0);

	int nbins = g_ctx.entry_excitation_bins
		? g_ctx.entry_excitation_bins->GetIntNumber() : 100;
	if (nbins < 10) nbins = 10;
	double bin_keV = 20000.0 / nbins;

	ec.h_14O_gate1 = new TH1D(
		"h_14O_gate1",
		TString::Format("^{14}O E_{x} (^{10}C E_{x} gate [%.0f, %.0f]);E_{x}(^{14}O) (MeV);Counts / %.0f keV",
			g_ctx.entry_x1 ? g_ctx.entry_x1->GetNumber() : -2,
			g_ctx.entry_x2 ? g_ctx.entry_x2->GetNumber() : 2,
			bin_keV),
		nbins, 5, 25);
	ec.h_14O_gate1->SetDirectory(0);

	ec.h_14O_gate2 = new TH1D(
		"h_14O_gate2",
		TString::Format("^{14}O E_{x} (^{10}C E_{x} gate [%.0f, %.0f]);E_{x}(^{14}O) (MeV);Counts / %.0f keV",
			g_ctx.entry_x2 ? g_ctx.entry_x2->GetNumber() : 2,
			g_ctx.entry_x3 ? g_ctx.entry_x3->GetNumber() : 7,
			bin_keV),
		nbins, 5, 25);
	ec.h_14O_gate2->SetDirectory(0);

	auto &ecf = g_ctx.ec_fixed;
	if (ecf.h_excitation_10C) delete ecf.h_excitation_10C;
	if (ecf.h_14O_gate1) delete ecf.h_14O_gate1;
	if (ecf.h_14O_gate2) delete ecf.h_14O_gate2;
	if (ecf.h_recoil_ep) delete ecf.h_recoil_ep;

	ecf.h_excitation_10C = new TH1D(
		"h_excitation_10C_fixed",
		"^{10}C E_{x} (^{14}O beam, E0 fixed);E_{x}(^{10}C) (MeV);Counts",
		100, -20, 80);
	ecf.h_excitation_10C->SetDirectory(0);

	ecf.h_recoil_ep = new TH2D(
		"h_recoil_ep_fixed",
		"Recoil E-P (^{14}O beam, E0 fixed);p^{2}_{D} / 2u (MeV);T_{beam} - T_{10C} - T_{4He} (MeV)",
		200, 0, 200, 200, 0, 200);
	ecf.h_recoil_ep->SetDirectory(0);

	ecf.h_14O_gate1 = new TH1D(
		"h_14O_gate1_fixed",
		TString::Format("^{14}O E_{x} (^{10}C gate [%.0f, %.0f], E0 fixed);E_{x}(^{14}O) (MeV);Counts / %.0f keV",
			g_ctx.entry_x1 ? g_ctx.entry_x1->GetNumber() : -2,
			g_ctx.entry_x2 ? g_ctx.entry_x2->GetNumber() : 2,
			bin_keV),
		nbins, 5, 25);
	ecf.h_14O_gate1->SetDirectory(0);

	ecf.h_14O_gate2 = new TH1D(
		"h_14O_gate2_fixed",
		TString::Format("^{14}O E_{x} (^{10}C gate [%.0f, %.0f], E0 fixed);E_{x}(^{14}O) (MeV);Counts / %.0f keV",
			g_ctx.entry_x2 ? g_ctx.entry_x2->GetNumber() : 2,
			g_ctx.entry_x3 ? g_ctx.entry_x3->GetNumber() : 7,
			bin_keV),
		nbins, 5, 25);
	ecf.h_14O_gate2->SetDirectory(0);
}

static void RebuildWindow6Histogram() {
	auto &w6 = g_ctx.w6;
	if (w6.h_missing_mass) delete w6.h_missing_mass;
	if (w6.h_inversion) delete w6.h_inversion;

	w6.h_missing_mass = new TH1D(
		"h_missing_mass",
		"^{14}O E_{x} (missing mass via P_{d}+P_{10C}+P_{4He});E_{x}(^{14}O) (MeV);Counts",
		200, -10, 40);
	w6.h_missing_mass->SetDirectory(0);

	w6.h_inversion = new TH1D(
		"h_inversion",
		"^{10}C E_{x} (inversion from missing-mass E_{x}(^{14}O));E_{x}(^{10}C) (MeV);Counts",
		200, -10, 40);
	w6.h_inversion->SetDirectory(0);
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

static void ProjectToZ(
	double target_x, double target_y,
	double d2_x, double d2_y, double d2_z,
	double z_layer,
	double &pred_x, double &pred_y
) {
	double t = z_layer / d2_z;
	pred_x = target_x + t * (d2_x - target_x);
	pred_y = target_y + t * (d2_y - target_y);
}

static bool PassTrackCut(const brill::C10He4Event &ev) {
	if (!ev.ppac_valid) return false;

	double pred_x, pred_y;

	ProjectToZ(ev.target_x, ev.target_y,
		ev.t0d2_10C_x, ev.t0d2_10C_y, ev.t0d2_10C_z,
		g_ctx.z_d1, pred_x, pred_y);
	if (std::abs(pred_x - ev.t0d1_10C_x) > 2.0 || std::abs(pred_y - ev.t0d1_10C_y) > 2.0)
		return false;

	ProjectToZ(ev.target_x, ev.target_y,
		ev.t0d2_10C_x, ev.t0d2_10C_y, ev.t0d2_10C_z,
		g_ctx.z_d3, pred_x, pred_y);
	if (std::abs(pred_x - ev.t0d3_10C_x) > 2.0 || std::abs(pred_y - ev.t0d3_10C_y) > 2.0)
		return false;

	ProjectToZ(ev.target_x, ev.target_y,
		ev.t0d2_4He_x, ev.t0d2_4He_y, ev.t0d2_4He_z,
		g_ctx.z_d1, pred_x, pred_y);
	bool has_d1_4He = (ev.t0d1_4He_x != 0.0 || ev.t0d1_4He_y != 0.0);
	if (has_d1_4He) {
		if (std::abs(pred_x - ev.t0d1_4He_x) > 2.0 || std::abs(pred_y - ev.t0d1_4He_y) > 2.0)
			return false;
	}

	ProjectToZ(ev.target_x, ev.target_y,
		ev.t0d2_4He_x, ev.t0d2_4He_y, ev.t0d2_4He_z,
		g_ctx.z_d3, pred_x, pred_y);
	if (std::abs(pred_x - ev.t0d3_4He_x) > 2.0 || std::abs(pred_y - ev.t0d3_4He_y) > 2.0)
		return false;

	ProjectToZ(ev.target_x, ev.target_y,
		ev.t0d2_4He_x, ev.t0d2_4He_y, ev.t0d2_4He_z,
		g_ctx.z_d4, pred_x, pred_y);
	if (std::abs(pred_x - ev.t0d4_4He_x) > 2.0 || std::abs(pred_y - ev.t0d4_4He_y) > 2.0)
		return false;

	return true;
}

static void FillFilteredAndAnalysis(TCutG *cut) {
	RebuildFilteredHistograms();
	RebuildAnalysisHistograms();

	if (g_beam_filter == 1) {
		RebuildExcitationHistogram();
		RebuildWindow6Histogram();
	}

	if (g_ctx.all_events.empty()) return;

	double m_10C = brill::GetMass(6, 10);
	double m_4He = brill::GetMass(2, 4);

	catima::Material cd2_mat;
	cd2_mat.add_element(12.0, 6.0, 1.0);
	cd2_mat.add_element(2.0, 1.0, 2.0);
	cd2_mat.density(0.79867);
	std::string cache_dir = brill::JoinPath(g_ctx.config.workspace, g_ctx.config.paths.energy_calculator);
	TString beam_cache_path = TString::Format("%s/target_z8_a14_cd2.root", cache_dir.c_str());
	brill::RangeEnergyCalculator beam_calc(8, 14, cd2_mat, beam_cache_path.Data());
	TString cache_10C = TString::Format("%s/target_z6_a10_cd2.root", cache_dir.c_str());
	brill::RangeEnergyCalculator calc_10C(6, 10, cd2_mat, cache_10C.Data());
	TString cache_4He = TString::Format("%s/target_z2_a4_cd2.root", cache_dir.c_str());
	brill::RangeEnergyCalculator calc_4He(2, 4, cd2_mat, cache_4He.Data());
	const double half_thick_um = 122.0 / 2.0;

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
		if (g_is_track && !PassTrackCut(ev)) continue;

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

		bool filled_14O_ex = false;

		if (g_beam_filter == 1) {
			static long long ex_cnt = 0;
			double E_beam_raw = ReadBeamEnergy(ev.run_number, ev.entry, g_ctx.beam_dir, g_ctx.is_t1);
			if (E_beam_raw > 0) {
				ex_cnt++;

				double norm_beam = std::sqrt(1.0 + ev.dir_x * ev.dir_x + ev.dir_y * ev.dir_y);
				double z = g_ctx.entry_z->GetNumber();
				double beam_depth = (half_thick_um + z) * norm_beam;
				double exit_thick = std::max(0.0, half_thick_um - z);

				double beam_range_before = beam_calc.Range(E_beam_raw);
				double beam_range_rxn = beam_range_before - beam_depth;
				double E_beam = beam_range_rxn > 0 ? beam_calc.Energy(beam_range_rxn) : E_beam_raw;

				double costh_10C = std::cos(ev.theta_10C * M_PI / 180.0);
				double costh_4He = std::cos(ev.theta_4He * M_PI / 180.0);
				double exit_10C = exit_thick / std::max(costh_10C, 0.001);
				double exit_4He = exit_thick / std::max(costh_4He, 0.001);

				double r_10C_before = calc_10C.Range(T_10C);
				double r_10C_rxn = r_10C_before + exit_10C;
				double T_10C_corr = calc_10C.Energy(r_10C_rxn);

				double r_4He_before = calc_4He.Range(T_4He);
				double r_4He_rxn = r_4He_before + exit_4He;
				double T_4He_corr = calc_4He.Energy(r_4He_rxn);

				double E_10C = T_10C_corr + m_10C;
				double E_4He = T_4He_corr + m_4He;

				double p_10C = std::sqrt(E_10C * E_10C - m_10C * m_10C);
				double p_4He = std::sqrt(E_4He * E_4He - m_4He * m_4He);

				double dx_10C = ev.t0d2_10C_x - ev.target_x;
				double dy_10C = ev.t0d2_10C_y - ev.target_y;
				double dz_10C = ev.t0d2_10C_z;
				double norm_10C = std::sqrt(dx_10C*dx_10C + dy_10C*dy_10C + dz_10C*dz_10C);

				double dx_4He = ev.t0d2_4He_x - ev.target_x;
				double dy_4He = ev.t0d2_4He_y - ev.target_y;
				double dz_4He = ev.t0d2_4He_z;
				double norm_4He = std::sqrt(dx_4He*dx_4He + dy_4He*dy_4He + dz_4He*dz_4He);

				TLorentzVector P_10C(p_10C * dx_10C / norm_10C,
					p_10C * dy_10C / norm_10C,
					p_10C * dz_10C / norm_10C, E_10C);
				TLorentzVector P_4He(p_4He * dx_4He / norm_4He,
					p_4He * dy_4He / norm_4He,
					p_4He * dz_4He / norm_4He, E_4He);

				double m_14O = brill::GetMass(8, 14);
				double m_d = brill::GetMass(1, 2);

				double E_A = E_beam + m_14O;
				double p_A_mag = std::sqrt(E_A * E_A - m_14O * m_14O);

				TVector3 beam_dir(ev.dir_x / norm_beam, ev.dir_y / norm_beam, 1.0 / norm_beam);

				TVector3 p_A_vec(p_A_mag * beam_dir.X(),
					p_A_mag * beam_dir.Y(),
					p_A_mag * beam_dir.Z());

				TVector3 p_D = p_A_vec - P_10C.Vect() - P_4He.Vect();
				double p_D_sq = p_D.Mag2();
				double T_D = std::sqrt(p_D_sq + m_d * m_d) - m_d;

				TLorentzVector P_A(p_A_vec, E_A);
				TLorentzVector P_d_target(0, 0, 0, m_d);
				TLorentzVector P_D_4(p_D, m_d + T_D);
				TLorentzVector P_14O_miss = P_A + P_d_target - P_D_4;
				double E_x_14O_mm = P_14O_miss.M() - m_14O;
				auto &w6 = g_ctx.w6;
				w6.h_missing_mass->Fill(E_x_14O_mm);

				double M_inv_mm = m_14O + E_x_14O_mm;
				double M0_sq = (P_10C + P_4He).M2();
				TVector3 p_10C_vec = P_10C.Vect();
				TVector3 p_4He_vec = P_4He.Vect();
				TVector3 p_10C_hat = p_10C_vec.Unit();
				double denom = 2.0 * (E_10C + E_4He - (p_10C_vec + p_4He_vec).Dot(p_10C_hat));
				double E_gamma_inv = (M_inv_mm * M_inv_mm - M0_sq) / denom;
				w6.h_inversion->Fill(E_gamma_inv);

				double Q = T_10C_corr + T_4He_corr + T_D - E_beam;
				double Q0 = m_14O - m_4He - m_10C;
				double E_x_10C = Q0 - Q;

				g_ctx.ec.h_excitation_10C->Fill(E_x_10C);
				g_ctx.ec.h_recoil_ep->Fill(
					p_D_sq / (2.0 * 931.494),
					E_beam - T_10C_corr - T_4He_corr);

				double x1 = g_ctx.entry_x1->GetNumber();
				double x2 = g_ctx.entry_x2->GetNumber();
				double x3 = g_ctx.entry_x3->GetNumber();

				bool in_gate1 = (E_x_10C >= x1 && E_x_10C <= x2);
				bool in_gate2 = (E_x_10C >= x2 && E_x_10C <= x3);

				double gamma_E = -1.0;
				if (in_gate1) {
					gamma_E = 0.0;
				} else if (in_gate2) {
					gamma_E = 3.35;
				}

				if (gamma_E >= 0.0) {
					TVector3 p_10C_hat = P_10C.Vect().Unit();
					TLorentzVector P_gamma(
						p_10C_hat.X() * gamma_E,
						p_10C_hat.Y() * gamma_E,
						p_10C_hat.Z() * gamma_E,
						gamma_E);
					TLorentzVector P_10C_star = P_10C + P_gamma;
					TLorentzVector P_14O_reco = P_10C_star + P_4He;
					double E_x_14O_corr = P_14O_reco.M() - m_14O;

					ac.h_excitation->Fill(E_x_14O_corr);
					filled_14O_ex = true;

					if (in_gate1)
						g_ctx.ec.h_14O_gate1->Fill(E_x_14O_corr);
					if (in_gate2)
						g_ctx.ec.h_14O_gate2->Fill(E_x_14O_corr);
				}

				double E_beam_fixed = ReadFixedBeamEnergy(ev.run_number, g_ctx.beam_dir);
				if (E_beam_fixed > 0) {
					double beam_range_before_f = beam_calc.Range(E_beam_fixed);
					double beam_range_rxn_f = beam_range_before_f - beam_depth;
					double E_beam_f = beam_range_rxn_f > 0
						? beam_calc.Energy(beam_range_rxn_f) : E_beam_fixed;

					double E_A_f = E_beam_f + m_14O;
					double p_A_mag_f = std::sqrt(E_A_f * E_A_f - m_14O * m_14O);
					TVector3 p_A_vec_f(p_A_mag_f * beam_dir.X(),
						p_A_mag_f * beam_dir.Y(),
						p_A_mag_f * beam_dir.Z());
					TVector3 p_D_f = p_A_vec_f - P_10C.Vect() - P_4He.Vect();
					double p_D_sq_f = p_D_f.Mag2();
					double T_D_f = std::sqrt(p_D_sq_f + m_d * m_d) - m_d;

					double Q_f = T_10C_corr + T_4He_corr + T_D_f - E_beam_f;
					double E_x_10C_f = Q0 - Q_f;

					auto &ecf = g_ctx.ec_fixed;
					ecf.h_excitation_10C->Fill(E_x_10C_f);
					ecf.h_recoil_ep->Fill(
						p_D_sq_f / (2.0 * 931.494),
						E_beam_f - T_10C_corr - T_4He_corr);

					bool in_gate1_f = (E_x_10C_f >= x1 && E_x_10C_f <= x2);
					bool in_gate2_f = (E_x_10C_f >= x2 && E_x_10C_f <= x3);

					double gamma_E_f = -1.0;
					if (in_gate1_f) {
						gamma_E_f = 0.0;
					} else if (in_gate2_f) {
						gamma_E_f = 3.35;
					}

					if (gamma_E_f >= 0.0) {
						TVector3 p_10C_hat_f = P_10C.Vect().Unit();
						TLorentzVector P_gamma_f(
							p_10C_hat_f.X() * gamma_E_f,
							p_10C_hat_f.Y() * gamma_E_f,
							p_10C_hat_f.Z() * gamma_E_f,
							gamma_E_f);
						TLorentzVector P_10C_star_f = P_10C + P_gamma_f;
						TLorentzVector P_14O_reco_f = P_10C_star_f + P_4He;
						double E_x_14O_corr_f = P_14O_reco_f.M() - m_14O;

						if (in_gate1_f)
							ecf.h_14O_gate1->Fill(E_x_14O_corr_f);
						if (in_gate2_f)
							ecf.h_14O_gate2->Fill(E_x_14O_corr_f);
					}
				}

				if (ex_cnt % 10000 == 0)
					printf("\r  [10CEx] event=%lld E_x(10C)=%.2f MeV",
						ex_cnt, E_x_10C);
			}
		}

		if (!filled_14O_ex) {
			ac.h_excitation->Fill(E_x);
		}
	}

	printf("\r  Filtering done: %d total, passed=%d, analyzed=%d        \n",
		total, passed, analyzed);
	g_ctx.status_bar->SetText(
		TString::Format("Beam=%s%s: %d passed, %d analyzed. File: %s",
			BeamLabel(), g_is_track ? "+Track" : "",
			passed, analyzed, g_ctx.current_file.c_str()));
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
	auto &ec = g_ctx.ec;
	if (ec.h_recoil_ep) {
		ec.h_recoil_ep->Draw("colz");
		TF1 *f1 = new TF1("f_y_eq_x", "x + 10.12", 0, 200);
		f1->SetLineColor(kRed);
		f1->SetLineWidth(2);
		f1->Draw("SAME");
		TF1 *f2 = new TF1("f_y_eq_halfx", "0.5*x + 10.12", 0, 200);
		f2->SetLineColor(kBlue);
		f2->SetLineWidth(2);
		f2->Draw("SAME");
	}
	ac.canvas->Modified();
	ac.canvas->Update();
}

static void DrawExcitationHistogram() {
	auto &ec = g_ctx.ec;
	if (!ec.canvas) return;
	if (!ec.h_excitation_10C) return;
	ec.canvas->Clear();
	ec.canvas->Divide(2, 2);

	ec.canvas->cd(1);
	ec.h_excitation_10C->Draw();

	ec.canvas->cd(2);
	if (ec.h_14O_gate1 && ec.h_14O_gate2) {
		ec.h_14O_gate1->SetLineColor(kRed);
		ec.h_14O_gate1->SetMarkerColor(kRed);
		ec.h_14O_gate2->SetLineColor(kBlue);
		ec.h_14O_gate2->SetMarkerColor(kBlue);
		double ymax = std::max(ec.h_14O_gate1->GetMaximum(),
			ec.h_14O_gate2->GetMaximum()) * 1.15;
		ec.h_14O_gate1->SetMaximum(ymax);
		ec.h_14O_gate1->Draw();
		ec.h_14O_gate2->Draw("SAME");
	} else if (ec.h_14O_gate1) {
		ec.h_14O_gate1->Draw();
	}

	ec.canvas->cd(3);
	if (ec.h_14O_gate1) ec.h_14O_gate1->Draw();

	ec.canvas->cd(4);
	if (ec.h_14O_gate2) ec.h_14O_gate2->Draw();

	ec.canvas->Modified();
	ec.canvas->Update();
}

static void DrawExcitationFixedHistogram() {
	auto &ecf = g_ctx.ec_fixed;
	if (!ecf.canvas) return;
	if (!ecf.h_excitation_10C) return;
	ecf.canvas->Clear();
	ecf.canvas->Divide(2, 2);

	ecf.canvas->cd(1);
	ecf.h_excitation_10C->Draw();

	ecf.canvas->cd(2);
	if (ecf.h_14O_gate1 && ecf.h_14O_gate2) {
		ecf.h_14O_gate1->SetLineColor(kRed);
		ecf.h_14O_gate1->SetMarkerColor(kRed);
		ecf.h_14O_gate2->SetLineColor(kBlue);
		ecf.h_14O_gate2->SetMarkerColor(kBlue);
		double ymax = std::max(ecf.h_14O_gate1->GetMaximum(),
			ecf.h_14O_gate2->GetMaximum()) * 1.15;
		ecf.h_14O_gate1->SetMaximum(ymax);
		ecf.h_14O_gate1->Draw();
		ecf.h_14O_gate2->Draw("SAME");
	} else if (ecf.h_14O_gate1) {
		ecf.h_14O_gate1->Draw();
	}

	ecf.canvas->cd(3);
	if (ecf.h_14O_gate1) ecf.h_14O_gate1->Draw();

	ecf.canvas->cd(4);
	if (ecf.h_14O_gate2) ecf.h_14O_gate2->Draw();

	ecf.canvas->Modified();
	ecf.canvas->Update();
}

static void DrawWindow6Histogram() {
	auto &w6 = g_ctx.w6;
	if (!w6.canvas) return;
	if (!w6.h_missing_mass) return;
	w6.canvas->Clear();
	w6.canvas->Divide(2, 2);

	w6.canvas->cd(1);
	w6.h_missing_mass->Draw();

	w6.canvas->cd(2);
	w6.h_inversion->Draw();

	w6.canvas->Modified();
	w6.canvas->Update();
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

	if (g_beam_filter == 1) {
		DrawExcitationHistogram();
		DrawExcitationFixedHistogram();
		DrawWindow6Histogram();
	}
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

	g_ctx.is_t1 = (filepath.find("_t1_") != std::string::npos);

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

	if (g_beam_filter == 1) {
		DrawExcitationHistogram();
		DrawExcitationFixedHistogram();
		DrawWindow6Histogram();
	}

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
	g_ctx.config = config;
	g_ctx.c10_he4_dir = brill::JoinPath(config.workspace, config.paths.c10_he4);
	g_ctx.beam_dir = brill::JoinPath(config.workspace, config.paths.beam);
	g_ctx.trigger = config.trigger;

	const auto *det_d1 = brill::FindDetectorConfig(config, "t0d1");
	const auto *det_d2 = brill::FindDetectorConfig(config, "t0d2");
	const auto *det_d3 = brill::FindDetectorConfig(config, "t0d3");
	const auto *det_d4 = brill::FindDetectorConfig(config, "t0d4");
	if (det_d1) g_ctx.z_d1 = det_d1->z_mm;
	if (det_d2) g_ctx.z_d2 = det_d2->z_mm;
	if (det_d3) g_ctx.z_d3 = det_d3->z_mm;
	if (det_d4) g_ctx.z_d4 = det_d4->z_mm;
	printf("T0 Z positions: d1=%.1f d2=%.1f d3=%.1f d4=%.1f mm\n",
		g_ctx.z_d1, g_ctx.z_d2, g_ctx.z_d3, g_ctx.z_d4);

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
	gInterpreter->Declare(
		TString::Format("volatile bool &g_is_track = *((volatile bool*)%lu);",
			(unsigned long)&g_is_track).Data()
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

	TGLabel *bins_label = new TGLabel(beam_frame, "  Bins: ");
	beam_frame->AddFrame(bins_label, new TGLayoutHints(kLHintsCenterY, 10, 2, 2, 2));

	TGNumberEntry *entry_excitation_bins = new TGNumberEntry(beam_frame, 100, 5, -1,
		TGNumberFormat::kNESInteger,
		TGNumberFormat::kNEANonNegative,
		TGNumberFormat::kNELLimitMinMax, 10, 1000);
	TGLabel *bins_label2 = new TGLabel(beam_frame, "/bin");
	beam_frame->AddFrame(entry_excitation_bins, new TGLayoutHints(kLHintsCenterY, 2, 2, 2, 2));
	beam_frame->AddFrame(bins_label2, new TGLayoutHints(kLHintsCenterY, 2, 2, 2, 2));
	g_ctx.entry_excitation_bins = entry_excitation_bins;

	TGTextButton *btn_draw = new TGTextButton(beam_frame, "Draw");
	btn_draw->SetCommand("g_redraw = true;");
	beam_frame->AddFrame(btn_draw, new TGLayoutHints(kLHintsCenterY, 10, 2, 2, 2));
	g_ctx.btn_draw = btn_draw;

	TGCheckButton *chk_track = new TGCheckButton(beam_frame, "Track");
	chk_track->SetCommand("g_is_track = !g_is_track; g_beam_changed = true;");
	beam_frame->AddFrame(chk_track, new TGLayoutHints(kLHintsCenterY, 10, 2, 2, 2));

	TGLabel *x1_label = new TGLabel(beam_frame, " x1: ");
	beam_frame->AddFrame(x1_label, new TGLayoutHints(kLHintsCenterY, 10, 2, 2, 2));
	g_ctx.entry_x1 = new TGNumberEntry(beam_frame, -2, 5, -1,
		TGNumberFormat::kNESReal,
		TGNumberFormat::kNEAAnyNumber,
		TGNumberFormat::kNELLimitMinMax, -20, 40);
	beam_frame->AddFrame(g_ctx.entry_x1, new TGLayoutHints(kLHintsCenterY, 2, 2, 2, 2));

	TGLabel *x2_label = new TGLabel(beam_frame, " x2: ");
	beam_frame->AddFrame(x2_label, new TGLayoutHints(kLHintsCenterY, 5, 2, 2, 2));
	g_ctx.entry_x2 = new TGNumberEntry(beam_frame, 2, 5, -1,
		TGNumberFormat::kNESReal,
		TGNumberFormat::kNEANonNegative,
		TGNumberFormat::kNELLimitMinMax, -20, 40);
	beam_frame->AddFrame(g_ctx.entry_x2, new TGLayoutHints(kLHintsCenterY, 2, 2, 2, 2));

	TGLabel *x3_label = new TGLabel(beam_frame, " x3: ");
	beam_frame->AddFrame(x3_label, new TGLayoutHints(kLHintsCenterY, 5, 2, 2, 2));
	g_ctx.entry_x3 = new TGNumberEntry(beam_frame, 7, 5, -1,
		TGNumberFormat::kNESReal,
		TGNumberFormat::kNEANonNegative,
		TGNumberFormat::kNELLimitMinMax, -20, 40);
	beam_frame->AddFrame(g_ctx.entry_x3, new TGLayoutHints(kLHintsCenterY, 2, 2, 2, 2));

	TGLabel *z_label = new TGLabel(beam_frame, " z(µm): ");
	beam_frame->AddFrame(z_label, new TGLayoutHints(kLHintsCenterY, 10, 2, 2, 2));
	g_ctx.entry_z = new TGNumberEntry(beam_frame, 0, 5, -1,
		TGNumberFormat::kNESReal,
		TGNumberFormat::kNEAAnyNumber,
		TGNumberFormat::kNELLimitMinMax, -61, 61);
	beam_frame->AddFrame(g_ctx.entry_z, new TGLayoutHints(kLHintsCenterY, 2, 2, 2, 2));

	main_frame->AddFrame(beam_frame, new TGLayoutHints(kLHintsTop | kLHintsLeft, 4, 4, 2, 2));

	TRootEmbeddedCanvas *embed = new TRootEmbeddedCanvas("embed_main", main_frame, 1200, 700);
	main_frame->AddFrame(embed, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));
	g_ctx.bcs.embed = embed;
	g_ctx.bcs.canvas = embed->GetCanvas();

	g_ctx.fcs.canvas = new TCanvas("canvas_filtered", "10C+4He Filtered", 1200, 800);
	g_ctx.ac.canvas = new TCanvas("canvas_analysis", "10C+4He Analysis", 1200, 800);
	g_ctx.ec.canvas = new TCanvas("canvas_excitation", "10C Excitation Energy (per-event E_beam)", 800, 600);
	g_ctx.ec_fixed.canvas = new TCanvas("canvas_excitation_fixed", "10C Excitation Energy (per-run E0)", 800, 600);
	g_ctx.w6.canvas = new TCanvas("canvas_window6", "Missing Mass + Inversion", 800, 600);

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