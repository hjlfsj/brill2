#include "include/config.h"
#include "include/d_6Li/d_6Li_event.h"
#include "include/energy_calculator/range_energy_calculator.h"
#include "include/Lise++/e_theta.h"
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
#include <fstream>
#include <TGraph.h>
#include <TH2D.h>
#include <TPad.h>
#include <TRootEmbeddedCanvas.h>
#include <TLorentzVector.h>
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
#include <utility>
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
	TH1D *h_excitation_mm = nullptr;
	TH1D *h_excitation_invmass = nullptr;
	TH1D *h_excitation_mm_fixed = nullptr;
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
	TGNumberEntry *entry_z = nullptr;

	BeamCanvases bcs;
	AnalysisCanvas ac;
	ExcitationCanvas ec;

	std::string current_file;
	std::string config_path;
	std::string d_Li6_dir;
	std::string beam_dir;
	std::string assets_dir;
	bool is_t1 = false;
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

struct BeamFileCache {
	TFile *file = nullptr;
	TTree *tree = nullptr;
	double E_beam = 0.0;
	int run = -1;
	~BeamFileCache() { if (file) { file->Close(); delete file; } }

	BeamFileCache() = default;
	BeamFileCache(BeamFileCache&& other) noexcept
		: file(other.file), tree(other.tree), E_beam(other.E_beam), run(other.run) {
		other.file = nullptr;
		other.tree = nullptr;
	}
	BeamFileCache& operator=(BeamFileCache&& other) noexcept {
		if (this != &other) {
			if (file) { file->Close(); delete file; }
			file = other.file;
			tree = other.tree;
			E_beam = other.E_beam;
			run = other.run;
			other.file = nullptr;
			other.tree = nullptr;
		}
		return *this;
	}
	BeamFileCache(const BeamFileCache&) = delete;
	BeamFileCache& operator=(const BeamFileCache&) = delete;
};
static std::map<std::pair<int, bool>, BeamFileCache> g_beam_cache;

static double ReadBeamEnergy(int run, Long64_t entry, const std::string &beam_dir, bool use_t1) {
	auto key = std::make_pair(run, use_t1);
	auto it = g_beam_cache.find(key);
	if (it == g_beam_cache.end()) {
		TString path = TString::Format("%s/beam_%s%04d.root", beam_dir.c_str(),
			use_t1 ? "t1_" : "", run);
		TFile *f = TFile::Open(path);
		if (!f || f->IsZombie()) {
			if (f) delete f;
			g_beam_cache[key] = {};
			return -1.0;
		}
		TTree *t = (TTree*)f->Get("tree");
		if (!t) {
			f->Close(); delete f;
			g_beam_cache[key] = {};
			return -1.0;
		}
		BeamFileCache cache;
		cache.file = f;
		cache.tree = t;
		cache.run = run;
		if (!t->GetBranch("E")) {
			f->Close(); delete f;
			g_beam_cache[key] = {};
			return -1.0;
		}
		g_beam_cache[key] = std::move(cache);
		it = g_beam_cache.find(key);
		it->second.tree->SetBranchAddress("E", &it->second.E_beam);
	}
	Long64_t max_entries = it->second.tree->GetEntries();
	if (!it->second.tree || entry < 0 || entry >= max_entries) return -1.0;
	it->second.tree->GetEntry(entry);
	return it->second.E_beam;
}

static std::map<int, double> g_fixed_beam_cache;

static double ReadFixedBeamEnergy(int run, const std::string &beam_dir) {
	auto it = g_fixed_beam_cache.find(run);
	if (it != g_fixed_beam_cache.end()) return it->second;

	TString path = TString::Format("%s/check_tof_and_beam_%04d.txt", beam_dir.c_str(), run);
	std::ifstream fin(path.Data());
	if (!fin.is_open()) {
		fprintf(stderr,
			"Error: per-run beam energy file not found: %s\n"
			"  Please run 'check_tof_and_beam -r %d' first to generate it.\n",
			path.Data(), run);
		exit(1);
	}

	std::string header;
	std::getline(fin, header);
	double E0 = -1.0, sigma_E = 0.0, t0 = 0.0, sigma_t0 = 0.0, dt = 0.0;
	fin >> E0 >> sigma_E >> t0 >> sigma_t0 >> dt;
	fin.close();

	if (E0 <= 0) {
		fprintf(stderr, "Error: invalid beam energy E0=%.4f in %s\n", E0, path.Data());
		exit(1);
	}

	g_fixed_beam_cache[run] = E0;
	return E0;
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
	if (ec.h_excitation_mm) delete ec.h_excitation_mm;
	if (ec.h_excitation_invmass) delete ec.h_excitation_invmass;
	if (ec.h_excitation_mm_fixed) delete ec.h_excitation_mm_fixed;
	int nbins = g_ctx.entry_excitation_bins
		? g_ctx.entry_excitation_bins->GetIntNumber() : 30;
	if (nbins < 5) nbins = 5;
	double bin_keV = 30000.0 / nbins;
	TString bin_label = TString::Format("Counts / %.0f keV", bin_keV);

	ec.h_excitation_mm = new TH1D(
		"h_excitation_mm",
		TString::Format("^{10}C Excitation (miss.mass, ^{14}O);E_{x} (MeV);%s",
			bin_label.Data()),
		nbins, -10, 20);
	ec.h_excitation_mm->SetDirectory(0);

	ec.h_excitation_invmass = new TH1D(
		"h_excitation_invmass",
		TString::Format("^{10}C Excitation (inv.mass, %s);E_{x} (MeV);%s",
			BeamLabel(), bin_label.Data()),
		nbins, -10, 20);
	ec.h_excitation_invmass->SetDirectory(0);

	ec.h_excitation_mm_fixed = new TH1D(
		"h_excitation_mm_fixed",
		TString::Format("^{10}C Excitation (miss.mass, per-run E_{b});E_{x} (MeV);%s",
			bin_label.Data()),
		nbins, -10, 20);
	ec.h_excitation_mm_fixed->SetDirectory(0);
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

	catima::Material cd2_mat;
	cd2_mat.add_element(12.0, 6.0, 1.0);
	cd2_mat.add_element(2.0, 1.0, 2.0);
	cd2_mat.density(0.79867);
	std::string cache_dir = brill::JoinPath(g_ctx.config.workspace, g_ctx.config.paths.energy_calculator);
	TString cache_10C = TString::Format("%s/target_z6_a10_cd2.root", cache_dir.c_str());
	TString cache_6Li = TString::Format("%s/target_z3_a6_cd2.root", cache_dir.c_str());
	if (!std::filesystem::exists(cache_6Li.Data())) {
		printf("  WARNING: 6Li cache %s not found. Energy-loss correction may fail.\n", cache_6Li.Data());
	}
	brill::RangeEnergyCalculator calc_10C(6, 10, cd2_mat, cache_10C.Data());
	brill::RangeEnergyCalculator calc_6Li(3, 6, cd2_mat, cache_6Li.Data());
	const double half_thick_um = 122.0 / 2.0;

	static std::unique_ptr<brill::RangeEnergyCalculator> s_beam_calc;
	static bool s_beam_calc_checked = false;
	static bool s_beam_calc_ok = false;
	TString beam_cache_path = TString::Format("%s/target_z8_a14_cd2.root", cache_dir.c_str());
	if (!s_beam_calc_checked) {
		s_beam_calc_checked = true;
		if (std::filesystem::exists(beam_cache_path.Data())) {
			s_beam_calc.reset(new brill::RangeEnergyCalculator(8, 14, cd2_mat, beam_cache_path.Data()));
			s_beam_calc_ok = true;
		} else {
			printf("  WARNING: beam cache %s not found, skipping beam-dependent methods.\n",
				beam_cache_path.Data());
			printf("  Run GUI_10C_4He first to generate it.\n");
			s_beam_calc_ok = false;
		}
	}

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
		fflush(stdout);
		passed++;

		double T_10C = ev.e1_10C + ev.e2_10C + ev.e3_10C + ev.e4_10C;
		double T_6Li = ev.e1_6Li + ev.e2_6Li;

		double z = g_ctx.entry_z ? g_ctx.entry_z->GetNumber() : 0.0;
		double exit_thick = std::max(0.0, half_thick_um - z);

		double costh_10C = std::cos(ev.theta_10C * M_PI / 180.0);
		double costh_6Li = std::cos(ev.theta_6Li * M_PI / 180.0);
		double exit_10C = exit_thick / std::max(costh_10C, 0.001);
		double exit_6Li = exit_thick / std::max(costh_6Li, 0.001);

		double r_10C = calc_10C.Range(T_10C);
		double r_10C_rxn = r_10C + exit_10C;
		double T_10C_corr = calc_10C.Energy(r_10C_rxn);

		double r_6Li = calc_6Li.Range(T_6Li);
		double r_6Li_rxn = r_6Li + exit_6Li;
		double T_6Li_corr = calc_6Li.Energy(r_6Li_rxn);

		g_ctx.ac.h_E_6Li_E_10C->Fill(T_10C_corr, T_6Li_corr);
		g_ctx.ac.h_10C_e_theta->Fill(ev.theta_10C, T_10C_corr);
		g_ctx.ac.h_6Li_e_theta->Fill(ev.theta_6Li, T_6Li_corr);
		g_ctx.ac.h_theta_theta->Fill(ev.theta_10C, ev.theta_6Li);

		// construct lab-frame 3D direction from t0d2 positions
		double dx_10C = ev.t0d2_10C_x - ev.target_x;
		double dy_10C = ev.t0d2_10C_y - ev.target_y;
		double dz_10C = ev.t0d2_10C_z;
		double norm_10C = std::sqrt(dx_10C*dx_10C + dy_10C*dy_10C + dz_10C*dz_10C);

		double dx_6Li = ev.t0d2_6Li_x - ev.target_x;
		double dy_6Li = ev.t0d2_6Li_y - ev.target_y;
		double dz_6Li = ev.t0d2_6Li_z;
		double norm_6Li = std::sqrt(dx_6Li*dx_6Li + dy_6Li*dy_6Li + dz_6Li*dz_6Li);

		double E_10C_tot = T_10C_corr + m_10C;
		double E_6Li_tot = T_6Li_corr + m_6Li;
		double p_10C = std::sqrt(E_10C_tot * E_10C_tot - m_10C * m_10C);
		double p_6Li = std::sqrt(E_6Li_tot * E_6Li_tot - m_6Li * m_6Li);

		TLorentzVector P_10C(p_10C * dx_10C / norm_10C,
			p_10C * dy_10C / norm_10C,
			p_10C * dz_10C / norm_10C, E_10C_tot);
		TLorentzVector P_6Li(p_6Li * dx_6Li / norm_6Li,
			p_6Li * dy_6Li / norm_6Li,
			p_6Li * dz_6Li / norm_6Li, E_6Li_tot);

		// === Method 2: momentum conservation + invariant mass (always) ===
		TVector3 p_14O_reco = P_10C.Vect() + P_6Li.Vect();
		double E_14O_reco = std::sqrt(p_14O_reco.Mag2() + m_14O * m_14O);
		TLorentzVector P_14O_reco(p_14O_reco, E_14O_reco);
		TLorentzVector P_d(0, 0, 0, m_d);
		TLorentzVector P_10C_star_inv = P_14O_reco + P_d - P_6Li;
		double E_x_inv = P_10C_star_inv.M() - m_10C;
		g_ctx.ec.h_excitation_invmass->Fill(E_x_inv);

		// === Method 1: missing mass with TOF beam energy (14O only) ===
		if (g_beam_filter == 1 && s_beam_calc_ok) {
			double E_beam_raw = ReadBeamEnergy(ev.run_number, ev.entry, g_ctx.beam_dir, g_ctx.is_t1);
			if (E_beam_raw > 0) {
				double norm_beam = std::sqrt(1.0 + ev.dir_x * ev.dir_x + ev.dir_y * ev.dir_y);
				double beam_depth = (half_thick_um + z) * norm_beam;
				double beam_range_before = s_beam_calc->Range(E_beam_raw);
				double beam_range_rxn = beam_range_before - beam_depth;
				if (beam_range_rxn > 0) {
					double E_beam = s_beam_calc->Energy(beam_range_rxn);
					double E_A = E_beam + m_14O;
					double p_A = std::sqrt(E_A * E_A - m_14O * m_14O);
					double beam_norm = std::sqrt(1.0 + ev.dir_x*ev.dir_x + ev.dir_y*ev.dir_y);
					TLorentzVector P_A(p_A * ev.dir_x / beam_norm,
						p_A * ev.dir_y / beam_norm,
						p_A * 1.0 / beam_norm, E_A);
					TLorentzVector P_B(0, 0, 0, m_d);
					TLorentzVector P_10C_star = P_A + P_B - P_6Li;
					double E_x_mm = P_10C_star.M() - m_10C;
					g_ctx.ec.h_excitation_mm->Fill(E_x_mm);
				}
			}
		}

		// === Method 3: missing mass with per-run fixed beam energy (all beams) ===
		if (s_beam_calc_ok) {
			double E_beam_fixed = ReadFixedBeamEnergy(ev.run_number, g_ctx.beam_dir);
			double norm_beam = std::sqrt(1.0 + ev.dir_x * ev.dir_x + ev.dir_y * ev.dir_y);
			double beam_depth = (half_thick_um + z) * norm_beam;
			double beam_range_before = s_beam_calc->Range(E_beam_fixed);
			double beam_range_rxn = beam_range_before - beam_depth;
			if (beam_range_rxn > 0) {
				double E_beam_at_rxn = s_beam_calc->Energy(beam_range_rxn);
				double E_A = E_beam_at_rxn + m_14O;
				double p_A = std::sqrt(E_A * E_A - m_14O * m_14O);
				double beam_norm = std::sqrt(1.0 + ev.dir_x*ev.dir_x + ev.dir_y*ev.dir_y);
				TLorentzVector P_A(p_A * ev.dir_x / beam_norm,
					p_A * ev.dir_y / beam_norm,
					p_A * 1.0 / beam_norm, E_A);
				TLorentzVector P_B(0, 0, 0, m_d);
				TLorentzVector P_10C_star = P_A + P_B - P_6Li;
				double E_x_mm_fixed = P_10C_star.M() - m_10C;
				g_ctx.ec.h_excitation_mm_fixed->Fill(E_x_mm_fixed);
			}
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
	ec.h_excitation_mm->Draw();
	ec.canvas->cd(2);
	ec.h_excitation_invmass->Draw();
	ec.canvas->cd(3);
	ec.h_excitation_mm_fixed->Draw();
	ec.canvas->Modified();
	ec.canvas->Update();
}

static void OnBeamChanged() {
	if (g_ctx.current_file.empty()) return;
	FillHistograms();
	DrawHistograms();

	std::string cut_path = "src/brill/Cut/cal_d1_d2_Li_large_cut.C";
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
	g_ctx.is_t1 = (filepath.find("_t1_") != std::string::npos);
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

	std::string cut_path = "src/brill/Cut/cal_d1_d2_Li_large_cut.C";
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
	g_ctx.beam_dir = brill::JoinPath(config.workspace, config.paths.beam);
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

	TGLabel *z_label = new TGLabel(beam_frame, "  z(µm): ");
	beam_frame->AddFrame(z_label, new TGLayoutHints(kLHintsCenterY, 10, 2, 2, 2));
	g_ctx.entry_z = new TGNumberEntry(beam_frame, 0, 5, -1,
		TGNumberFormat::kNESReal,
		TGNumberFormat::kNEAAnyNumber,
		TGNumberFormat::kNELLimitMinMax, -61, 61);
	beam_frame->AddFrame(g_ctx.entry_z, new TGLayoutHints(kLHintsCenterY, 2, 2, 2, 2));

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
				OnBeamChanged();
			}
		}
		g_menu_action = 0;
	}

	return 0;
}