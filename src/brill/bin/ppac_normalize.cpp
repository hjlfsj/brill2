#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>
#include <TChain.h>
#include <TString.h>
#include <TF1.h>
#include <TGraph.h>
#include <TSpectrum.h>

#include "external/cxxopts.hpp"
#include "include/config.h"
#include "include/event/ingot/ppac_event.h"
#include "include/utils.h"

namespace brill {

struct PpacOffsetResult {
	double delay_x_peak = 0.0;
	double delay_x_sigma = 0.0;
	double delay_y_peak = 0.0;
	double delay_y_sigma = 0.0;
	double position_x_p0 = 0.0;
	double position_x_p1 = 0.0;
	double position_x_chi2 = 0.0;
	double position_x_ndf = 1.0;
	double position_y_p0 = 0.0;
	double position_y_p1 = 0.0;
	double position_y_chi2 = 0.0;
	double position_y_ndf = 1.0;
	double anode_efficiency = 0.0;
	double x_cathode_efficiency = 0.0;
	double y_cathode_efficiency = 0.0;
};

void PrintUsage(const cxxopts::Options &options) {
	std::cout << options.help() << "\n";
}

bool FitDelayPeak(TH1D *hist, TSpectrum *spectrum, double &peak, double &sigma, const char *name) {
	int nfound = spectrum->Search(hist, 0.7, "", 0.9);
	if (nfound < 1) {
		std::cerr << "Warning: No peak found in " << name << "\n";
		return false;
	}
	double *pos = spectrum->GetPositionX();
	hist->Fit("gaus", "Q", "", pos[0]-1.2, pos[0]+1.2);
	TF1 *fit = hist->GetFunction("gaus");
	if (!fit) {
		std::cerr << "Warning: Fit failed for " << name << "\n";
		return false;
	}
	peak = fit->GetParameter(1);
	sigma = fit->GetParameter(2);
	return true;
}

bool FitChipPeaks(
	TH1D *hist, TSpectrum *spectrum, 
	std::vector<double> &chip_peak_positions,
	double search_sigma, double search_threshold,
	const char *direction
) {
	int nfound = spectrum->Search(hist, search_sigma, "", search_threshold);
	if (nfound < 10) {
		std::cerr << "Warning: " << direction << " only " << nfound << " peaks found, need at least 10\n";
		return false;
	}
	
	double *peak_positions = new double[nfound];
	for (int i = 0; i < nfound; ++i) {
		peak_positions[i] = spectrum->GetPositionX()[i];
	}
	std::sort(peak_positions, peak_positions + nfound);
	
	int pair_idx = -1;
	double min_gap = 1e9;
	for (int i = 1; i < nfound; ++i) {
		double gap = peak_positions[i] - peak_positions[i-1];
		if (gap < min_gap) {
			min_gap = gap;
			pair_idx = i;
		}
	}
	
	int left_count = pair_idx - 1;
	int start_idx;
	if (left_count >= 10) {
		start_idx = pair_idx - 9;
	} else {
		start_idx = pair_idx + 1;
	}
	
	if (start_idx < 0 || start_idx + 8 > nfound) {
		std::cerr << "Warning: " << direction << " not enough peaks (start_idx=" << start_idx << ")\n";
		delete[] peak_positions;
		return false;
	}
	
	std::cout << "  " << direction << " wire pair at peak " << pair_idx << " (gap=" << min_gap 
	          << " ns), left_count=" << left_count << ", start_idx=" << start_idx << "\n";
	
	chip_peak_positions.clear();
	for (int i = 0; i < 8; ++i) {
		int idx = start_idx + i;
		TF1 *f = new TF1(TString::Format("f%d", i), "gaus");
		hist->Fit(f, "SQ+", "", peak_positions[idx]-0.6, peak_positions[idx]+0.6);
		chip_peak_positions.push_back(f->GetParameter(1));
		delete f;
	}
	
	delete[] peak_positions;
	return true;
}

bool FitPositionLinear(TGraph *gr, double &p0, double &p1, double &chi2, double &ndf) {
	gr->Fit("pol1", "Q");
	TF1 *fit = gr->GetFunction("pol1");
	if (!fit) return false;
	p0 = fit->GetParameter(0);
	p1 = fit->GetParameter(1);
	chi2 = fit->GetChisquare();
	ndf = fit->GetNDF();
	return true;
}

int CalibratePpacOffset(
	TChain &chain,
	PpacOffsetResult results[3]
) {
	PpacEvent event;
	SetupInput(&chain, event, "");
	
	long long total_entries = chain.GetEntries();
	
	for (int ppac_idx = 0; ppac_idx < 3; ++ppac_idx) {
		std::cout << "\nProcessing PPAC" << (ppac_idx + 1) << "...\n";
		
		TString prefix = TString::Format("ppac%d_", ppac_idx + 1);
		
		int ch_x1 = PpacChannelMap::ch_x1[ppac_idx];
		int ch_x2 = PpacChannelMap::ch_x2[ppac_idx];
		int ch_y1 = PpacChannelMap::ch_y1[ppac_idx];
		int ch_y2 = PpacChannelMap::ch_y2[ppac_idx];
		int ch_anode = PpacChannelMap::ch_anode[ppac_idx];
		
		TH1D *h_delay_x = new TH1D(prefix + "Tx_Delay", "Tx Delay", 6000, -200, 400);
		TH1D *h_delay_y = new TH1D(prefix + "Ty_Delay", "Ty Delay", 6000, -200, 400);
		
		TSpectrum *spectrum_delay_x = new TSpectrum(5);
		TSpectrum *spectrum_delay_y = new TSpectrum(5);
		
		long long progress_percentage = 0;
		
		long long total_events = 0;
		long long anode_valid_count = 0;
		long long x_cathode_valid_count = 0;
		long long y_cathode_valid_count = 0;
		
		for (long long entry = 0; entry < total_entries; ++entry) {
			if (entry * 100ll / total_entries > progress_percentage) {
				progress_percentage = entry * 100ll / total_entries;
				printf("\r  Pass1: %3lld%%", progress_percentage);
				fflush(stdout);
			}
			chain.GetEntry(entry);
			
			total_events++;
			
			bool anode_valid = event.valid[ch_anode] && event.time[ch_anode] > 0;
			bool x1_valid = event.valid[ch_x1] && event.time[ch_x1] > 0;
			bool x2_valid = event.valid[ch_x2] && event.time[ch_x2] > 0;
			bool y1_valid = event.valid[ch_y1] && event.time[ch_y1] > 0;
			bool y2_valid = event.valid[ch_y2] && event.time[ch_y2] > 0;
			
			if (anode_valid) anode_valid_count++;
			if (anode_valid && x1_valid && x2_valid) {
				x_cathode_valid_count++;
				double delay_x = event.time[ch_x1] + event.time[ch_x2] - 2 * event.time[ch_anode];
				h_delay_x->Fill(delay_x);
			}
			if (anode_valid && y1_valid && y2_valid) {
				y_cathode_valid_count++;
				double delay_y = event.time[ch_y1] + event.time[ch_y2] - 2 * event.time[ch_anode];
				h_delay_y->Fill(delay_y);
			}
		}
		printf("\r  Pass1: 100%%\n");
		
		results[ppac_idx].anode_efficiency = total_events > 0 ? double(anode_valid_count) / total_events * 100 : 0;
		results[ppac_idx].x_cathode_efficiency = anode_valid_count > 0 ? double(x_cathode_valid_count) / anode_valid_count * 100 : 0;
		results[ppac_idx].y_cathode_efficiency = anode_valid_count > 0 ? double(y_cathode_valid_count) / anode_valid_count * 100 : 0;
		
		FitDelayPeak(h_delay_x, spectrum_delay_x, results[ppac_idx].delay_x_peak, results[ppac_idx].delay_x_sigma, "Tx Delay");
		FitDelayPeak(h_delay_y, spectrum_delay_y, results[ppac_idx].delay_y_peak, results[ppac_idx].delay_y_sigma, "Ty Delay");
		
		std::cout << "  Tx_Delay: peak=" << results[ppac_idx].delay_x_peak << ", sigma=" << results[ppac_idx].delay_x_sigma << "\n";
		std::cout << "  Ty_Delay: peak=" << results[ppac_idx].delay_y_peak << ", sigma=" << results[ppac_idx].delay_y_sigma << "\n";
		
		h_delay_x->Write();
		h_delay_y->Write();
		
		TH1D *h_dtx = new TH1D(prefix + "dtx12", "dtx12", 4000, -200, 200);
		TH1D *h_dty = new TH1D(prefix + "dty21", "dty21", 4000, -200, 200);
		
		TSpectrum *spectrum_dtx = new TSpectrum(50);
		TSpectrum *spectrum_dty = new TSpectrum(50);
		
		progress_percentage = 0;
		for (long long entry = 0; entry < total_entries; ++entry) {
			if (entry * 100ll / total_entries > progress_percentage) {
				progress_percentage = entry * 100ll / total_entries;
				printf("\r  Pass2: %3lld%%", progress_percentage);
				fflush(stdout);
			}
			chain.GetEntry(entry);
			
			bool anode_valid = event.valid[ch_anode] && event.time[ch_anode] > 0;
			bool x1_valid = event.valid[ch_x1] && event.time[ch_x1] > 0;
			bool x2_valid = event.valid[ch_x2] && event.time[ch_x2] > 0;
			bool y1_valid = event.valid[ch_y1] && event.time[ch_y1] > 0;
			bool y2_valid = event.valid[ch_y2] && event.time[ch_y2] > 0;
			
			if (anode_valid && x1_valid && x2_valid) {
				double delay_x = event.time[ch_x1] + event.time[ch_x2] - 2 * event.time[ch_anode];
				if (std::fabs(delay_x - results[ppac_idx].delay_x_peak) < 3 * results[ppac_idx].delay_x_sigma) {
					double dtx = event.time[ch_x1] - event.time[ch_x2];
					h_dtx->Fill(dtx);
				}
			}
			if (anode_valid && y1_valid && y2_valid) {
				double delay_y = event.time[ch_y1] + event.time[ch_y2] - 2 * event.time[ch_anode];
				if (std::fabs(delay_y - results[ppac_idx].delay_y_peak) < 3 * results[ppac_idx].delay_y_sigma) {
					double dty = event.time[ch_y1] - event.time[ch_y2];
					h_dty->Fill(dty);
				}
			}
		}
		printf("\r  Pass2: 100%%\n");
		
		std::vector<double> chip_peaks_x, chip_peaks_y;
		FitChipPeaks(h_dtx, spectrum_dtx, chip_peaks_x, 2, 0.2, "X");
		FitChipPeaks(h_dty, spectrum_dty, chip_peaks_y, 3, 0.001, "Y");
		
		h_dtx->Write();
		h_dty->Write();
		
		if (chip_peaks_x.size() == 8) {
			double chip_positions[8] = {-4, -3, -2, -1, 0, 1, 2, 3};
			TGraph *graph_x = new TGraph(8, chip_peaks_x.data(), chip_positions);
			graph_x->SetName(prefix + "x_graph");
			graph_x->SetDrawOption("AP*");
			FitPositionLinear(graph_x, results[ppac_idx].position_x_p0, results[ppac_idx].position_x_p1, 
			                  results[ppac_idx].position_x_chi2, results[ppac_idx].position_x_ndf);
			std::cout << "  X fit: p0=" << results[ppac_idx].position_x_p0 << ", p1=" << results[ppac_idx].position_x_p1 
					  << ", chi2/ndf=" << (results[ppac_idx].position_x_ndf > 0 ? results[ppac_idx].position_x_chi2 / results[ppac_idx].position_x_ndf : 0) << "\n";
			graph_x->Write();
			delete graph_x;
		}
		
		if (chip_peaks_y.size() == 8) {
			double chip_positions[8] = {-4, -3, -2, -1, 0, 1, 2, 3};
			TGraph *graph_y = new TGraph(8, chip_peaks_y.data(), chip_positions);
			graph_y->SetName(prefix + "y_graph");
			graph_y->SetDrawOption("AP*");
			FitPositionLinear(graph_y, results[ppac_idx].position_y_p0, results[ppac_idx].position_y_p1, 
			                  results[ppac_idx].position_y_chi2, results[ppac_idx].position_y_ndf);
			std::cout << "  Y fit: p0=" << results[ppac_idx].position_y_p0 << ", p1=" << results[ppac_idx].position_y_p1 
					  << ", chi2/ndf=" << (results[ppac_idx].position_y_ndf > 0 ? results[ppac_idx].position_y_chi2 / results[ppac_idx].position_y_ndf : 0) << "\n";
			graph_y->Write();
			delete graph_y;
		}
		
		delete h_delay_x;
		delete h_delay_y;
		delete h_dtx;
		delete h_dty;
		delete spectrum_delay_x;
		delete spectrum_delay_y;
		delete spectrum_dtx;
		delete spectrum_dty;
	}
	
	return 0;
}

} // namespace brill

int main(int argc, char **argv) {
	cxxopts::Options options("ppac_normalize", "PPAC offset calibration.");
	options.add_options()
		("h,help", "Print help information.")
		("r,run", "Run number.", cxxopts::value<int>(), "run")
		("t,trigger", "Trigger type.", cxxopts::value<std::string>()->default_value("main"), "trigger")
		("c,config", "Config file path.", cxxopts::value<std::string>()->default_value("config.toml"), "file");
	
	auto result = options.parse(argc, argv);
	if (result.count("help")) {
		brill::PrintUsage(options);
		return 0;
	}
	if (!result.count("run")) {
		std::cerr << "Error: Missing required option --run.\n";
		brill::PrintUsage(options);
		return 1;
	}
	
	brill::AppConfig config;
	if (brill::LoadConfig(result["config"].as<std::string>(), config)) {
		return 1;
	}
	
	std::string trigger = result["trigger"].as<std::string>();
	config.trigger = trigger;
	
	const int run = result["run"].as<int>();
	
	std::string input_dir = brill::JoinPath(config.workspace, config.paths.ingot);
	std::string output_dir = brill::JoinPath(config.workspace, config.paths.normalize);
	std::filesystem::create_directories(output_dir);
	
	TChain chain("tree");
	TString file_path = TString::Format(
		"%s/ppac_%s%04d.root",
		input_dir.c_str(),
		brill::TriggerInfix(config.trigger).c_str(),
		run
	);
	if (std::filesystem::exists(file_path.Data())) {
		chain.Add(file_path);
		std::cout << "Add run " << run << ".\n";
	} else {
		std::cerr << "Error: Input file not found: " << file_path << "\n";
		return 1;
	}
	
	if (chain.GetEntries() == 0) {
		std::cerr << "Error: No input files found.\n";
		return 1;
	}
	
	std::string prefix = (trigger == "main") ? "" : "t1_";
	TString root_path = TString::Format(
		"%s/ppac_offset_%s%04d.root",
		output_dir.c_str(),
		prefix.c_str(),
		run
	);
	TFile root_file(root_path, "recreate");
	
	brill::PpacOffsetResult offset_results[3];
	brill::CalibratePpacOffset(chain, offset_results);
	
	TString txt_path = TString::Format(
		"%s/ppac_offset_%s%04d.txt",
		output_dir.c_str(),
		prefix.c_str(),
		run
	);
	std::ofstream txt_file(txt_path.Data());
	txt_file << "ppac_id\tdelay_x_peak\tdelay_x_sigma\tdelay_y_peak\tdelay_y_sigma\tposition_x_p0\tposition_x_p1\tposition_x_chi2/ndf\tposition_y_p0\tposition_y_p1\tposition_y_chi2/ndf\tanode_eff(%)\tx_cathode_eff(%)\ty_cathode_eff(%)\n";
	
	for (int i = 0; i < 3; ++i) {
		txt_file << i+1 << "\t"
				 << offset_results[i].delay_x_peak << "\t"
				 << offset_results[i].delay_x_sigma << "\t"
				 << offset_results[i].delay_y_peak << "\t"
				 << offset_results[i].delay_y_sigma << "\t"
				 << offset_results[i].position_x_p0 << "\t"
				 << offset_results[i].position_x_p1 << "\t"
				 << (offset_results[i].position_x_ndf > 0 ? offset_results[i].position_x_chi2 / offset_results[i].position_x_ndf : 0) << "\t"
				 << offset_results[i].position_y_p0 << "\t"
				 << offset_results[i].position_y_p1 << "\t"
				 << (offset_results[i].position_y_ndf > 0 ? offset_results[i].position_y_chi2 / offset_results[i].position_y_ndf : 0) << "\t"
				 << offset_results[i].anode_efficiency << "\t"
				 << offset_results[i].x_cathode_efficiency << "\t"
				 << offset_results[i].y_cathode_efficiency << "\n";
	}
	txt_file.close();
	
	std::cout << "\nResults saved to:\n";
	std::cout << "  " << root_path << "\n";
	std::cout << "  " << txt_path << "\n";
	
	return 0;
}