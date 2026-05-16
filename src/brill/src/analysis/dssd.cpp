// #include "include/analysis/dssd.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>

// #include "include/analysis/file_name.h"
// #include "include/analysis/geometry.h"

namespace brill {

struct StripAccumulator {
	double sum = 0.0;
	int count = 0;
};

void ResetMergeEvent(DssdMergeEvent &event) {
	event.num = 0;
}

double MeanOrDefault(const StripAccumulator &accumulator, double fallback) {
	return accumulator.count > 0 ? accumulator.sum / accumulator.count : fallback;
}


// int RunDssdNormalize(
// 	const AppConfig &config,
// 	const std::string &detector,
// 	const std::string &trigger,
// 	int run,
// 	int end_run
// ) {
// 	const SquareDetectorConfig *detector_config = FindDetectorConfig(config, detector);
// 	if (!detector_config) {
// 		std::cerr << "Error: Unknown detector " << detector << ".\n";
// 		return -1;
// 	}

// 	std::vector<StripAccumulator> front(detector_config->front_strips);
// 	std::vector<StripAccumulator> back(detector_config->back_strips);
// 	double total_reference = 0.0;
// 	int total_count = 0;
// 	TH2F before("before", "front-back before normalize", 1000, 0, 60000, 1000, 0, 60000);
// 	for (int current = run; current <= end_run; ++current) {
// 		std::string input_path = ForgeFileName(config, detector, trigger, current);
// 		TFile input(input_path.c_str(), "read");
// 		TTree *tree = dynamic_cast<TTree*>(input.Get("tree"));
// 		if (!tree) {
// 			std::cerr << "Error: Open DSSD forge file " << input_path << " failed.\n";
// 			return -1;
// 		}
// 		DssdEvent event;
// 		SetupInput(tree, event, "");
// 		for (long long entry = 0; entry < tree->GetEntriesFast(); ++entry) {
// 			tree->GetEntry(entry);
// 			if (event.front_num != 1 || event.back_num != 1) continue;
// 			int fs = event.front_strip[0];
// 			int bs = event.back_strip[0];
// 			if (fs < 0 || fs >= detector_config->front_strips) continue;
// 			if (bs < 0 || bs >= detector_config->back_strips) continue;
// 			before.Fill(event.back_energy[0], event.front_energy[0]);
// 			double mean = 0.5 * (event.front_energy[0] + event.back_energy[0]);
// 			front[fs].sum += event.front_energy[0];
// 			front[fs].count += 1;
// 			back[bs].sum += event.back_energy[0];
// 			back[bs].count += 1;
// 			total_reference += mean;
// 			++total_count;
// 		}
// 	}
// 	double global_mean = total_count > 0 ? total_reference / total_count : 1.0;

// 	DssdNormalizeParameters parameters;
// 	parameters.front_strips = detector_config->front_strips;
// 	parameters.back_strips = detector_config->back_strips;
// 	for (int i = 0; i < detector_config->front_strips; ++i) {
// 		parameters.front_offset[i] = 0.0;
// 		double local_mean = MeanOrDefault(front[i], global_mean);
// 		parameters.front_scale[i] = std::fabs(local_mean) > 1e-9 ? global_mean / local_mean : 1.0;
// 	}
// 	for (int i = 0; i < detector_config->back_strips; ++i) {
// 		parameters.back_offset[i] = 0.0;
// 		double local_mean = MeanOrDefault(back[i], global_mean);
// 		parameters.back_scale[i] = std::fabs(local_mean) > 1e-9 ? global_mean / local_mean : 1.0;
// 	}

// 	std::string param_path = NormalizeParamFileName(config, detector, trigger, run, end_run);
// 	TFile output(param_path.c_str(), "recreate");
// 	TTree params_tree("params", "dssd normalize parameters");
// 	params_tree.Branch("front_strips", &parameters.front_strips, "front_strips/I");
// 	params_tree.Branch("back_strips", &parameters.back_strips, "back_strips/I");
// 	params_tree.Branch("front_offset", parameters.front_offset, "front_offset[128]/D");
// 	params_tree.Branch("front_scale", parameters.front_scale, "front_scale[128]/D");
// 	params_tree.Branch("back_offset", parameters.back_offset, "back_offset[128]/D");
// 	params_tree.Branch("back_scale", parameters.back_scale, "back_scale[128]/D");
// 	params_tree.Fill();
// 	params_tree.Write();
// 	before.Write();
// 	output.Close();

// 	for (int current = run; current <= end_run; ++current) {
// 		std::string input_path = ForgeFileName(config, detector, trigger, current);
// 		std::string normalize_path = NormalizeFileName(config, detector, trigger, current);
// 		TFile input(input_path.c_str(), "read");
// 		TTree *tree = dynamic_cast<TTree*>(input.Get("tree"));
// 		if (!tree) return -1;
// 		DssdEvent raw;
// 		SetupInput(tree, raw, "");
// 		TFile normalized_file(normalize_path.c_str(), "recreate");
// 		TTree opt("tree", "normalized dssd");
// 		DssdNormalizedEvent normalized;
// 		SetupOutput(&opt, normalized);
// 		TH2F after("after", "front-back after normalize", 1000, 0, 60000, 1000, 0, 60000);
// 		for (long long entry = 0; entry < tree->GetEntriesFast(); ++entry) {
// 			tree->GetEntry(entry);
// 			ApplyDssdNormalize(raw, parameters, normalized);
// 			if (normalized.front_num == 1 && normalized.back_num == 1) {
// 				after.Fill(normalized.back_energy[0], normalized.front_energy[0]);
// 			}
// 			opt.Fill();
// 		}
// 		opt.Write();
// 		after.Write();
// 		normalized_file.Close();
// 	}
// 	return 0;
// }

int ReadDssdNormalizeParameters(const std::string &path, DssdNormalizeParameters &parameters) {
	TFile input(path.c_str(), "read");
	TTree *tree = dynamic_cast<TTree*>(input.Get("params"));
	if (!tree) {
		std::cerr << "Error: Open DSSD parameter file " << path << " failed.\n";
		return -1;
	}
	tree->SetBranchAddress("front_strips", &parameters.front_strips);
	tree->SetBranchAddress("back_strips", &parameters.back_strips);
	tree->SetBranchAddress("front_offset", parameters.front_offset);
	tree->SetBranchAddress("front_scale", parameters.front_scale);
	tree->SetBranchAddress("back_offset", parameters.back_offset);
	tree->SetBranchAddress("back_scale", parameters.back_scale);
	tree->GetEntry(0);
	return 0;
}

void ApplyDssdNormalize(const DssdEvent &input, const DssdNormalizeParameters &parameters, DssdNormalizedEvent &output) {
	output.front_num = input.front_num;
	output.back_num = input.back_num;
	for (int i = 0; i < input.front_num && i < kMaxDssdHit; ++i) {
		output.front_strip[i] = input.front_strip[i];
		output.front_energy[i] =
			input.front_energy[i] * parameters.front_scale[input.front_strip[i]]
			+ parameters.front_offset[input.front_strip[i]];
		output.front_time[i] = input.front_time[i];
	}
	for (int i = 0; i < input.back_num && i < kMaxDssdHit; ++i) {
		output.back_strip[i] = input.back_strip[i];
		output.back_energy[i] =
			input.back_energy[i] * parameters.back_scale[input.back_strip[i]]
			+ parameters.back_offset[input.back_strip[i]];
		output.back_time[i] = input.back_time[i];
	}
}

void MergeDssdEvent(const DssdNormalizedEvent &input, const SquareDetectorConfig &detector, DssdMergeEvent &output) {
	ResetMergeEvent(output);
	bool back_used[kMaxDssdHit] = {false};
	for (int i = 0; i < input.front_num && output.num < kMaxDssdHit; ++i) {
		double best_diff = detector.merge_tolerance;
		int best = -1;
		for (int j = 0; j < input.back_num; ++j) {
			if (back_used[j]) continue;
			double diff = std::fabs(input.front_energy[i] - input.back_energy[j]);
			if (diff < best_diff) {
				best_diff = diff;
				best = j;
			}
		}
		if (best < 0) continue;
		int index = output.num++;
		back_used[best] = true;
		output.front_strip[index] = input.front_strip[i];
		output.back_strip[index] = input.back_strip[best];
		output.front_energy[index] = input.front_energy[i];
		output.back_energy[index] = input.back_energy[best];
		output.energy[index] = 0.5 * (input.front_energy[i] + input.back_energy[best]);
		output.time[index] = input.front_time[i];
		DssdStripToPosition(
			detector,
			input.front_strip[i],
			input.back_strip[best],
			output.x[index],
			output.y[index],
			output.z[index]
		);
		output.merge_tag[index] = 0;
	}
}

int RunDssdMerge(
	const AppConfig &config,
	const std::string &detector,
	const std::string &trigger,
	int run
) {
	const SquareDetectorConfig *detector_config = FindDetectorConfig(config, detector);
	if (!detector_config) return -1;
	DssdNormalizeParameters parameters;
	if (ReadDssdNormalizeParameters(NormalizeParamFileName(config, detector, trigger, run, run), parameters)) {
		return -1;
	}
	std::string input_path = NormalizeFileName(config, detector, trigger, run);
	std::string output_path = MergeFileName(config, detector, trigger, run);
	TFile input(input_path.c_str(), "read");
	TTree *tree = dynamic_cast<TTree*>(input.Get("tree"));
	if (!tree) {
		std::cerr << "Error: Open normalized DSSD file failed.\n";
		return -1;
	}
	DssdNormalizedEvent normalized;
	SetupInput(tree, normalized, "");
	TFile output(output_path.c_str(), "recreate");
	TTree opt("tree", "merged dssd");
	DssdMergeEvent merged;
	SetupOutput(&opt, merged);
	TH2F pid("pid", "single detector front-back", 1000, 0, 60000, 1000, 0, 60000);
	for (long long entry = 0; entry < tree->GetEntriesFast(); ++entry) {
		tree->GetEntry(entry);
		MergeDssdEvent(normalized, *detector_config, merged);
		for (int i = 0; i < merged.num; ++i) {
			pid.Fill(merged.back_energy[i], merged.front_energy[i]);
		}
		opt.Fill();
	}
	opt.Write();
	pid.Write();
	output.Close();
	return 0;
}

void FillDssdPidHistogram(const DssdMergeEvent &left, const DssdMergeEvent &right, TH2F &hist) {
	if (left.num <= 0 || right.num <= 0) return;
	hist.Fill(right.energy[0], left.energy[0]);
}

} // namespace glimmer
