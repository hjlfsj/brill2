// check_match - Interactive DSSD matching checker for debugging
// Usage: ./check_match -r <run> -n <start_event> <detector>
//        Press <Enter> for next event, 'q' + <Enter> to quit
//        Or use GUI: enter entry number, click "Read"

#include "include/t0/dssd.h"

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

#include <TGClient.h>
#include <TGFrame.h>
#include <TGLabel.h>
#include <TGNumberEntry.h>
#include <TGButton.h>
#include <TGLayout.h>

#include <filesystem>
#include <iostream>
#include <set>
#include <string>
#include <vector>
#include <cstdio>
#include <sys/select.h>

#include "external/cxxopts.hpp"
#include "include/config.h"
#include "include/event/ingot/dssd_event.h"
#include "include/event/t0/dssd_match_event.h"
#include "include/utils.h"

void PrintUsage(const cxxopts::Options &options) {
	std::cout << options.help() << "\n";
}

struct CheckContext {
	TCanvas *canvas = nullptr;
	TPad *pad_input = nullptr;
	TPad *pad_output = nullptr;
	TPad *pad_ininfo = nullptr;
	TPad *pad_outinfo = nullptr;
	TH2D *h_input = nullptr;
	TH2D *h_output = nullptr;
	TPaveText *info_input = nullptr;
	TPaveText *info_output = nullptr;
	const brill::SquareDetectorConfig *detector = nullptr;
	TTree *ipt = nullptr;
	brill::DssdEvent *raw_event = nullptr;
	brill::DssdEvent *normalized_event = nullptr;
	brill::DssdMatchEvent *match_event = nullptr;
	brill::DssdNormalizeParameters *parameters = nullptr;
	brill::SquareDetectorConfig *working_detector = nullptr;
	TGNumberEntry *entry_entry = nullptr;
	TGNumberEntry *entry_min_front = nullptr;
	TGNumberEntry *entry_min_back = nullptr;
	long long *current_entry = nullptr;
	long long total = 0;
	const char **merge_tag_names = nullptr;
};

static CheckContext g_ctx;

void UpdateCanvas() {
	const brill::DssdEvent &norm = *g_ctx.normalized_event;
	const brill::DssdMatchEvent &match = *g_ctx.match_event;
	long long current_entry = *g_ctx.current_entry;
	long long total = g_ctx.total;
	const char **merge_tag_names = g_ctx.merge_tag_names;

	// Top-left: Input 2D (cross pattern)
	g_ctx.pad_input->cd();
	g_ctx.h_input->Reset();
	g_ctx.h_input->SetTitle(
		TString::Format("Event %lld/%lld  Input: %df %db",
			current_entry, total, norm.front_num, norm.back_num)
	);
	for (int i = 0; i < norm.front_num; ++i) {
		int fs = norm.front_strip[i];
		double fe = norm.front_energy[i];
		for (int bs = 0; bs < g_ctx.detector->back_strips; ++bs) {
			g_ctx.h_input->Fill(fs, bs, fe);
		}
	}
	for (int i = 0; i < norm.back_num; ++i) {
		int bs = norm.back_strip[i];
		double be = norm.back_energy[i];
		for (int fs = 0; fs < g_ctx.detector->front_strips; ++fs) {
			g_ctx.h_input->Fill(fs, bs, be);
		}
	}
	g_ctx.h_input->Draw("colz");

	// Bottom-left: Output 2D (matched pairs)
	g_ctx.pad_output->cd();
	g_ctx.h_output->Reset();
	g_ctx.h_output->SetTitle(
		TString::Format("Output: %d matched clusters", match.num)
	);
	for (int i = 0; i < match.num; ++i) {
		g_ctx.h_output->Fill(match.front_strip[i], match.back_strip[i], match.energy[i]);
	}
	g_ctx.h_output->Draw("colz");

	// Top-right: Input event details
	g_ctx.pad_ininfo->cd();
	g_ctx.info_input->Clear();
	int min_front = g_ctx.entry_min_front->GetIntNumber();
	int min_back = g_ctx.entry_min_back->GetIntNumber();
	if (min_front > 0 || min_back > 0) {
		g_ctx.info_input->AddText(TString::Format("INPUT EVENT %lld/%lld  Filter: f>=%d b>=%d",
			current_entry, total, min_front, min_back));
	} else {
		g_ctx.info_input->AddText(TString::Format("INPUT EVENT %lld/%lld", current_entry, total));
	}
	g_ctx.info_input->AddText(" ");
	g_ctx.info_input->AddText(TString::Format("Front hits: %d", norm.front_num));
	for (int i = 0; i < norm.front_num; ++i) {
		g_ctx.info_input->AddText(TString::Format("  strip=%2d  E=%.1f  t=%.1f",
			norm.front_strip[i], norm.front_energy[i], norm.front_time[i]));
	}
	g_ctx.info_input->AddText(" ");
	g_ctx.info_input->AddText(TString::Format("Back hits: %d", norm.back_num));
	for (int i = 0; i < norm.back_num; ++i) {
		g_ctx.info_input->AddText(TString::Format("  strip=%2d  E=%.1f  t=%.1f",
			norm.back_strip[i], norm.back_energy[i], norm.back_time[i]));
	}
	g_ctx.info_input->Draw();

	// Bottom-right: Output event details
	g_ctx.pad_outinfo->cd();
	g_ctx.info_output->Clear();
	g_ctx.info_output->AddText(TString::Format("OUTPUT MATCHES: %d clusters", match.num));
	g_ctx.info_output->AddText(" ");
	if (match.num > 0) {
		for (int i = 0; i < match.num; ++i) {
			g_ctx.info_output->AddText(TString::Format("[%d] %s", i,
				(match.merge_tag[i] >= 0 && match.merge_tag[i] <= 3)
					? merge_tag_names[match.merge_tag[i]] : "?"));
			g_ctx.info_output->AddText(TString::Format("  x=%.1f  y=%.1f  E=%.1f  |diff|=%.1f",
				match.front_strip[i], match.back_strip[i],
				match.energy[i], match.energy_diff[i]));
			g_ctx.info_output->AddText(" ");
		}
	} else {
		g_ctx.info_output->AddText("NO MATCH");
	}
	g_ctx.info_output->Draw();

	g_ctx.canvas->Modified();
	g_ctx.canvas->Update();
}

void PrintEvent() {
	const brill::DssdEvent &norm = *g_ctx.normalized_event;
	const brill::DssdMatchEvent &match = *g_ctx.match_event;
	long long current_entry = *g_ctx.current_entry;
	long long total = g_ctx.total;
	const char **merge_tag_names = g_ctx.merge_tag_names;
	int min_front = g_ctx.entry_min_front->GetIntNumber();
	int min_back = g_ctx.entry_min_back->GetIntNumber();

	printf("\n=== Event %lld/%lld", current_entry, total);
	if (min_front > 0 || min_back > 0) {
		printf("  Filter: front >= %d, back >= %d", min_front, min_back);
	}
	printf(" ===\n");
	printf("Input: %d front hits, %d back hits\n", norm.front_num, norm.back_num);
	printf("Front: ");
	for (int i = 0; i < norm.front_num; ++i) {
		printf("%d(%.1f) ", norm.front_strip[i], norm.front_energy[i]);
	}
	printf("\nBack:  ");
	for (int i = 0; i < norm.back_num; ++i) {
		printf("%d(%.1f) ", norm.back_strip[i], norm.back_energy[i]);
	}
	printf("\nOutput: %d matched clusters\n", match.num);
	for (int i = 0; i < match.num; ++i) {
		printf("  [%d] x=%.1f y=%.1f E=%.1f |diff|=%.1f %s\n",
			i, match.front_strip[i], match.back_strip[i],
			match.energy[i], match.energy_diff[i],
			(match.merge_tag[i] >= 0 && match.merge_tag[i] <= 3)
				? merge_tag_names[match.merge_tag[i]] : "?");
	}
}

bool PassFilter() {
	int min_front = g_ctx.entry_min_front->GetIntNumber();
	int min_back = g_ctx.entry_min_back->GetIntNumber();
	if (min_front <= 0 && min_back <= 0) return true;
	return (g_ctx.normalized_event->front_num >= min_front &&
		g_ctx.normalized_event->back_num >= min_back);
}

void ProcessEvent() {
	long long &current_entry = *g_ctx.current_entry;
	if (current_entry >= g_ctx.total) {
		printf("\nReached end of file.\n");
		return;
	}
	while (current_entry < g_ctx.total) {
		g_ctx.ipt->GetEntry(current_entry);
		brill::ApplyDssdNormalize(*g_ctx.raw_event, *g_ctx.parameters, *g_ctx.normalized_event);
		if (PassFilter()) break;
		current_entry++;
	}
	if (current_entry >= g_ctx.total) {
		printf("\nNo more events matching filter.\n");
		return;
	}
	brill::MatchDssdEvent(*g_ctx.normalized_event, *g_ctx.working_detector,
		*g_ctx.match_event, nullptr);
	UpdateCanvas();
	PrintEvent();
}

int main(int argc, char **argv) {
	cxxopts::Options options("check_match", "Interactive DSSD matching checker.");
	options.add_options()
		("h,help", "Print help information.")
		("r,run", "Run number.", cxxopts::value<int>(), "run")
		("n,start", "Starting event number (default: 0)",
			cxxopts::value<long long>()->default_value("0"), "event")
		("w,window", "Override matching tolerance.", cxxopts::value<double>(), "window")
		("t,trigger", "Trigger type.", cxxopts::value<std::string>(), "trigger")
		("c,config", "Config file path.",
			cxxopts::value<std::string>()->default_value("config.toml"), "file")
		("detector", "Detector to check.", cxxopts::value<std::string>(), "detector");
	options.parse_positional({"detector"});

	auto result = options.parse(argc, argv);
	if (result.count("help")) { PrintUsage(options); return 0; }
	if (!result.count("run")) {
		std::cerr << "Error: Missing required option --run.\n";
		PrintUsage(options); return 1;
	}
	if (!result.count("detector")) {
		std::cerr << "Error: Missing detector positional argument.\n";
		PrintUsage(options); return 1;
	}

	const std::set<std::string> allowed = {"t0d1", "t0d2", "t0d3", "t0d4"};
	std::string detector_name = result["detector"].as<std::string>();
	if (allowed.count(detector_name) == 0) {
		std::cerr << "Error: Unsupported detector " << detector_name << ".\n";
		return 1;
	}

	brill::AppConfig config;
	if (brill::LoadConfig(result["config"].as<std::string>(), config)) return 1;
	if (result.count("trigger")) config.trigger = result["trigger"].as<std::string>();
	const int run = result["run"].as<int>();
	if (brill::IsJumpRun(config, run)) {
		std::cout << "Skipping jump run " << run << ".\n"; return 0;
	}
	long long start_event = result["start"].as<long long>();

	int normalize_file_run = 0;
	for (const auto &[s, u] : config.normalize.runs) {
		if (run >= s) normalize_file_run = u;
	}

	const brill::SquareDetectorConfig *detector =
		brill::FindDetectorConfig(config, detector_name);
	if (!detector) {
		std::cerr << "Error: Detector " << detector_name << " not found.\n";
		return 1;
	}
	double match_tolerance = detector->match_tolerance;
	if (result.count("window")) match_tolerance = result["window"].as<double>();
	brill::SquareDetectorConfig working_detector = *detector;
	working_detector.match_tolerance = match_tolerance;

	brill::DssdNormalizeParameters parameters;
	parameters.front_strips = detector->front_strips;
	parameters.back_strips = detector->back_strips;
	std::string normalize_dir = brill::JoinPath(config.workspace, config.paths.normalize);
	TString fpath = TString::Format("%s/%s_front_t1_%04d.txt",
		normalize_dir.c_str(), detector_name.c_str(), normalize_file_run);
	TString bpath = TString::Format("%s/%s_back_t1_%04d.txt",
		normalize_dir.c_str(), detector_name.c_str(), normalize_file_run);
	if (brill::ReadDssdNormalizeParameters(fpath.Data(), bpath.Data(), parameters)) return 1;

	TString ipath = TString::Format("%s/%s_%s%04d.root",
		brill::JoinPath(config.workspace, config.paths.ingot).c_str(),
		detector_name.c_str(), brill::TriggerInfix(config.trigger).c_str(), run);
	TFile ipf(ipath, "read");
	TTree *ipt = (TTree*)ipf.Get("tree");
	if (!ipt) { std::cerr << "Error: Get tree failed.\n"; return 1; }

	brill::DssdEvent raw_event;
	brill::DssdEvent normalized_event;
	brill::DssdMatchEvent match_event;
	brill::SetupInput(ipt, raw_event);

	long long total = ipt->GetEntries();
	printf("Total events in file: %lld\n", total);
	if (start_event >= total) {
		std::cerr << "Error: start_event " << start_event << " >= " << total << "\n";
		return 1;
	}

	TApplication app("check_match", nullptr, nullptr);

	gInterpreter->ProcessLine("bool g_check_match_read = false;");
	gInterpreter->ProcessLine("bool g_check_match_next = false;");

	TString canvas_title = TString::Format(
		"check_match - Run %d, Trigger %s, Detector %s, Tol %.0f",
		run, config.trigger.c_str(), detector_name.c_str(), match_tolerance);

	TGMainFrame *main_frame = new TGMainFrame(gClient->GetRoot(), 1400, 950);
	main_frame->SetWindowName(canvas_title.Data());

	TGHorizontalFrame *ctrl_frame = new TGHorizontalFrame(main_frame, 1400, 40);
	ctrl_frame->AddFrame(new TGLabel(ctrl_frame, "Entry:"),
		new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 5, 2, 5, 5));

	TGNumberEntry *entry_entry = new TGNumberEntry(ctrl_frame, start_event, 8, -1,
		TGNumberFormat::kNESInteger,
		TGNumberFormat::kNEANonNegative,
		TGNumberFormat::kNELLimitMinMax, 0, total - 1);
	ctrl_frame->AddFrame(entry_entry,
		new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 5, 5, 5));

	ctrl_frame->AddFrame(new TGLabel(ctrl_frame, "  f≥:"),
		new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 8, 2, 5, 5));

	TGNumberEntry *entry_min_front = new TGNumberEntry(ctrl_frame, 0, 2, -1,
		TGNumberFormat::kNESInteger,
		TGNumberFormat::kNEANonNegative,
		TGNumberFormat::kNELLimitMinMax, 0, 32);
	ctrl_frame->AddFrame(entry_min_front,
		new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 5, 5, 5));

	ctrl_frame->AddFrame(new TGLabel(ctrl_frame, "  b≥:"),
		new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 5, 2, 5, 5));

	TGNumberEntry *entry_min_back = new TGNumberEntry(ctrl_frame, 0, 2, -1,
		TGNumberFormat::kNESInteger,
		TGNumberFormat::kNEANonNegative,
		TGNumberFormat::kNELLimitMinMax, 0, 32);
	ctrl_frame->AddFrame(entry_min_back,
		new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 8, 5, 5));

	TGTextButton *read_btn = new TGTextButton(ctrl_frame, "Read");
	read_btn->SetCommand("g_check_match_read = true;");
	ctrl_frame->AddFrame(read_btn,
		new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 8, 5, 5, 5));

	TGTextButton *next_btn = new TGTextButton(ctrl_frame, "Next");
	next_btn->SetCommand("g_check_match_next = true;");
	ctrl_frame->AddFrame(next_btn,
		new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 5, 5, 5, 5));

	ctrl_frame->AddFrame(new TGLabel(ctrl_frame, "  q=quit"),
		new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 15, 5, 5, 5));

	main_frame->AddFrame(ctrl_frame,
		new TGLayoutHints(kLHintsTop | kLHintsExpandX, 2, 2, 2, 2));

	TRootEmbeddedCanvas *embed = new TRootEmbeddedCanvas("embed", main_frame, 1400, 880);
	main_frame->AddFrame(embed,
		new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));

	main_frame->MapSubwindows();
	main_frame->Resize(main_frame->GetDefaultSize());
	main_frame->MapWindow();

	TCanvas *canvas = embed->GetCanvas();
	canvas->cd();

	// Create 4 pads once: left-top, left-bottom, right-top, right-bottom
	TPad *pad_input = new TPad("pad_input", "Input", 0.0, 0.45, 0.55, 0.99);
	pad_input->Draw();
	TPad *pad_output = new TPad("pad_output", "Output", 0.0, 0.0, 0.55, 0.45);
	pad_output->Draw();
	TPad *pad_ininfo = new TPad("pad_ininfo", "InputInfo", 0.55, 0.45, 1.0, 0.99);
	pad_ininfo->Draw();
	TPad *pad_outinfo = new TPad("pad_outinfo", "OutputInfo", 0.55, 0.0, 1.0, 0.45);
	pad_outinfo->Draw();

	TH2D *h_input = new TH2D("h_input",
		"Input (cross);Front Strip;Back Strip",
		detector->front_strips, -0.5, detector->front_strips - 0.5,
		detector->back_strips, -0.5, detector->back_strips - 0.5);
	h_input->SetStats(0);

	TH2D *h_output = new TH2D("h_output",
		"Output (matched);Front Strip;Back Strip",
		detector->front_strips, -0.5, detector->front_strips - 0.5,
		detector->back_strips, -0.5, detector->back_strips - 0.5);
	h_output->SetStats(0);

	TPaveText *info_input = new TPaveText(0.02, 0.02, 0.98, 0.98, "brNDC");
	info_input->SetFillColor(kWhite);
	info_input->SetBorderSize(1);

	TPaveText *info_output = new TPaveText(0.02, 0.02, 0.98, 0.98, "brNDC");
	info_output->SetFillColor(kWhite);
	info_output->SetBorderSize(1);

	const char *merge_tag_names[] = {"1f-1b", "1f-2b", "2f-1b", "2f-2b"};

	long long current_entry = start_event;

	g_ctx.canvas = canvas;
	g_ctx.pad_input = pad_input;
	g_ctx.pad_output = pad_output;
	g_ctx.pad_ininfo = pad_ininfo;
	g_ctx.pad_outinfo = pad_outinfo;
	g_ctx.h_input = h_input;
	g_ctx.h_output = h_output;
	g_ctx.info_input = info_input;
	g_ctx.info_output = info_output;
	g_ctx.detector = detector;
	g_ctx.ipt = ipt;
	g_ctx.raw_event = &raw_event;
	g_ctx.normalized_event = &normalized_event;
	g_ctx.match_event = &match_event;
	g_ctx.parameters = &parameters;
	g_ctx.working_detector = &working_detector;
	g_ctx.entry_entry = entry_entry;
	g_ctx.entry_min_front = entry_min_front;
	g_ctx.entry_min_back = entry_min_back;
	g_ctx.current_entry = &current_entry;
	g_ctx.total = total;
	g_ctx.merge_tag_names = merge_tag_names;

	ProcessEvent();

	printf("\nControls: <Enter> = next matching, 'q' = quit, 'n<N>' = goto N, or use GUI\n");
	printf("Press <Enter> for next event, 'q' + <Enter> to quit: ");
	fflush(stdout);

	bool should_exit = false;

	while (!should_exit) {
		gSystem->ProcessEvents();

		Longptr_t val_read = gInterpreter->Calc("g_check_match_read");
		Longptr_t val_next = gInterpreter->Calc("g_check_match_next");

		if (val_read != 0) {
			gInterpreter->ProcessLine("g_check_match_read = false;");
			long long requested = entry_entry->GetIntNumber();
			if (requested >= 0 && requested < total) {
				current_entry = requested;
			}
			ProcessEvent();
			if (current_entry >= total) {
				should_exit = true;
			} else {
				entry_entry->SetIntNumber(current_entry);
				printf("\nPress <Enter> for next event, 'q' + <Enter> to quit: ");
				fflush(stdout);
			}
		}

		if (val_next != 0) {
			gInterpreter->ProcessLine("g_check_match_next = false;");
			current_entry++;
			ProcessEvent();
			if (current_entry >= total) {
				should_exit = true;
			} else {
				entry_entry->SetIntNumber(current_entry);
				printf("\nPress <Enter> for next event, 'q' + <Enter> to quit: ");
				fflush(stdout);
			}
		}

		fd_set fds;
		FD_ZERO(&fds);
		FD_SET(STDIN_FILENO, &fds);
		struct timeval tv = {0, 100000};
		int ret = select(STDIN_FILENO + 1, &fds, nullptr, nullptr, &tv);

		if (ret > 0 && FD_ISSET(STDIN_FILENO, &fds)) {
			char line[256];
			if (!fgets(line, sizeof(line), stdin)) {
				should_exit = true;
			} else if (line[0] == 'q' || line[0] == 'Q') {
				should_exit = true;
			} else if (line[0] == 'n' || line[0] == 'N') {
				long long req = atoll(line + 1);
				if (req >= 0 && req < total) {
					current_entry = req;
					ProcessEvent();
					if (current_entry >= total) {
						should_exit = true;
					} else {
						entry_entry->SetIntNumber(current_entry);
						printf("\nPress <Enter> for next event, 'q' + <Enter> to quit: ");
						fflush(stdout);
					}
				} else {
					printf("Invalid entry %lld. Range: 0-%lld\n", req, total - 1);
					printf("\nPress <Enter> for next event, 'q' + <Enter> to quit: ");
					fflush(stdout);
				}
				continue;
			} else {
				current_entry++;
				ProcessEvent();
				if (current_entry >= total) {
					should_exit = true;
				} else {
					entry_entry->SetIntNumber(current_entry);
					printf("\nPress <Enter> for next event, 'q' + <Enter> to quit: ");
					fflush(stdout);
				}
			}
		}
	}

	ipf.Close();
	printf("Done.\n");
	gApplication->Terminate(0);
}