#include "include/analysis/cli.h"

#include <cstdlib>
#include <iostream>

namespace glimmer {

void PrintGeneralUsage(const char *name, bool need_detector, bool allow_end_run) {
	std::cout
		<< "Usage: " << name << " --run <run> [options]\n"
		<< "Options:\n"
		<< "  --config <path>          Config file path, default is config.toml.\n"
		<< "  -t, --trigger <trigger>  Trigger name.\n";
	if (allow_end_run) {
		std::cout << "  --end-run <run>          Last run for chained normalization.\n";
	}
	if (need_detector) {
		std::cout << "  --detector <name>        Detector name.\n";
	}
	std::cout << "  -h, --help               Print this help information.\n";
}

int ParseProgramOptions(
	int argc,
	char **argv,
	ProgramOptions &options,
	bool need_detector,
	bool allow_end_run
) {
	for (int i = 1; i < argc; ++i) {
		std::string arg = argv[i];
		if (arg == "-h" || arg == "--help") {
			options.help = true;
			return 0;
		}
		if ((arg == "--config" || arg == "-c") && i + 1 < argc) {
			options.config_path = argv[++i];
		} else if ((arg == "--run" || arg == "-r") && i + 1 < argc) {
			options.run = std::atoi(argv[++i]);
		} else if (arg == "--end-run" && i + 1 < argc) {
			options.end_run = std::atoi(argv[++i]);
		} else if ((arg == "--trigger" || arg == "-t") && i + 1 < argc) {
			options.trigger = argv[++i];
		} else if (arg == "--detector" && i + 1 < argc) {
			options.detector = argv[++i];
		} else {
			std::cerr << "Error: Unknown or incomplete option " << arg << ".\n";
			return -1;
		}
	}
	if (options.end_run < 0) options.end_run = options.run;
	if (options.run < 0 && !options.help) {
		std::cerr << "Error: Missing --run.\n";
		return -1;
	}
	if (need_detector && options.detector.empty() && !options.help) {
		std::cerr << "Error: Missing --detector.\n";
		return -1;
	}
	if (!allow_end_run) options.end_run = options.run;
	return 0;
}

} // namespace glimmer
