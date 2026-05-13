#pragma once

#include <string>

namespace glimmer {

struct ProgramOptions {
	std::string config_path = "config.toml";
	std::string detector;
	std::string trigger;
	int run = -1;
	int end_run = -1;
	bool help = false;
};

int ParseProgramOptions(
	int argc,
	char **argv,
	ProgramOptions &options,
	bool need_detector = false,
	bool allow_end_run = false
);

void PrintGeneralUsage(
	const char *name,
	bool need_detector = false,
	bool allow_end_run = false
);

} // namespace glimmer
