#include "include/analysis/cli.h"
#include "include/analysis/config.h"
#include "include/analysis/t0.h"

int main(int argc, char **argv) {
	glimmer::ProgramOptions options;
	if (glimmer::ParseProgramOptions(argc, argv, options) != 0 || options.help) {
		glimmer::PrintGeneralUsage(argv[0]);
		return options.help ? 0 : 1;
	}
	glimmer::AppConfig config;
	if (glimmer::LoadConfig(options.config_path, config) != 0) return 1;
	if (options.trigger.empty()) options.trigger = config.trigger_default;
	return glimmer::RunEstimatePpacCenter(config, options.trigger, options.run) == 0 ? 0 : 1;
}
