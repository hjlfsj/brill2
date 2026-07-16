add_executable(normalize normalize.cpp)
target_link_libraries(normalize PRIVATE config dssd_event t0_dssd ROOT::Graf)

add_executable(ppac_normalize ppac_normalize.cpp)
target_link_libraries(
	ppac_normalize
	PRIVATE
	config
	ppac_event
	ROOT::Spectrum
	ROOT::RIO
	ROOT::Tree
	ROOT::Hist
	ROOT::Graf
	ROOT::Gpad
)

add_executable(track_ppac track_ppac.cpp)
target_link_libraries(
	track_ppac
	PRIVATE
	config
	ppac_track
	ROOT::RIO
	ROOT::Tree
)

add_executable(match_dssd match_dssd.cpp)
target_link_libraries(
	match_dssd
	PRIVATE
	config
	t0_dssd
	ROOT::Hist
)

add_executable(check_match check_match.cpp)
target_link_libraries(
	check_match
	PRIVATE
	config
	t0_dssd
	ROOT::Hist
	ROOT::Graf
	ROOT::Gpad
	ROOT::Gui
)

add_executable(check_ppac_track check_ppac_track.cpp)
target_link_libraries(
	check_ppac_track
	PRIVATE
	config
	ppac_track
	ROOT::RIO
	ROOT::Tree
	ROOT::Hist
	ROOT::Graf
	ROOT::Gpad
	ROOT::Gui
)

add_executable(track_t0 track_t0.cpp)
target_link_libraries(
	track_t0
	PRIVATE
	config
	dssd_match_event
	t0_event
	silicon_event
	ROOT::Graf
)

add_executable(rebuild_t0 rebuild_t0.cpp)
target_link_libraries(
	rebuild_t0
	PRIVATE
	config
	t0_event
	particle_event
	energy_calculator
)

add_executable(rebuild_2alpha_2p rebuild_2alpha_2p.cpp)
target_link_libraries(
	rebuild_2alpha_2p
	PRIVATE
	config
	dssd_match_event
	energy_calculator
	rebuild_physics
	ROOT::Graf
)

add_executable(calibrate_t0 calibrate_t0.cpp)
target_link_libraries(
	calibrate_t0
	PRIVATE
	config
	t0_event
	energy_calculator
	ROOT::Minuit
)

add_subdirectory(estimate)
add_subdirectory(rebuild)