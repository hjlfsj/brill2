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

add_executable(GUI_track GUI_track.cpp)
target_link_libraries(
	GUI_track
	PRIVATE
	config
	ppac_track
	dssd_match_event
	ROOT::RIO
	ROOT::Tree
	ROOT::Hist
	ROOT::Graf
	ROOT::Graf3d
	ROOT::Gpad
	ROOT::Gui
)

add_executable(adjust_ppac_position adjust_ppac_position.cpp)
target_link_libraries(
	adjust_ppac_position
	PRIVATE
	config
	ppac_track
	ROOT::RIO
	ROOT::Tree
)

add_executable(adjust_t0_step1 adjust_t0_step1.cpp)
target_link_libraries(
	adjust_t0_step1
	PRIVATE
	config
	adjust_t0
	ROOT::RIO
	ROOT::Tree
)

add_executable(adjust_t0_step2 adjust_t0_step2.cpp)
target_link_libraries(
	adjust_t0_step2
	PRIVATE
	config
	adjust_t0
	ROOT::RIO
	ROOT::Tree
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

add_executable(extract_d_Li6 d_6Li/extract_d_Li6.cpp)
target_link_libraries(
	extract_d_Li6
	PRIVATE
	config
	dssd_match_event
	ppac_track
	beam_event
	beam_sort
	d_6Li_event
	d_6Li_extract
	physics
	ROOT::RIO
	ROOT::Tree
	ROOT::Hist
)

add_subdirectory(estimate)
add_subdirectory(analysis)
add_subdirectory(rebuild)

add_executable(extract_10C_4He 10C+4He/extract_10C_4He.cpp)
target_link_libraries(
	extract_10C_4He
	PRIVATE
	config
	dssd_match_event
	ppac_track
	beam_event
	beam_sort
	silicon_event
	C10He4_event
	C10He4_extract
	physics
	ROOT::RIO
	ROOT::Tree
	ROOT::Hist
)

add_executable(extract_6Li_two4He C12_d_6Li/extract_6Li_two4He.cpp)
target_link_libraries(
	extract_6Li_two4He
	PRIVATE
	config
	dssd_match_event
	ppac_track
	beam_event
	beam_sort
	silicon_event
	C12D6Li_event
	C12D6Li_extract
	physics
	ROOT::RIO
	ROOT::Tree
	ROOT::Hist
)

add_executable(GUI_12C_d_6Li GUI_12C_d_6Li.cpp)
target_link_libraries(
	GUI_12C_d_6Li
	PRIVATE
	config
	C12D6Li_event
	lise_physics
	ROOT::RIO
	ROOT::Tree
	ROOT::Hist
	ROOT::Graf
	ROOT::Gui
)

add_executable(GUI_10C_4He GUI_10C_4He.cpp)
target_link_libraries(
	GUI_10C_4He
	PRIVATE
	config
	C10He4_event
	rebuild_physics
	physics
	ROOT::RIO
	ROOT::Tree
	ROOT::Hist
	ROOT::Graf
	ROOT::Gui
)

add_executable(GUI_d_Li6 GUI_d_Li6.cpp)
target_link_libraries(
	GUI_d_Li6
	PRIVATE
	config
	d_6Li_event
	rebuild_physics
	physics
	lise_physics
	ROOT::RIO
	ROOT::Tree
	ROOT::Hist
	ROOT::Graf
	ROOT::Gpad
	ROOT::Gui
)

add_executable(sort_beam sort_beam.cpp)
target_link_libraries(
	sort_beam
	PRIVATE
	config
	beam_event
	beam_sort
	ROOT::RIO
	ROOT::Tree
	ROOT::Hist
	ROOT::Graf
	ROOT::Gpad
)

add_executable(GUI_pid GUI_pid.cpp)
target_link_libraries(
	GUI_pid
	PRIVATE
	config
	dssd_match_event
	silicon_event
	d_6Li_extract
	ROOT::RIO
	ROOT::Tree
	ROOT::Hist
	ROOT::Graf
	ROOT::Gpad
	ROOT::Gui
)