add_executable(estimate_normalize estimate_normalize.cpp)
target_link_libraries(
	estimate_normalize
	PRIVATE
	config
	t0_dssd
	ROOT::Graf
)

add_executable(estimate_t0_center estimate_t0_center.cpp)
target_link_libraries(
	estimate_t0_center
	PRIVATE
	config
	t0_dssd
	ROOT::Graf
)

add_executable(estimate_t0_pid estimate_t0_pid.cpp)
target_link_libraries(
	estimate_t0_pid
	PRIVATE
	config
	silicon_event
	t0_dssd
	ROOT::Graf
)

add_executable(estimate_t0_straight estimate_t0_straight.cpp)
target_link_libraries(
	estimate_t0_straight
	PRIVATE
	config
	silicon_event
	t0_dssd
	ROOT::Graf
)

add_executable(estimate_t0_csi_pid estimate_t0_csi_pid.cpp)
target_link_libraries(
	estimate_t0_csi_pid
	PRIVATE
	config
	silicon_event
	ROOT::Graf
)