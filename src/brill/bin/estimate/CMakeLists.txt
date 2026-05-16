add_executable(estimate_normalize estimate_normalize.cpp)
target_link_libraries(
	estimate_normalize
	PRIVATE
	estimate_support
)

add_executable(estimate_t0_center estimate_t0_center.cpp)
target_link_libraries(
	estimate_t0_center
	PRIVATE
	estimate_support
)

add_executable(estimate_t0_pid estimate_t0_pid.cpp)
target_link_libraries(
	estimate_t0_pid
	PRIVATE
	estimate_support
	silicon_event
)

add_executable(estimate_t0_straight estimate_t0_straight.cpp)
target_link_libraries(
	estimate_t0_straight
	PRIVATE
	estimate_support
	silicon_event
	ROOT::Graf
)
