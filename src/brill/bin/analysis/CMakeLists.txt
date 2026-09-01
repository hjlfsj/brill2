add_executable(calc_velocity calc_velocity.cpp)
target_include_directories(
	calc_velocity
	PRIVATE
	${PROJECT_SOURCE_DIR}
)
target_link_libraries(
	calc_velocity
	PRIVATE
	rebuild_physics
)