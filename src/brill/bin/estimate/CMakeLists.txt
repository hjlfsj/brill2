add_executable(estimate_normalize estimate_normalize.cpp)
target_link_libraries(
	estimate_normalize
	PRIVATE
	estimate_support
)
