add_library(three_alpha_event three_alpha_event.cpp)
target_include_directories(three_alpha_event PUBLIC ${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
target_link_libraries(three_alpha_event PUBLIC config ROOT::RIO ROOT::Tree)

add_library(three_alpha_extract extract.cpp)
target_include_directories(three_alpha_extract PUBLIC ${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
target_link_libraries(three_alpha_extract PUBLIC dssd_match_event energy_calculator ROOT::Graf)