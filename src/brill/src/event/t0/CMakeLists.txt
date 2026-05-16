add_library(dssd_match_event dssd_match_event.cpp)
target_include_directories(dssd_match_event ${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
target_link_libraries(dssd_match_event PUBLIC ROOT::Tree ROOT::RIO ROOT::Core)