add_library(dssd_match_event dssd_match_event.cpp)
target_include_directories(dssd_match_event PUBLIC ${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
target_link_libraries(dssd_match_event PUBLIC ROOT::Tree ROOT::RIO ROOT::Core)

add_library(dssd_match dssd_match.cpp)
target_include_directories(dssd_match PUBLIC ${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
target_link_libraries(dssd_match PUBLIC dssd_match_event ROOT::Tree ROOT::RIO ROOT::Core)

add_library(t0_event t0_event.cpp)
target_include_directories(t0_event PUBLIC ${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
target_link_libraries(t0_event PUBLIC ROOT::Tree ROOT::RIO ROOT::Core)