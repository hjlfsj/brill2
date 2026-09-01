add_library(d_6Li_event d_6Li_event.cpp)
target_include_directories(d_6Li_event PUBLIC ${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
target_link_libraries(d_6Li_event PUBLIC config ROOT::RIO ROOT::Tree)

add_library(d_6Li_extract extract.cpp)
target_include_directories(d_6Li_extract PUBLIC ${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
target_link_libraries(d_6Li_extract PUBLIC dssd_match_event ROOT::Graf)