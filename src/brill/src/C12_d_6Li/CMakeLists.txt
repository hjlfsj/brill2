add_library(C12D6Li_event C12_d_6Li_event.cpp)
target_include_directories(C12D6Li_event PUBLIC ${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
target_link_libraries(C12D6Li_event PUBLIC config ROOT::RIO ROOT::Tree)

add_library(C12D6Li_extract extract.cpp)
target_include_directories(C12D6Li_extract PUBLIC ${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
target_link_libraries(C12D6Li_extract PUBLIC dssd_match_event ROOT::Graf)