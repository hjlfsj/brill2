add_library(C10He4_event 10C_4He_event.cpp)
target_include_directories(C10He4_event PUBLIC ${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
target_link_libraries(C10He4_event PUBLIC config ROOT::RIO ROOT::Tree)

add_library(C10He4_tafd tafd.cpp)
target_include_directories(C10He4_tafd PUBLIC ${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
target_link_libraries(C10He4_tafd PUBLIC config tafd_event ROOT::RIO ROOT::Tree)

add_library(C10He4_tafcsi tafcsi.cpp)
target_include_directories(C10He4_tafcsi PUBLIC ${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
target_link_libraries(C10He4_tafcsi PUBLIC ROOT::Tree ROOT::RIO)

add_library(C10He4_extract extract.cpp)
target_include_directories(C10He4_extract PUBLIC ${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
target_link_libraries(C10He4_extract PUBLIC dssd_match_event silicon_event ROOT::Graf)