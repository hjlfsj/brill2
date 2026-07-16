add_library(ppac_track ppac_track.cpp ppac_track_event.cpp)
target_include_directories(ppac_track PUBLIC ${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
target_link_libraries(ppac_track PUBLIC config ppac_event ROOT::Graf ROOT::Core)