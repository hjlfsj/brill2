add_subdirectory(ingot)
add_subdirectory(t0)
add_subdirectory(ppac)

add_library(particle_event particle_event.cpp)
target_include_directories(particle_event PUBLIC ${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
target_link_libraries(particle_event PUBLIC ROOT::Tree ROOT::RIO ROOT::Core)