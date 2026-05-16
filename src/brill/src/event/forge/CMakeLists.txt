add_library(ppac_event ppac_event.cpp)
target_include_directories(ppac_event PUBLIC ${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
target_link_libraries(ppac_event PUBLIC ROOT::Core ROOT::Tree ROOT::RIO)

add_library(tafd_event tafd_event.cpp)
target_include_directories(tafd_event PUBLIC ${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
target_link_libraries(tafd_event PUBLIC ROOT::Core ROOT::Tree ROOT::RIO)

add_library(dssd_event dssd_event.cpp)
target_include_directories(dssd_event PUBLIC ${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
target_link_libraries(dssd_event PUBLIC ROOT::Core ROOT::Tree ROOT::RIO)

add_library(beam_event beam_event.cpp)
target_include_directories(beam_event PUBLIC ${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
target_link_libraries(beam_event PUBLIC ROOT::Core ROOT::Tree ROOT::RIO)

add_library(trigger_event trigger_event.cpp)
target_include_directories(trigger_event PUBLIC ${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
target_link_libraries(trigger_event PUBLIC ROOT::Core ROOT::Tree ROOT::RIO)

add_library(silicon_event silicon_event.cpp)
target_include_directories(silicon_event PUBLIC ${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
target_link_libraries(silicon_event PUBLIC ROOT::Core ROOT::Tree ROOT::RIO)

add_library(csi_event csi_event.cpp)
target_include_directories(csi_event PUBLIC ${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
target_link_libraries(csi_event PUBLIC ROOT::Core ROOT::Tree ROOT::RIO)

add_library(vme_event vme_event.cpp)
target_include_directories(vme_event PUBLIC ${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
target_link_libraries(vme_event PUBLIC ROOT::Core ROOT::Tree ROOT::RIO)

add_library(align_event align_event.cpp)
target_include_directories(align_event PUBLIC ${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
target_link_libraries(align_event PUBLIC ROOT::Core ROOT::Tree ROOT::RIO)
