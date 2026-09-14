#----------------------------------------------------------------
# Generated CMake target import file for configuration "Debug".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "catima::catima" for configuration "Debug"
set_property(TARGET catima::catima APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(catima::catima PROPERTIES
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/libcatima.so"
  IMPORTED_SONAME_DEBUG "libcatima.so"
  )

list(APPEND _cmake_import_check_targets catima::catima )
list(APPEND _cmake_import_check_files_for_catima::catima "${_IMPORT_PREFIX}/lib/libcatima.so" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
