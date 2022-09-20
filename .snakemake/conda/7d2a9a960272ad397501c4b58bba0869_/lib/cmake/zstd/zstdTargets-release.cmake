#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "zstd::libzstd_shared" for configuration "Release"
set_property(TARGET zstd::libzstd_shared APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(zstd::libzstd_shared PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libzstd.so.1.5.2"
  IMPORTED_SONAME_RELEASE "libzstd.so.1"
  )

list(APPEND _IMPORT_CHECK_TARGETS zstd::libzstd_shared )
list(APPEND _IMPORT_CHECK_FILES_FOR_zstd::libzstd_shared "${_IMPORT_PREFIX}/lib/libzstd.so.1.5.2" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
