# - Find progress
# Find the native PROGRESS libraries.
#
#  PROGRESS_LIBRARIES    - List of libraries when using progress.
#  PROGRESS_FOUND        - True if progress found.
#

find_package(PkgConfig)

pkg_check_modules(PC_PROGRESS progress)
find_path(PROGRESS_INCLUDE_DIR prg_pulaymixer_mod.mod HINTS ${PC_PROGRESS_INCLUDE_DIRS})

find_library(PROGRESS_LIBRARY NAMES progress HINTS ${PC_PROGRESS_LIBRARY_DIRS})

set(PROGRESS_LIBRARIES ${PROGRESS_LIBRARY})
set(PROGRESS_INCLUDE_DIRS ${PROGRESS_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set PROGRESS_FOUND to TRUE
# if all listed variables are TRUE

find_package_handle_standard_args(PROGRESS DEFAULT_MSG PROGRESS_LIBRARY PROGRESS_INCLUDE_DIR)

mark_as_advanced(PROGRESS_LIBRARY PROGRESS_INCLUDE_DIR)
