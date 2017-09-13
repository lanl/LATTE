# - Find bml
# Find the native BML libraries.
#
#  BML_LIBRARIES    - List of libraries when using bml.
#  BML_FOUND        - True if bml found.
#

find_package(PkgConfig)

pkg_check_modules(PC_BML bml)
find_path(BML_INCLUDE_DIR bml.h HINTS ${PC_BML_INCLUDE_DIRS})

find_library(BML_LIBRARY NAMES bml HINTS ${PC_BML_LIBRARY_DIRS})

set(BML_LIBRARIES ${BML_LIBRARY})
set(BML_INCLUDE_DIRS ${BML_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set BML_FOUND to TRUE
# if all listed variables are TRUE

find_package_handle_standard_args(BML DEFAULT_MSG BML_LIBRARY BML_INCLUDE_DIR)

mark_as_advanced(BML_LIBRARY BML_INCLUDE_DIR)
