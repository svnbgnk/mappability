# Minimum cmake version
cmake_minimum_required (VERSION 3.0.0)

# Name of project and that it is C++ only.
project (mappability CXX)

#set (CMAKE_BUILD_TYPE "Release")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

set (CXX11_FOUND "1")

# Search for zlib as a dependency for SeqAn.
find_package (ZLIB)

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_PREFIX_PATH "~/Dev/seqan/util/cmake")
set (SEQAN_INCLUDE_PATH "~/Dev/seqan/include")

# Load the SeqAn module and fail if not found.
find_package (OpenMP REQUIRED)
find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "-std=c++14 ${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# Add executable and link against SeqAn dependencies.
add_executable (my_project my_project.cpp)
target_link_libraries (my_project ${SEQAN_LIBRARIES})
