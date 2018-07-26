# Minimum cmake version
cmake_minimum_required (VERSION 3.0.0)

# Name of project and that it is C++ only.
project (mappability CXX)

#set (CMAKE_BUILD_TYPE "Debug") #Release
set (CUR_DIR "${CMAKE_CURRENT_LIST_DIR}/..")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

set (CXX11_FOUND "1")

# Search for zlib as a dependency for SeqAn.
find_package (ZLIB)

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_PREFIX_PATH "/home/sven/devel/seqan/util/cmake")
set (SEQAN_INCLUDE_PATH "/home/sven/devel/seqan/include")

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
set (CMAKE_CXX_FLAGS "-std=c++14 -pthread -DSEQAN_DISABLE_VERSION_CHECK=YES -I ${CUR_DIR}/sdsl-lite/include ${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")
# set (CMAKE_CXX_FLAGS "-std=c++14 ${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} -I/home/sven/devel/sdsl-lite/include") # -ldivsufsort -ldivsufsort64

# Add executable and link against SeqAn dependencies.
add_executable (mappability mappability.cpp)
target_link_libraries (mappability ${SEQAN_LIBRARIES})

add_executable (tests tests.cpp)
target_link_libraries (tests ${SEQAN_LIBRARIES})

add_executable (create_index create_index.cpp)
target_link_libraries (create_index ${SEQAN_LIBRARIES})

add_executable (create_bit_v create_bit_v.cpp)
target_link_libraries (create_bit_v ${SEQAN_LIBRARIES})

add_executable (search search.cpp)
target_link_libraries (search ${SEQAN_LIBRARIES})

add_executable (heatmap heatmap.cpp)
target_link_libraries (heatmap ${SEQAN_LIBRARIES})

add_executable (benchmark benchmark.cpp)
target_link_libraries (benchmark ${SEQAN_LIBRARIES})


