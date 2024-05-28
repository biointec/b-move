cmake_minimum_required(VERSION 3.0)

project(bmove)

# Include custom modules
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Check for SDSL library
find_package(SDSL REQUIRED)
if (NOT SDSL_FOUND)
    message(FATAL_ERROR "SDSL library not found. Please install it and try again.")
else()
    include_directories(${SDSL_INCLUDE_DIRS})
    link_directories(${SDSL_LIBRARY_DIRS})
endif()

# Check for SparseHash library
find_package(SparseHash REQUIRED)
if (NOT SPARSEHASH_FOUND)
    message(FATAL_ERROR "SparseHash library not found. Please install it and try again.")
else()
    include_directories(${SPARSEHASH_INCLUDE_DIR})
endif()

# Auxiliary classes
add_library(AuxiliaryClasses src/bitparallelmatrix.cpp src/nucleotide.cpp src/search.cpp)

add_executable(bmove-build src/build.cpp src/buildMove.cpp)

add_executable(bmove-locate src/locate.cpp src/move.cpp src/rindex.cpp src/rsearchstrategy.cpp)
target_link_libraries(bmove-locate AuxiliaryClasses)

add_executable(bmove-benchmarkCharExt src/benchmark.cpp src/move.cpp src/rindex.cpp src/rsearchstrategy.cpp)
target_link_libraries(bmove-benchmarkCharExt AuxiliaryClasses)
target_compile_definitions(bmove-benchmarkCharExt PRIVATE ENABLE_BENCHMARK_FUNCTIONALITY)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -mpopcnt -std=gnu++11")

option(THIRTY_TWO "Use 32-bit types" OFF)
# Pass the value to the source code as a preprocessor definition.
if (THIRTY_TWO)
    add_definitions(-DTHIRTY_TWO)
endif()

# add_compile_definitions(ALPHABET=5)
add_definitions(-DALPHABET=5)
# ALPHABET = alphabet size

set(default_build_type "Release")
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

target_link_libraries(bmove-build ${SDSL_LIBRARIES})
target_link_libraries(bmove-locate ${SDSL_LIBRARIES})
target_link_libraries(bmove-benchmarkCharExt ${SDSL_LIBRARIES})

install(TARGETS bmove-build DESTINATION bin)
install(TARGETS bmove-locate DESTINATION bin)
install(TARGETS bmove-benchmarkCharExt DESTINATION bin)

# add_subdirectory(unittest)
