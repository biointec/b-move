cmake_minimum_required(VERSION 3.14)

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

# Auxiliary classes
add_library(AuxiliaryClasses src/bitparallelmatrix.cpp src/nucleotide.cpp src/search.cpp)

add_executable(bmove-build src/build.cpp src/buildhelpers.cpp src/move.cpp)

# Define the submodule directory
set(SUBMODULE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/Big-BWT)

# Add a custom command to build the submodule executable
add_custom_command(
    OUTPUT ${SUBMODULE_DIR}/bwtcheck
    COMMAND ${CMAKE_COMMAND} -E chdir ${SUBMODULE_DIR} make
    WORKING_DIRECTORY ${SUBMODULE_DIR}
    COMMENT "Building Big-BWT executable"
)

# Add a custom target that depends on the custom command
add_custom_target(bigbwt_target ALL DEPENDS ${SUBMODULE_DIR}/bwtcheck)

# Ensure bmove-build depends on bigbwt_target
add_dependencies(bmove-build bigbwt_target)

# Define the source and destination paths
set(SCRIPT_SOURCE_DIR ${CMAKE_SOURCE_DIR}/src)
set(SCRIPT_DEST_DIR ${CMAKE_BINARY_DIR})

# Copy the script to the build directory
add_custom_command(
    TARGET bmove-build POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${SCRIPT_SOURCE_DIR}/bmove-build-pfp.sh ${SCRIPT_DEST_DIR}/bmove-build-pfp.sh
    COMMENT "Copying bmove-build-pfp.sh to the build directory"
)

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

# add libsais to bmove-build
# add dependencies
include(cmake/CPM.cmake)
CPMAddPackage(
  NAME libsais
  GITHUB_REPOSITORY IlyaGrebnov/libsais
  GIT_TAG v2.8.3
  OPTIONS
    "LIBSAIS_USE_OPENMP OFF"
    "LIBSAIS_BUILD_SHARED_LIB OFF"
)

target_link_libraries(bmove-build libsais)

target_link_libraries(bmove-build ${SDSL_LIBRARIES})
target_link_libraries(bmove-locate ${SDSL_LIBRARIES})
target_link_libraries(bmove-benchmarkCharExt ${SDSL_LIBRARIES})

install(TARGETS bmove-build DESTINATION bin)
install(TARGETS bmove-locate DESTINATION bin)
install(TARGETS bmove-benchmarkCharExt DESTINATION bin)

# add_subdirectory(unittest)
