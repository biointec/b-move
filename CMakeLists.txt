# CMake minimum version and project name
cmake_minimum_required(VERSION 3.14)
project(bmove)



# ---- Configuration Options ----

# Option to enable developer mode (default: OFF)
option(DEVELOPER_MODE "Enable developer mode" OFF)
if(DEVELOPER_MODE)
    add_definitions(-DDEVELOPER_MODE)
    message(STATUS "Developer mode enabled.")
else()
    message(STATUS "Developer mode disabled.")
endif()

# Compiler flags for warnings and C++ standard
set(COMMON_FLAGS "-Wall -pedantic -mpopcnt -std=gnu++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COMMON_FLAGS}")

# Option to enable bit packing for the LF move structures (default: ON)
option(LF_MOVE_BIT_PACKED "Enable bit packing the LF move structures" ON)
if(LF_MOVE_BIT_PACKED)
    add_definitions(-DLF_MOVE_BIT_PACKED)
    message(STATUS "Bit packing the LF move structures enabled.")
else()
    message(STATUS "Bit packing the LF move structures disabled.")
endif()

# Option to enable move structures for the phi operations (default: ON)
option(PHI_MOVE "Enable phi move structures" ON)
if(PHI_MOVE)
    add_definitions(-DPHI_MOVE)
    message(STATUS "Phi move structures enabled.")

    # Option to enable balancing the phi move structures (default: ON)
    option(PHI_MOVE_BALANCED "Enable balancing the phi move structures" ON)
    if(PHI_MOVE_BALANCED)
        add_definitions(-DPHI_MOVE_BALANCED)
        message(STATUS "Balancing the phi move structures enabled.")
    else()
        message(STATUS "Balancing the phi move structures disabled.")
    endif()

    # Option to enable bit packing for the phi move structures (default: ON)
    option(PHI_MOVE_BIT_PACKED "Enable bit packing the phi move structures" ON)
    if(PHI_MOVE_BIT_PACKED)
        add_definitions(-DPHI_MOVE_BIT_PACKED)
        message(STATUS "Bit packing the phi move structures enabled.")
    else()
        message(STATUS "Bit packing the phi move structures disabled.")
    endif()

else()
    message(STATUS "Phi move structures disabled.")
endif()

# Always define ALPHABET as 5
add_definitions(-DALPHABET=5)

# ---- Source Files ----

set(SOURCE_FILES
    src/rindexhelpers.cpp
    src/bitparallelmatrix.cpp
    src/rsearchstrategy.cpp
    src/nucleotide.cpp
    src/search.cpp
    src/logger.cpp
    src/parameters.cpp
    src/rindex.cpp
    src/bmove.cpp
    src/move.cpp
    src/movephirepr.cpp
)

# Create the auxiliary library
add_library(AuxiliaryClasses ${SOURCE_FILES})
add_library(AuxiliaryClassesNoReport ${SOURCE_FILES})
add_library(AuxiliaryClassesBenchmarkLF ${SOURCE_FILES})
add_library(AuxiliaryClassesBenchmarkPhi ${SOURCE_FILES})

# ---- Update Git Submodules ----

# Update the Big-BWT submodule
find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
    message(STATUS "Updating Big-BWT submodule")
    execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init external/Big-BWT
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                    RESULT_VARIABLE GIT_SUBMOD_RESULT)
    if(NOT GIT_SUBMOD_RESULT EQUAL "0")
        message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
    endif()
endif()

# ---- Executables and Linking ----

add_executable(bmove-locate src/locate.cpp src/logger.cpp)
target_link_libraries(bmove-locate AuxiliaryClasses)

add_executable(bmove-locate-no-report src/locate.cpp src/logger.cpp)
target_link_libraries(bmove-locate-no-report AuxiliaryClassesNoReport)
target_compile_definitions(bmove-locate-no-report PRIVATE NO_REPORT)
target_compile_definitions(AuxiliaryClassesNoReport PRIVATE NO_REPORT)

add_executable(bmove-benchmarkCharExt src/locate.cpp src/logger.cpp)
target_link_libraries(bmove-benchmarkCharExt AuxiliaryClassesBenchmarkLF)
target_compile_definitions(bmove-benchmarkCharExt PRIVATE LF_BENCHMARK_FUNCTIONALITY)
target_compile_definitions(AuxiliaryClassesBenchmarkLF PRIVATE LF_BENCHMARK_FUNCTIONALITY)
add_executable(bmove-benchmarkPhi src/locate.cpp src/logger.cpp)
target_link_libraries(bmove-benchmarkPhi AuxiliaryClassesBenchmarkPhi)
target_compile_definitions(bmove-benchmarkPhi PRIVATE PHI_BENCHMARK_FUNCTIONALITY)
target_compile_definitions(AuxiliaryClassesBenchmarkPhi PRIVATE PHI_BENCHMARK_FUNCTIONALITY)

add_executable(bmove-build src/build.cpp src/logger.cpp src/buildhelpers.cpp src/move.cpp src/movephirepr.cpp)
# Build the Big-BWT submodule
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/external/Big-BWT/bwtcheck
    COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_CURRENT_SOURCE_DIR}/external/Big-BWT make
    COMMENT "Building Big-BWT executable"
)
add_custom_target(bigbwt_target ALL DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/external/Big-BWT/bwtcheck)
add_dependencies(bmove-build bigbwt_target)

# Copy the script to the build directory after bmove-build is built
add_custom_command(
    TARGET bmove-build POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/src/bmove-build-pfp.sh ${CMAKE_BINARY_DIR}/bmove-build-pfp.sh
    COMMENT "Copying bmove-build_pfp.sh to the build directory"
)

# ---- External Dependencies ----

# FetchContent module for including external projects
include(FetchContent)

# libsais library (added to bmove-build)

# Set specific build options for libsais
message(STATUS "Fetching libsais library")
set(LIBSAIS_USE_OPENMP OFF CACHE BOOL "Disable OpenMP support in libsais")
set(LIBSAIS_BUILD_SHARED_LIB OFF CACHE BOOL "Build libsais as a static library")

FetchContent_Declare(
  libsais
  GIT_REPOSITORY https://github.com/IlyaGrebnov/libsais.git
  GIT_TAG v2.8.3
)
FetchContent_MakeAvailable(libsais)
target_link_libraries(bmove-build libsais)

# parallel-hashmap library
message(STATUS "Fetching parallel-hashmap library")
FetchContent_Declare(
    parallel-hashmap
    GIT_REPOSITORY https://github.com/greg7mdp/parallel-hashmap.git
    GIT_TAG e5b892baed478513adcb6425773cae1eda033057
)
FetchContent_MakeAvailable(parallel-hashmap)
include_directories(${parallel-hashmap_SOURCE_DIR})

# fmt library
message(STATUS "Fetching fmt library")
FetchContent_Declare(
  fmtlib
  GIT_REPOSITORY https://github.com/fmtlib/fmt.git
  GIT_TAG 10.2.1
)
FetchContent_MakeAvailable(fmtlib)
target_link_libraries(AuxiliaryClasses fmt::fmt)
target_link_libraries(AuxiliaryClassesNoReport fmt::fmt)
target_link_libraries(bmove-locate fmt::fmt)
target_link_libraries(bmove-locate-no-report fmt::fmt)
target_link_libraries(bmove-benchmarkCharExt fmt::fmt)
target_link_libraries(bmove-benchmarkPhi fmt::fmt)
target_link_libraries(AuxiliaryClassesBenchmarkLF fmt::fmt)
target_link_libraries(AuxiliaryClassesBenchmarkPhi fmt::fmt)

# FetchContent module for including external projects
include(FetchContent)

# Fetch fmtlib from GitHub
FetchContent_Declare(
  fmtlib
  GIT_REPOSITORY https://github.com/fmtlib/fmt.git
  GIT_TAG 10.2.1 # Adjust the version as needed
)

FetchContent_MakeAvailable(fmtlib)

target_link_libraries(AuxiliaryClasses fmt::fmt)
target_link_libraries(AuxiliaryClassesNoReport fmt::fmt)
target_link_libraries(bmove-locate fmt::fmt)
target_link_libraries(bmove-locate-no-report fmt::fmt)
target_link_libraries(bmove-build fmt::fmt)
target_link_libraries(bmove-benchmarkCharExt fmt::fmt)
target_link_libraries(bmove-benchmarkPhi fmt::fmt)
target_link_libraries(AuxiliaryClassesBenchmarkLF fmt::fmt)
target_link_libraries(AuxiliaryClassesBenchmarkPhi fmt::fmt)

# ---- Build Configuration ----

# Set default build type to Release if not specified
set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# ----Installed Dependencies and Linking ----

# Set the path to your custom FindSDSL.cmake module
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
find_package(SDSL REQUIRED)

if(SDSL_FOUND)
    message(STATUS "Found SDSL library at ${SDSL_LIBRARY}")
    message(STATUS "SDSL include directory: ${SDSL_INCLUDE_DIR}")

    # Include directories
    include_directories(${SDSL_INCLUDE_DIR})
    # Link the SDSL library
    target_link_libraries(bmove-build ${SDSL_LIBRARY})
    target_link_libraries(bmove-locate ${SDSL_LIBRARY})
    target_link_libraries(bmove-locate-no-report ${SDSL_LIBRARY})
    target_link_libraries(bmove-benchmarkCharExt ${SDSL_LIBRARY})
    target_link_libraries(bmove-benchmarkPhi ${SDSL_LIBRARY})
else()
    message(FATAL_ERROR "\nFATAL ERROR: SDSL library not found. Please install it and try again. If "
        "you installed SDSL in a non-standard location, you can point cmake to the "
        "installation location as follows: \ncmake -DSDSL_INCLUDE_DIR=<path-to-sdsl>/include .")
endif()

# ---- Installation ----

# Install targets to the bin directory
install(TARGETS bmove-build bmove-locate bmove-locate-no-report DESTINATION bin)
install(TARGETS bmove-benchmarkCharExt DESTINATION bin)
install(TARGETS bmove-benchmarkPhi DESTINATION bin)

# ---- Unit Tests ----

# Add the unittest subdirectory
add_subdirectory(unittest)