cmake_minimum_required(VERSION 3.0)

project(bmove)

# Include SDSL headers from ~/include
include_directories($ENV{HOME}/include)
# Set SDSL library directory to ~/lib
link_directories($ENV{HOME}/lib)

# auxiliary classes, used by both fm and r index
add_library(AuxiliaryClasses src/bitparallelmatrix.cpp src/nucleotide.cpp src/search.cpp)

add_executable(bmove-build src/build.cpp src/buildMove.cpp)

add_executable(bmove-locate src/locate.cpp src/move.cpp src/rindex.cpp src/rsearchstrategy.cpp)
target_link_libraries(bmove-locate AuxiliaryClasses )

add_executable(bmove-benchmarkCharExt src/benchmark.cpp src/move.cpp src/rindex.cpp src/rsearchstrategy.cpp)
target_link_libraries(bmove-benchmarkCharExt AuxiliaryClasses )
target_compile_definitions(bmove-benchmarkCharExt PRIVATE ENABLE_BENCHMARK_FUNCTIONALITY)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -mpopcnt -std=gnu++11")

option(THIRTY_TWO "Use 32-bit types" OFF)
# Pass the value to the source code as a preprocessor definition.
if (THIRTY_TWO)
    add_definitions(-DTHIRTY_TWO)
endif ()

#add_compile_definitions(ALPHABET=5)
add_definitions(-DALPHABET=5)
# ALPHABET = alphabet size


# set the module path
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# set include path for Google's sparse hash table
find_package(SparseHash)
if (SPARSEHASH_FOUND)
    include_directories(${SPARSEHASH_INCLUDE_DIR})
else (SPARSEHASH_FOUND)
    message(FATAL_ERROR "\nFATAL ERROR: The required Google SparseHash package"
            " could not be found on this system.  Please refer to the Detox "
            "manual for the Google Sparsehash installation instructions.  If "
            "you installed Google Sparsehash in a non-standard location "
            "(e.g. somewhere in your homedir), you can point cmake to the "
            "installation location as follows: \ncmake "
            "-DSPARSEHASH_INCLUDE_DIR=<path-to-sparsehash>/include .")
endif(SPARSEHASH_FOUND)

target_link_libraries(bmove-build sdsl)
target_link_libraries(bmove-locate sdsl)
target_link_libraries(bmove-benchmarkCharExt sdsl)

install(TARGETS bmove-build DESTINATION bin)
install(TARGETS bmove-locate DESTINATION bin)
install(TARGETS bmove-benchmarkCharExt DESTINATION bin)

# add_subdirectory(unittest)


