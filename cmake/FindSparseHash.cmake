# FindSparseHash.cmake

# Locate SparseHash library headers
# Once done, this will define
# SPARSEHASH_FOUND - system has SparseHash
# SPARSEHASH_INCLUDE_DIR - the SparseHash include directories

find_path(SPARSEHASH_INCLUDE_DIR
    NAMES google/sparsehash/sparsehashtable.h
    HINTS ENV CPATH ENV HOME ENV C_INCLUDE_PATH ENV CPLUS_INCLUDE_PATH /usr/include /usr/local/include
    PATH_SUFFIXES include
    NO_DEFAULT_PATH)

if (SPARSEHASH_INCLUDE_DIR)
    set(SPARSEHASH_FOUND TRUE)
else ()
    set(SPARSEHASH_FOUND FALSE)
endif ()

if (SPARSEHASH_FOUND)
    message(STATUS "Found SparseHash include directory: ${SPARSEHASH_INCLUDE_DIR}")
else ()
    message(STATUS "Could not find SparseHash library. Ensure that the paths are set correctly.")
    message(STATUS "Searched in:")
    message(STATUS "  INCLUDE_DIR: ${SPARSEHASH_INCLUDE_DIR}")
endif ()

mark_as_advanced(SPARSEHASH_INCLUDE_DIR)
