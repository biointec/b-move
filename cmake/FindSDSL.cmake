# FindSDSL.cmake

# Locate SDSL library and headers

find_path(SDSL_INCLUDE_DIR
    NAMES sdsl
    HINTS ENV CPATH ENV HOME ENV C_INCLUDE_PATH ENV CPLUS_INCLUDE_PATH /usr/include /usr/local/include
    PATH_SUFFIXES include
    NO_DEFAULT_PATH)

find_library(SDSL_LIBRARY
    NAMES sdsl
    HINTS ENV LD_LIBRARY_PATH ENV HOME ENV LIBRARY_PATH /usr/lib /usr/local/lib
    PATH_SUFFIXES lib
    NO_DEFAULT_PATH)

if (SDSL_INCLUDE_DIR AND SDSL_LIBRARY)
    set(SDSL_FOUND TRUE)
    set(SDSL_LIBRARIES ${SDSL_LIBRARY})
    set(SDSL_INCLUDE_DIRS ${SDSL_INCLUDE_DIR})
else ()
    set(SDSL_FOUND FALSE)
endif ()

if (SDSL_FOUND)
    message(STATUS "Found SDSL library: ${SDSL_LIBRARY}")
    message(STATUS "Found SDSL include directory: ${SDSL_INCLUDE_DIR}")
else ()
    message(STATUS "Could not find SDSL library. Ensure that the paths are set correctly.")
    message(STATUS "Searched in:")
    message(STATUS "  INCLUDE_DIR: ${SDSL_INCLUDE_DIR}")
    message(STATUS "  LIBRARY: ${SDSL_LIBRARY}")
endif ()

mark_as_advanced(SDSL_INCLUDE_DIR SDSL_LIBRARY)
