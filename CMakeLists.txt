# minimum cmake version required
cmake_minimum_required(VERSION 3.9)

# name, version, description and language of the project
project(McDIPPER VERSION 1.0
        DESCRIPTION "Monte-Carlo Dipole Event-geneRator"
        LANGUAGES CXX)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
    message(FATAL_ERROR "You don't want to configure in the source directory!")
endif()

if(NOT "{CMAKE_CXX_STANDARD}")
    set(CMAKE_CXX_STANDARD 17)    
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fPIC -O3")
message(STATUS "CXX_FLAGS = " ${CMAKE_CXX_FLAGS})

# Option for specifying GSL library path via command line
option(GSL_PATH "Path to GSL library")

# Look for GSL package
find_package(GSL REQUIRED)
if(GSL_FOUND)
    message(STATUS "GSL found: ${GSL_INCLUDE_DIRS}")
    include_directories(${GSL_INCLUDE_DIRS})
    list(APPEND LIBRARIES_TO_LINK ${GSL_LIBRARIES})
else()
    message(STATUS "GSL not found. Using provided library path or defaults.")
    if(NOT GSL_PATH)
        message(FATAL_ERROR "GSL not found. Please specify the path to GSL library using -DGSL_PATH=<path>")
    endif()
endif()


# Specifications for the Cuba library
# Retrieve the CUBA_PATH variable if provided; otherwise, default to an empty string
if(NOT DEFINED CUBA_PATH)
    message(FATAL_ERROR "-DCUBA_PATH=<path> not given")
endif()

# Add Cuba library and source directories to the include path using the provided path
include_directories(${CUBA_PATH})
include_directories(${CUBA_PATH}/src)
list(APPEND LIBRARIES_TO_LINK ${CUBA_PATH}/libcuba.a)




# Option for specifying LHAPDF path via command line
option(LHAPDF_PATH "Path to LHAPDF library")
if(NOT LHAPDF_FOUND)
    if(NOT LHAPDF_PATH)
        message(FATAL_ERROR "LHAPDF not found. Please specify the path to LHAPDF library using -DLHAPDF_PATH=<path>")
    endif()

    # Include directories and link LHAPDF library if not found by find_package
    include_directories(${LHAPDF_PATH}/include)
    list(APPEND LIBRARIES_TO_LINK ${LHAPDF_PATH}/lib/libLHAPDF.a) # Adjust library file if needed
endif()


# HDF5 (C API, used by h5_out.cpp)
find_package(HDF5 REQUIRED COMPONENTS C)
include_directories(${HDF5_INCLUDE_DIRS})
list(APPEND LIBRARIES_TO_LINK ${HDF5_LIBRARIES})

# OpenMP. macOS: Apple clang needs libomp (brew install libomp) and -Xpreprocessor
if(APPLE AND NOT OpenMP_CXX_FLAGS)
    execute_process(COMMAND brew --prefix libomp
                    OUTPUT_VARIABLE LIBOMP_PREFIX
                    OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    if(LIBOMP_PREFIX)
        set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I${LIBOMP_PREFIX}/include")
        set(OpenMP_CXX_LIB_NAMES "omp")
        set(OpenMP_omp_LIBRARY "${LIBOMP_PREFIX}/lib/libomp.dylib")
    endif()
endif()
find_package(OpenMP REQUIRED)

message(STATUS "LIBRARIES_TO_LINK = " ${LIBRARIES_TO_LINK})


# Source files
set(SRC
    src/main.cpp
    src/config.cpp
    src/config_mv.cpp
    src/pdfs_lhapdf.cpp
    src/nucleus.cpp
    src/event.cpp
    src/charges.cpp
    src/routines.cpp
    src/dipole_ipsat.cpp
    src/dipole_mv.cpp
    src/model_gbw.cpp
    src/model_ipsat.cpp
    src/model_mv.cpp
    src/random.cpp
    src/h5_out.cpp
    # Add other source files here
)

file(COPY ${CMAKE_SOURCE_DIR}/ipsat DESTINATION  ${CMAKE_BINARY_DIR})

# Create the executable
add_executable(mcDipper ${SRC})



# Link libraries
target_link_libraries(mcDipper ${LIBRARIES_TO_LINK} OpenMP::OpenMP_CXX) # Add more libraries here if needed







