cmake_minimum_required(VERSION 3.18 FATAL_ERROR)

if( CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR )
    message( FATAL_ERROR "Please select another Build Directory ! (and give it a clever name, like bin_Visual2012_64bits/)" )
endif()
if( CMAKE_SOURCE_DIR MATCHES " " )
	message( "Your Source Directory contains spaces. If you experience problems when compiling, this can be the cause." )
endif()
if( CMAKE_BINARY_DIR MATCHES " " )
	message( "Your Build Directory contains spaces. If you experience problems when compiling, this can be the cause." )
endif()

#project name
project(RatRace LANGUAGES CXX CUDA)

#c++ standard
set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED True)
set (CMAKE_POSITION_INDEPENDENT_CODE ON)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

#cuda
if (CMAKE_CUDA_COMPILER)
    MESSAGE(STATUS "Generate with CUDA")
    enable_language(CUDA)
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-extended-lambda")
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-relaxed-constexpr")
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --use_fast_math")
    set(CMAKE_CUDA_ARCHITECTURES 61)    # adapt CUDA architecture to your GPU
else()
    MESSAGE(FATAL_ERROR "Cuda not found")
endif()

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH} )

########## required libraries ##########

#glog
find_package(glog REQUIRED)
set(LIBS ${LIBS} glog::glog)

#OpenMesh
find_package(OpenMesh REQUIRED)
set(LIBS ${LIBS} ${OPENMESH_LIBRARIES})
add_definitions(-D_USE_MATH_DEFINES)

#OMP
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    set(LIBS ${LIBS} OpenMP::OpenMP_CXX)
endif()


########## executables ##########

set(SRC ${PROJECT_SOURCE_DIR}/src)

add_subdirectory(src)