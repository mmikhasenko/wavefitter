PROJECT(POLE-SEARCH)

cmake_minimum_required(VERSION 3.5)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmakeModules")

# ROOT
find_package(ROOT REQUIRED)
include_directories(${ROOT_INCLUDE_DIR})
LINK_DIRECTORIES(${ROOT_LIBRARY_DIR})

# ublas
FIND_PACKAGE(Boost 1.53.0 REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})

# setup libconfig
# environment variable $LIBCONFIG is expected to point to libconfig location
message(STATUS "")
message(STATUS "Looking for libconfig library...")
FIND_PACKAGE(Libconfig 1.4 REQUIRED)

# compilation flags:
#   -Wall to print all warnings
#   -D_GNU_SOURCE
#   -fexceptions to enable exception handling, while I don't use it
#   -pthread to support multithreading, while I don't use it
#   -std=gnu++11 to use c++11 
#   -Wno-c99-extensions to remove warning about c++11
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -D_GNU_SOURCE -fexceptions -pthread -std=gnu++11 -Wno-c99-extensions" )

# include directories
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/src/rpwa)
include_directories(${ROOT_INCLUDE_DIR})


add_subdirectory(src)
add_subdirectory(kmatrix)
