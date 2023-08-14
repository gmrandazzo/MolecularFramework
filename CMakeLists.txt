cmake_minimum_required(VERSION 3.18)
project(MolecularFramework C CXX)
set(AUTHOR "Giuseppe Marco Randazzo")
set(VERSION "0.1.0")

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules/")
#find_package(RDKit)
find_package(Libscientific)
#message(STATUS "${LIBSCIENTIFIC_INCLUDE_DIR}")

set(CMAKE_MACOSX_RPATH 1)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build, options are: Debug Release" FORCE)
endif(NOT CMAKE_BUILD_TYPE)

if(NOT CMAKE_INSTALL_PREFIX)
  message(SEND_ERROR "Please specify the installation prefix -DCMAKE_INSTALL_PREFIX=...")
endif()

if(WIN32)
  message(STATUS "Build Release: WIN32")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu99 -O2 -pipe -mwindows -pedantic -Wall -static -DWIN32")
  set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS}")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -s")
  set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "")
else(NOT WIN32)
  string(TOLOWER "${CMAKE_BUILD_TYPE}" cmake_build_type_tolower)
  if(cmake_build_type_tolower STREQUAL "debug")
    set(CMAKE_BUILD_TYPE "Debug")
    #set(CMAKE_C_FLAGS "-O0 -pedantic -Werror -g -DDEBUG")
    set(CMAKE_C_FLAGS "-O0  -Werror -g -fsanitize=address -fprofile-arcs -ftest-coverage ")
    set(CMAKE_CXX_FLAGS "-O0  -Werror -g -fsanitize=address -fprofile-arcs -ftest-coverage")
    set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "")
  else()
    set(CMAKE_BUILD_TYPE "Release")
      if(NOT CMAKE_C_FLAGS)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2 -Werror")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -Werror")
      endif()
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,-O1")
  endif()
endif(WIN32)

enable_testing()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu99 -Wall")

message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")

link_directories("/usr/lib" "/usr/local/lib" ${LIBSCIENTIFIC_LIB_DIR})
include_directories("/usr/lib" "/usr/local/include" ${LIBSCIENTIFIC_INCLUDE_DIR})

#include_directories(${RDKIT_INCLUDE_DIR})

if(cmake_build_type_tolower STREQUAL "debug")
  include(Dart)
  cmake_policy(SET CMP0145 OLD)
  enable_testing()
endif()


add_subdirectory(src)

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake" IMMEDIATE @ONLY)
add_custom_target(uninstall "${CMAKE_COMMAND}" -P
        "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake")
