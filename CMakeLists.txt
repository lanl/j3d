# (c) 2025. Triad National Security, LLC. All rights reserved.
# This program was produced under U.S. Government contract 89233218CNA000001
# for Los Alamos National Laboratory (LANL), which is operated by Triad National
# Security, LLC for the U.S. Department of Energy/National Nuclear Security
# Administration. All rights in the program are reserved by Triad National Security,
# LLC, and the U.S. Department of Energy/National Nuclear Security Administration.
# The Government is granted for itself and others acting on its behalf a nonexclusive,
# paid-up, irrevocable worldwide license in this material to reproduce, prepare.
# derivative works, distribute copies to the public, perform publicly and display
# publicly, and to permit others to do so.

cmake_minimum_required(VERSION 3.13)

project(j3d CXX)   # CXX needed for one of the tests

# version
set(J3D_VERSION_MAJOR 1)
set(J3D_VERSION_MINOR 2)
set(J3D_VERSION_PATCH 2)
set(CMAKE_CXX_EXTENSIONS OFF)

# J3D compile time variables
set(J3D_MAX_VERTS 256 CACHE STRING "Maximum number of vertices in a polyhedron")
set(J3D_MAX_EDGES 512 CACHE STRING "Maximum number of edges in a polyhedron")
set(J3D_MAX_ORDER 5 CACHE STRING "Maximum order for moments calculation of polyhedra")
# J2D compile time variables
set(J2D_MAX_VERTS 256 CACHE STRING "Maximum number of vertices in a polygon")
set(J2D_MAX_ORDER 5 CACHE STRING "Maximum order for moments calculation of polygon")
# Other options
option(J3D_ENABLE_Kokkos "Enable Kokkos" OFF)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
message(STATUS "CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")
set(CMAKE_CXX_FLAGS_COVERAGE "-g -Og --coverage")

if (CMAKE_BUILD_TYPE STREQUAL "Debug"
  OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo"
  OR CMAKE_BUILD_TYPE STREQUAL "Coverage")
  add_definitions(-DDEBUG)
endif()

# Unit testing
if (ENABLE_UNIT_TESTS)
  enable_testing()
endif ()

if (USE_CLANG_FORMAT)
  include(ClangFormat)
  file(GLOB_RECURSE ALL_SOURCE_FILES CONFIGURE_DEPENDS "*.h" "*.c")
  clangformat_setup(${ALL_SOURCE_FILES})
endif()

add_library(j3d)

# Kokkos
if (J3D_ENABLE_Kokkos)
  find_package(Kokkos REQUIRED)
  target_link_libraries(j3d Kokkos::kokkos)
endif()

# Alias (Daniel Pfeiffer, Effective CMake) - this allows other
# projects that use Pkg as a subproject to find_package(Nmspc::Pkg)
# which does nothing because Pkg is already part of the project

add_library(j3d::j3d ALIAS j3d)
set(j3d_LIBRARIES j3d::j3d CACHE STRING "J3D library target")


add_subdirectory(src)


# Wonton targets

install(TARGETS j3d
  EXPORT j3d_LIBRARIES
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include
  INCLUDES DESTINATION include
  )

#-----------------------------------------------------------------------------
# Prepare output for configuration files to be used by projects importing Wonton
#-----------------------------------------------------------------------------

# Write a configuration file from template replacing only variables enclosed
# by the @ sign.

configure_file(${PROJECT_SOURCE_DIR}/cmake/j3dConfig.cmake.in
  j3dConfig.cmake @ONLY)
install(FILES ${PROJECT_BINARY_DIR}/j3dConfig.cmake
  DESTINATION lib/cmake/j3d)


# export targets

install(EXPORT j3d_LIBRARIES
  FILE j3dTargets.cmake
  NAMESPACE j3d::
  EXPORT_LINK_INTERFACE_LIBRARIES
  DESTINATION lib/cmake/j3d)



# Dynamically configured header files that contains defines like
# J3D_MAX_VERTS if enabled.

configure_file(${PROJECT_SOURCE_DIR}/cmake/j3d-config.h.in
  ${PROJECT_BINARY_DIR}/j3d-config.h @ONLY)
install(FILES ${PROJECT_BINARY_DIR}/j3d-config.h
  DESTINATION include)

