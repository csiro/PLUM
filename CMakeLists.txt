cmake_minimum_required(VERSION 3.0.2)

set(CMAKE_CXX_COMPILER "g++")

project(Plum
  VERSION 0.1.0
  LANGUAGES CXX)

#------------------------------------------------------------------------------#
# Compilation options

option(STATIC "compile with the -static linker flag" OFF)

#------------------------------------------------------------------------------#
# Add libraries

add_subdirectory(lime)

#------------------------------------------------------------------------------#
# Compiler Configuration

## Static build
if(STATIC)
  set(CMAKE_EXE_LINKER_FLAGS "-static")
endif()

## Build Type
## Note that Visual Studio can build multiple build types at once and
## thus the build type need to be specified with the parameter
## --config when building, e.g.,
##   cmake --build . --config Release
set(DEFAULT_BUILD_TYPE "Release")
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${DEFAULT_BUILD_TYPE}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}" CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
	"Debug"
	"Release"
	"MinSizeRel"
	"RelWithDebInfo"
	)
endif()
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to '${DEFAULT_BUILD_TYPE}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}")
endif()

## C++ Standard
set(CMAKE_CXX_STANDARD 20)

add_compile_definitions (PLUM_GUROBI)
add_compile_definitions (PLUM_HIGHS)
add_compile_definitions (PLUM_LPSOLVE)

message(STATUS "  Gurobi home is $ENV{GUROBI_HOME}")
message(STATUS "lp_solve home is $ENV{LP_SOLVE_HOME}")
set(SOLVER_INCLUDE_DIRS "")
set(SOLVER_LIB_DIRS "")
if(NOT $ENV{GUROBI_HOME} STREQUAL "")
  message(STATUS "gurobi_home is set")
  set(SOLVER_INCLUDE_DIRS ${SOLVER_INCLUDE_DIRS} $ENV{GUROBI_HOME}/include)
  set(SOLVER_LIB_DIRS ${SOLVER_LIB_DIRS} $ENV{GUROBI_HOME}/lib)
  set(SOLVER_LIB ${SOLVER_LIB} gurobi_g++5.2 gurobi100)
endif()

if(NOT $ENV{LP_SOLVE_HOME} STREQUAL "")
  message(STATUS "lp_solver_home is set")
  set(SOLVER_INCLUDE_DIRS ${SOLVER_INCLUDE_DIRS} $ENV{LP_SOLVE_HOME})
  set(SOLVER_LIB_DIRS ${SOLVER_LIB_DIRS} $ENV{LP_SOLVE_HOME}/lpsolve55/bin/ux64)
  set(SOLVER_LIB ${SOLVER_LIB} lpsolve55)
endif()

if(NOT $ENV{HIGHS_HOME} STREQUAL "")
  message(STATUS "highs_home is set")
  set(SOLVER_INCLUDE_DIRS ${SOLVER_INCLUDE_DIRS} $ENV{HIGHS_HOME} $ENV{HIGHS_HOME}/src $ENV{HIGHS_HOME}/build)
  set(SOLVER_LIB_DIRS ${SOLVER_LIB_DIRS} $ENV{HIGHS_HOME}/build/lib)
  set(SOLVER_LIB ${SOLVER_LIB} highs)
endif()

message(STATUS "Using solver include dir(s) ${SOLVER_INCLUDE_DIRS}")
message(STATUS "Using solver lib dir(s) ${SOLVER_LIB_DIRS}")


#------------------------------------------------------------------------------#
# Set directories

set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
set(PROJECT_INCLUDE_REL_DIR include)
#file(RELATIVE_PATH PROJECT_INCLUDE_REL_DIR
#     ${CMAKE_CURRENT_SOURCE_DIR}
#     ${PROJECT_SOURCE_DIR}/.)
   
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

include_directories(
  ${PROJECT_INCLUDE_REL_DIR} ${LIME_INCLUDE_DIRS} ${SOLVER_INCLUDE_DIRS}
  )

message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
message(STATUS "CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "PROJECT_INCLUDE_DIR: ${PROJECT_INCLUDE_DIR}")
message(STATUS "PROJECT_INCLUDE_REL_DIR: ${PROJECT_INCLUDE_REL_DIR}")

#------------------------------------------------------------------------------#
# Plum library

set(PLUM_SOURCE_CFILES
  ${PROJECT_SOURCE_DIR}/src/dummysolver.cpp
  ${PROJECT_SOURCE_DIR}/src/grblpsolverimp.cpp
  ${PROJECT_SOURCE_DIR}/src/lnsmxsolver.cpp
  ${PROJECT_SOURCE_DIR}/src/highslpsolverimp.cpp
  ${PROJECT_SOURCE_DIR}/src/highsutil.cpp
  ${PROJECT_SOURCE_DIR}/src/incrsolver.cpp
  ${PROJECT_SOURCE_DIR}/src/intsolver.cpp
  ${PROJECT_SOURCE_DIR}/src/lpslpsolverimp.cpp
  ${PROJECT_SOURCE_DIR}/src/lpsutil.cpp
  ${PROJECT_SOURCE_DIR}/src/lpsolver.cpp
  ${PROJECT_SOURCE_DIR}/src/multisol.cpp
  ${PROJECT_SOURCE_DIR}/src/params.cpp
#  ${PROJECT_SOURCE_DIR}/src/mathheur.cpp
  ${PROJECT_SOURCE_DIR}/src/pathfinder.cpp
  ${PROJECT_SOURCE_DIR}/src/scenario.cpp
  ${PROJECT_SOURCE_DIR}/src/solution.cpp
  )

add_library(plumlib 
  ${PLUM_SOURCE_CFILES}
  )
target_link_libraries(plumlib lime)

#------------------------------------------------------------------------------#
# Executable
add_custom_target(
  show_which_make ALL
  )
add_custom_command(
  OUTPUT show_which_make
  COMMAND cmake -E cmake_echo_color --cyan "Build type is ${CMAKE_BUILD_TYPE}"
  )

add_custom_target(
  touch_builddate ALL
  )
add_custom_command(
  OUTPUT touch_builddate # Never exists, so this always runs...
  COMMAND cmake -E touch ${PROJECT_SOURCE_DIR}/include/mosh/builddate.h
  )

add_custom_command (
  OUTPUT __tags # Never exists, so this always runs...
  COMMAND find src include -name "*.cpp" -o -name "*.h" | etags -Q -o tags -
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  VERBATIM
  )  

add_executable(plum
  show_which_make touch_builddate __tags
  ${PROJECT_SOURCE_DIR}/src/plum.cpp
  )

add_executable(plumsp
  show_which_make touch_builddate __tags
  ${PROJECT_SOURCE_DIR}/src/plumsp.cpp
  )

add_executable(plumchk
  show_which_make touch_builddate __tags
  ${PROJECT_SOURCE_DIR}/src/plumchk.cpp
  )

add_executable(plummx
  show_which_make touch_builddate __tags
  ${PROJECT_SOURCE_DIR}/src/plummx.cpp
  )

#add_executable(plumseq
#  show_which_make touch_builddate __tags
#  ${PROJECT_SOURCE_DIR}/src/plumseq.cpp
#  )

add_executable(plumcmp
  show_which_make touch_builddate __tags
  ${PROJECT_SOURCE_DIR}/src/plumcmp.cpp
  )

add_executable(plummerge
  show_which_make touch_builddate __tags
  ${PROJECT_SOURCE_DIR}/src/plummerge.cpp
  )

target_link_directories(plum PUBLIC ${SOLVER_LIB_DIRS} )
target_link_libraries(plum plumlib ${SOLVER_LIB} )

target_link_directories(plummx PUBLIC ${SOLVER_LIB_DIRS} )
target_link_libraries(plummx plumlib ${SOLVER_LIB} )

#target_link_directories(plumseq PUBLIC ${SOLVER_LIB_DIRS} )
#target_link_libraries(plumseq plumlib gurobi_g++5.2 gurobi81 )

target_link_libraries(plumsp plumlib )
target_link_libraries(plumchk plumlib )
target_link_libraries(plumcmp plumlib )
target_link_libraries(plummerge plumlib )

#------------------------------------------------------------------------------#
# Unit Tests

# Option to enable/disable unit tests (enabled by default)
option(BUILD_TESTS "Build unit tests" ON)

if(BUILD_TESTS)
  message(STATUS "Building unit tests")
  add_subdirectory(unit_tests)
endif()

#------------------------------------------------------------------------------#
# Install

install(
  TARGETS plumsp plumchk plumcmp plum plummerge plummx
  RUNTIME DESTINATION bin
  )

install(CODE "MESSAGE(\"\n\
Plum was successfully installed in '${CMAKE_INSTALL_PREFIX}'.\")")

message(STATUS "Plum will be installed in '${CMAKE_INSTALL_PREFIX}'.")
