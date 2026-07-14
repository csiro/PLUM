# Unit Tests for PLUM using Google Test
cmake_minimum_required(VERSION 3.14)

# GoogleTest requires at least C++14
set(CMAKE_CXX_STANDARD 20)

include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        v1.14.0
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# Enable testing
enable_testing()

# Include directories
include_directories(
  ${PROJECT_INCLUDE_DIR}
  ${LIME_INCLUDE_DIRS}
  ${SOLVER_INCLUDE_DIRS}
)

# Test executables
add_executable(metabolite_test
  metabolite_test.cpp
)

add_executable(reaction_test
  reaction_test.cpp
)

add_executable(scenario_test
  scenario_test.cpp
)

add_executable(solution_test
  solution_test.cpp
)

add_executable(multisol_test
  multisol_test.cpp
)

add_executable(params_test
  params_test.cpp
)

# Link against gtest, gtest_main, and plumlib
target_link_libraries(metabolite_test
  GTest::gtest_main
  plumlib
)

target_link_libraries(reaction_test
  GTest::gtest_main
  plumlib
)

target_link_libraries(scenario_test
  GTest::gtest_main
  plumlib
)

target_link_libraries(solution_test
  GTest::gtest_main
  plumlib
)

target_link_libraries(multisol_test
  GTest::gtest_main
  plumlib
)

target_link_libraries(params_test
  GTest::gtest_main
  plumlib
)

# Discover tests
include(GoogleTest)
gtest_discover_tests(metabolite_test)
gtest_discover_tests(reaction_test)
gtest_discover_tests(scenario_test)
gtest_discover_tests(solution_test)
gtest_discover_tests(multisol_test)
gtest_discover_tests(params_test)
