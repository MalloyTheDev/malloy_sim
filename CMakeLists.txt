cmake_minimum_required(VERSION 3.21)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
    message(FATAL_ERROR "In-source builds are not allowed. Use: cmake --preset windows-msvc-debug")
endif()

project(MalloySim
    VERSION 0.1.0
    DESCRIPTION "MalloySim: a simulation-first C++20 engine project"
    LANGUAGES CXX
)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(MALLOY_BUILD_APPS "Build MalloySim applications" ON)
option(MALLOY_BUILD_TESTS "Build MalloySim tests" ON)

add_library(malloy_project_options INTERFACE)
add_library(malloy::project_options ALIAS malloy_project_options)
target_compile_features(malloy_project_options INTERFACE cxx_std_20)

add_library(malloy_project_warnings INTERFACE)
add_library(malloy::project_warnings ALIAS malloy_project_warnings)

if(MSVC)
    target_compile_options(malloy_project_warnings INTERFACE
        /W4
        /permissive-
        /fp:precise
    )
else()
    target_compile_options(malloy_project_warnings INTERFACE
        -Wall
        -Wextra
        -Wpedantic
    )
endif()

# ----------------------------------------------------------------------------
# Libraries
# ----------------------------------------------------------------------------

# malloy_math (M2): header-only Real + Vec2 vocabulary.
add_library(malloy_math INTERFACE)
add_library(malloy::math ALIAS malloy_math)
target_include_directories(malloy_math INTERFACE include)
target_link_libraries(malloy_math INTERFACE malloy::project_options)

if(MALLOY_BUILD_APPS)
    add_executable(malloy_nbody_terminal
        apps/nbody_terminal/main.cpp
    )

    target_link_libraries(malloy_nbody_terminal
        PRIVATE
            malloy::project_options
            malloy::project_warnings
    )
endif()

if(MALLOY_BUILD_TESTS)
    enable_testing()

    add_executable(malloy_smoke_tests
        tests/smoke/smoke_tests.cpp
    )

    target_include_directories(malloy_smoke_tests
        PRIVATE
            tests
    )

    target_link_libraries(malloy_smoke_tests
        PRIVATE
            malloy::project_options
            malloy::project_warnings
    )

    add_test(NAME malloy_smoke_tests COMMAND malloy_smoke_tests)

    add_executable(malloy_math_tests
        tests/math/vec2_tests.cpp
    )

    target_include_directories(malloy_math_tests
        PRIVATE
            tests
    )

    target_link_libraries(malloy_math_tests
        PRIVATE
            malloy::math
            malloy::project_options
            malloy::project_warnings
    )

    add_test(NAME malloy_math_tests COMMAND malloy_math_tests)
endif()
