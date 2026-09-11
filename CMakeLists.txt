cmake_minimum_required(VERSION 3.21)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
    message(FATAL_ERROR "In-source builds are not allowed. Use: cmake --preset windows-msvc-debug")
endif()

project(MalloySim
    VERSION 0.1.0
    DESCRIPTION "MalloySim: an all-in-one science physics simulation in C++20"
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
        # The MSVC branch sets /fp:precise. GCC and Clang default to
        # -ffp-contract=fast, which lets the compiler fuse a*b+c into an FMA and
        # silently change results. docs/04 forbids fast-math, so match MSVC.
        -ffp-contract=off
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

# malloy_time (M3): header-only fixed timestep + tick counter.
add_library(malloy_time INTERFACE)
add_library(malloy::time ALIAS malloy_time)
target_include_directories(malloy_time INTERFACE include)
target_link_libraries(malloy_time INTERFACE malloy::math malloy::project_options)

# malloy_sim_core (M3): tiny shared simulation vocabulary (compiled STATIC).
add_library(malloy_sim_core STATIC
    src/sim_core/sim_core.cpp
)
add_library(malloy::sim_core ALIAS malloy_sim_core)
target_include_directories(malloy_sim_core PUBLIC include)
target_link_libraries(malloy_sim_core
    PUBLIC
        malloy::math
        malloy::project_options
    PRIVATE
        malloy::project_warnings
)

# malloy_nbody (M4): concrete 2D N-body gravity simulation (compiled STATIC).
add_library(malloy_nbody STATIC
    src/nbody/nbody.cpp
    src/nbody/nbody3d.cpp
)
add_library(malloy::nbody ALIAS malloy_nbody)
target_include_directories(malloy_nbody PUBLIC include)
target_link_libraries(malloy_nbody
    PUBLIC
        malloy::math
        malloy::sim_core
        malloy::time
        malloy::project_options
    PRIVATE
        malloy::project_warnings
)

# malloy_scenario (M7): minimal text scenario/config loading (compiled STATIC).
add_library(malloy_scenario STATIC
    src/scenario/scenario.cpp
)
add_library(malloy::scenario ALIAS malloy_scenario)
target_include_directories(malloy_scenario PUBLIC include)
target_link_libraries(malloy_scenario
    PUBLIC
        malloy::charges
        malloy::collide
        malloy::nbody
        malloy::particles
        malloy::rigid
        malloy::springs
        malloy::sim_core
        malloy::math
        malloy::project_options
    PRIVATE
        malloy::project_warnings
)

# malloy_ascii (M8): simple 2D ASCII debug visualization (compiled STATIC).
add_library(malloy_ascii STATIC
    src/ascii/grid.cpp
)
add_library(malloy::ascii ALIAS malloy_ascii)
target_include_directories(malloy_ascii PUBLIC include)
target_link_libraries(malloy_ascii
    PUBLIC
        malloy::math
        malloy::project_options
    PRIVATE
        malloy::project_warnings
)

# malloy_collide (M9): 2D collision primitives and contacts (compiled STATIC).
add_library(malloy_collide STATIC
    src/collide/collide.cpp
    src/collide/collide3d.cpp
)
add_library(malloy::collide ALIAS malloy_collide)
target_include_directories(malloy_collide PUBLIC include)
target_link_libraries(malloy_collide
    PUBLIC
        malloy::math
        malloy::project_options
    PRIVATE
        malloy::project_warnings
)

# malloy_particles (M10): 2D colliding particles (compiled STATIC).
add_library(malloy_particles STATIC
    src/particles/particles.cpp
)
add_library(malloy::particles ALIAS malloy_particles)
target_include_directories(malloy_particles PUBLIC include)
target_link_libraries(malloy_particles
    PUBLIC
        malloy::collide
        malloy::math
        malloy::sim_core
        malloy::time
        malloy::project_options
    PRIVATE
        malloy::project_warnings
)

# malloy_rigid (M11): 2D rigid bodies, pose and mass distribution (STATIC).
add_library(malloy_rigid STATIC
    src/rigid/rigid.cpp
    src/rigid/rigid3d.cpp
    src/rigid/mass_properties3d.cpp
)
add_library(malloy::rigid ALIAS malloy_rigid)
target_include_directories(malloy_rigid PUBLIC include)
target_link_libraries(malloy_rigid
    PUBLIC
        malloy::collide
        malloy::math
        malloy::sim_core
        malloy::time
        malloy::project_options
    PRIVATE
        malloy::project_warnings
)

# malloy_springs (M13): spring networks and deterministic force accumulation.
add_library(malloy_springs STATIC
    src/springs/springs.cpp
)
add_library(malloy::springs ALIAS malloy_springs)
target_include_directories(malloy_springs PUBLIC include)
target_link_libraries(malloy_springs
    PUBLIC
        malloy::math
        malloy::sim_core
        malloy::time
        malloy::project_options
    PRIVATE
        malloy::project_warnings
)

# malloy_charges (M18): charged particles in electric and magnetic fields,
# in two dimensions and (M29) three.
add_library(malloy_charges STATIC
    src/charges/charges.cpp
    src/charges/charges3d.cpp
)
add_library(malloy::charges ALIAS malloy_charges)
target_include_directories(malloy_charges PUBLIC include)
target_link_libraries(malloy_charges
    PUBLIC
        malloy::math
        malloy::sim_core
        malloy::time
        malloy::project_options
    PRIVATE
        malloy::project_warnings
)

if(MALLOY_BUILD_APPS)
    add_executable(malloy_nbody_terminal
        apps/nbody_terminal/main.cpp
    )

    target_link_libraries(malloy_nbody_terminal
        PRIVATE
            malloy::ascii
            malloy::charges
            malloy::nbody
            malloy::scenario
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

    add_executable(malloy_quat_tests
        tests/math/quat_tests.cpp
    )

    target_include_directories(malloy_quat_tests
        PRIVATE
            tests
    )

    target_link_libraries(malloy_quat_tests
        PRIVATE
            malloy::math
            malloy::project_options
            malloy::project_warnings
    )

    add_test(NAME malloy_quat_tests COMMAND malloy_quat_tests)

    add_executable(malloy_time_tests
        tests/time/fixed_step_tests.cpp
    )

    target_include_directories(malloy_time_tests
        PRIVATE
            tests
    )

    target_link_libraries(malloy_time_tests
        PRIVATE
            malloy::time
            malloy::project_options
            malloy::project_warnings
    )

    add_test(NAME malloy_time_tests COMMAND malloy_time_tests)

    add_executable(malloy_sim_core_tests
        tests/sim_core/sim_core_tests.cpp
    )

    target_include_directories(malloy_sim_core_tests
        PRIVATE
            tests
    )

    target_link_libraries(malloy_sim_core_tests
        PRIVATE
            malloy::sim_core
            malloy::project_options
            malloy::project_warnings
    )

    add_test(NAME malloy_sim_core_tests COMMAND malloy_sim_core_tests)

    add_executable(malloy_nbody_tests
        tests/nbody/nbody_world_tests.cpp
    )

    target_include_directories(malloy_nbody_tests
        PRIVATE
            tests
    )

    target_link_libraries(malloy_nbody_tests
        PRIVATE
            malloy::nbody
            malloy::project_options
            malloy::project_warnings
    )

    add_test(NAME malloy_nbody_tests COMMAND malloy_nbody_tests)

    add_executable(malloy_scenario_tests
        tests/scenario/scenario_tests.cpp
    )

    target_include_directories(malloy_scenario_tests
        PRIVATE
            tests
    )

    target_link_libraries(malloy_scenario_tests
        PRIVATE
            malloy::scenario
            malloy::project_options
            malloy::project_warnings
    )

    # The shipped templates are a deliverable (CLAUDE.md rule 16), so the
    # scenario tests parse every one of them rather than trusting them by sight.
    target_compile_definitions(malloy_scenario_tests
        PRIVATE
            MALLOY_SCENARIO_DIR="${CMAKE_CURRENT_SOURCE_DIR}/scenarios"
    )

    add_test(NAME malloy_scenario_tests COMMAND malloy_scenario_tests)

    add_executable(malloy_ascii_tests
        tests/ascii/grid_tests.cpp
    )

    target_include_directories(malloy_ascii_tests
        PRIVATE
            tests
    )

    target_link_libraries(malloy_ascii_tests
        PRIVATE
            malloy::ascii
            malloy::project_options
            malloy::project_warnings
    )

    add_test(NAME malloy_ascii_tests COMMAND malloy_ascii_tests)

    add_executable(malloy_collide_tests
        tests/collide/collide_tests.cpp
    )

    target_include_directories(malloy_collide_tests
        PRIVATE
            tests
    )

    target_link_libraries(malloy_collide_tests
        PRIVATE
            malloy::collide
            malloy::project_options
            malloy::project_warnings
    )

    add_test(NAME malloy_collide_tests COMMAND malloy_collide_tests)

    add_executable(malloy_particles_tests
        tests/particles/particle_world_tests.cpp
    )

    target_include_directories(malloy_particles_tests
        PRIVATE
            tests
    )

    target_link_libraries(malloy_particles_tests
        PRIVATE
            malloy::particles
            malloy::project_options
            malloy::project_warnings
    )

    add_test(NAME malloy_particles_tests COMMAND malloy_particles_tests)

    add_executable(malloy_rigid_tests
        tests/rigid/rigid_world_tests.cpp
    )

    target_include_directories(malloy_rigid_tests
        PRIVATE
            tests
    )

    target_link_libraries(malloy_rigid_tests
        PRIVATE
            malloy::rigid
            malloy::project_options
            malloy::project_warnings
    )

    add_test(NAME malloy_rigid_tests COMMAND malloy_rigid_tests)

    add_executable(malloy_springs_tests
        tests/springs/spring_world_tests.cpp
    )

    target_include_directories(malloy_springs_tests
        PRIVATE
            tests
    )

    target_link_libraries(malloy_springs_tests
        PRIVATE
            malloy::springs
            malloy::project_options
            malloy::project_warnings
    )

    add_test(NAME malloy_springs_tests COMMAND malloy_springs_tests)

    add_executable(malloy_charges_tests
        tests/charges/charge_world_tests.cpp
    )

    target_include_directories(malloy_charges_tests
        PRIVATE
            tests
    )

    target_link_libraries(malloy_charges_tests
        PRIVATE
            malloy::charges
            malloy::project_options
            malloy::project_warnings
    )

    add_test(NAME malloy_charges_tests COMMAND malloy_charges_tests)

    # Issue #13: the module graph is enforced at link time and nowhere else,
    # because every module puts the whole include/ tree on its consumers'
    # search path. This target reads the build files and the sources and checks
    # they agree, so it needs the repository root.
    add_executable(malloy_layering_tests
        tests/layering/layering_tests.cpp
    )

    target_include_directories(malloy_layering_tests
        PRIVATE
            tests
    )

    target_compile_definitions(malloy_layering_tests
        PRIVATE
            MALLOY_REPO_ROOT="${CMAKE_CURRENT_SOURCE_DIR}"
    )

    target_link_libraries(malloy_layering_tests
        PRIVATE
            malloy::project_options
            malloy::project_warnings
    )

    add_test(NAME malloy_layering_tests COMMAND malloy_layering_tests)

    # The number of test executables is stated in README.md and
    # docs/00_START_HERE.md, and such counts go stale. This is the only
    # place that knows the real number, so it is handed to the scenario
    # tests, which check the documents against it.
    get_property(malloy_registered_tests DIRECTORY PROPERTY TESTS)
    list(LENGTH malloy_registered_tests malloy_test_count)
    target_compile_definitions(malloy_scenario_tests
        PRIVATE
            MALLOY_TEST_EXECUTABLE_COUNT=${malloy_test_count}
    )
endif()
