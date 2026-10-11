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
# Off by default: the fuzzers need a fuzzer/sanitizer-capable toolchain (clang),
# not the default MSVC presets. See fuzz/parse_scenario_fuzzer.cpp.
option(MALLOY_BUILD_FUZZERS "Build MalloySim fuzz targets (requires clang)" OFF)

add_library(malloy_project_options INTERFACE)
add_library(malloy::project_options ALIAS malloy_project_options)
target_compile_features(malloy_project_options INTERFACE cxx_std_20)

# Issue #37: exploit mitigations for every MSVC target, in every configuration.
# They sit on project_options because every library links it PUBLIC, so every
# object file and every executable gets them. Control Flow Guard needs both the
# compile and the link flag; /CETCOMPAT marks the image shadow-stack compatible.
# Neither changes a result or costs measurable time.
#
# Two further hardenings were measured and declined, because the Release pair
# loop is the product (docs/06 has the timings). _MSVC_STL_HARDENING=1 costs
# about 20%, and Debug already runs the same checks through iterator
# debugging; to turn it on, add
#   target_compile_definitions(malloy_project_options INTERFACE _MSVC_STL_HARDENING=1)
# here. /sdl costs about 40%, and every warning it promotes is already an error
# under /W4 /WX below, so it would add only runtime checks; to turn it on, add
# /sdl to the target_compile_options line here.
if(MSVC)
    target_compile_options(malloy_project_options INTERFACE /guard:cf)
    target_link_options(malloy_project_options INTERFACE /guard:cf /CETCOMPAT)
endif()

add_library(malloy_project_warnings INTERFACE)
add_library(malloy::project_warnings ALIAS malloy_project_warnings)

if(MSVC)
    # /WX (issue #37): the tree builds warning-free at /W4, so a new warning
    # fails the build instead of accumulating unnoticed.
    target_compile_options(malloy_project_warnings INTERFACE
        /W4
        /WX
        /permissive-
        /fp:precise
    )
else()
    target_compile_options(malloy_project_warnings INTERFACE
        -Wall
        -Wextra
        -Wpedantic
        # As /WX above (issue #66): a warning fails the build under clang and
        # GCC as well, instead of accumulating unnoticed.
        -Werror
        # The MSVC branch sets /fp:precise. GCC and Clang default to
        # -ffp-contract=fast, which lets the compiler fuse a*b+c into an FMA and
        # silently change results. docs/04 forbids fast-math, so match MSVC.
        -ffp-contract=off
    )
endif()

# ----------------------------------------------------------------------------
# Libraries
# ----------------------------------------------------------------------------

# malloy_math (M2): header-only Real, Vec2, Vec3, Quat and Mat3 vocabulary.
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

# malloy_nbody (M4): concrete N-body gravity, 2D and (M19) 3D (compiled STATIC).
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

# malloy_collide (M9): 2D and 3D collision primitives and contacts (STATIC).
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

# malloy_particles (M10): colliding particles in 2D and (M32) 3D (STATIC).
add_library(malloy_particles STATIC
    src/particles/particles.cpp
    src/particles/particles3d.cpp
    src/particles/thermo.cpp
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

# malloy_rigid (M11): 2D and (M20) 3D rigid bodies, contacts, mass properties
# (STATIC).
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

# malloy_springs (M13): spring networks and deterministic force accumulation,
# in two dimensions and (M31) three.
add_library(malloy_springs STATIC
    src/springs/springs.cpp
    src/springs/springs3d.cpp
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
    # The manifest makes UTF-8 the process code page (Windows 10 1903 and
    # later), so the narrow argv and the narrow path the scenario parser opens
    # both carry UTF-8, and a path outside the system code page opens (issue
    # #33). CMake merges it into the linker's own manifest; toolchains that do
    # not embed manifests do not compile it.
    add_executable(malloy_nbody_terminal
        apps/nbody_terminal/main.cpp
        apps/nbody_terminal/nbody_terminal.manifest
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

# ----------------------------------------------------------------------------
# App checks: run the real terminal app on inputs it must refuse cleanly, with a
# clear message rather than a crash or a "successful" empty run (issues #24 and
# #31), on one it must open (issue #33), and on one shipped template of every
# scenario type (issue #65). They run the app rather than adding a test
# executable, so they are registered after the count above, which is unchanged.
# The refused inputs are written at configure time, so no fixture files live in
# the repository.
#
# Every check goes through cmake/run_app_check.cmake, which requires the exit
# code AND the output (issue #65). PASS_REGULAR_EXPRESSION alone makes CTest
# ignore the exit code, so a refusal printed with exit 0, the failure #31 is
# about, used to pass. A crash has no exit code, so it fails too.
# ----------------------------------------------------------------------------
if(MALLOY_BUILD_TESTS AND MALLOY_BUILD_APPS)
    function(malloy_add_app_check name scenario exit_code pattern)
        add_test(NAME ${name}
                 COMMAND ${CMAKE_COMMAND}
                     "-DAPP=$<TARGET_FILE:malloy_nbody_terminal>"
                     "-DSCENARIO=${scenario}"
                     "-DEXPECT_EXIT=${exit_code}"
                     "-DEXPECT_OUTPUT=${pattern}"
                     -P "${PROJECT_SOURCE_DIR}/cmake/run_app_check.cmake")
    endfunction()

    set(malloy_app_inputs "${CMAKE_CURRENT_BINARY_DIR}/app_inputs")
    file(WRITE "${malloy_app_inputs}/empty.scn" "")
    file(WRITE "${malloy_app_inputs}/comment_only.scn" "# nothing but a comment\n")
    file(WRITE "${malloy_app_inputs}/charges_empty.scn" "type charges\nsteps 1\n")
    file(WRITE "${malloy_app_inputs}/charges3d_empty.scn" "type charges3d\nsteps 1\n")
    foreach(input IN ITEMS empty comment_only charges_empty charges3d_empty)
        malloy_add_app_check(malloy_app_refuses_${input}
            "${malloy_app_inputs}/${input}.scn" 1 "defines nothing to simulate")
    endforeach()

    # Issue #33: a scenario whose path lies outside the ANSI code page must
    # open. The name holds a Greek psi (U+03C8, missing from the Western code
    # pages) and a rocket (U+1F680, outside every ANSI code page), so without
    # the app's UTF-8 manifest it reaches the parser as "?" and the open fails
    # on any system. Only a run that loaded the file and stepped it prints the
    # step 10 line.
    file(WRITE "${malloy_app_inputs}/two_body_ψ🚀.scn"
        "dt 0.001\nsteps 10\noutput_every 10\n"
        "body 1.0 0.0 0.0 0.0 0.0\nbody 0.000001 1.0 0.0 0.0 1.0\n")
    malloy_add_app_check(malloy_app_opens_unicode_path
        "${malloy_app_inputs}/two_body_ψ🚀.scn" 0 "step +10 +sep")

    # Issue #65: one shipped template per scenario type, run to the end. Each
    # must exit 0 and print its final report, the `step <n>` line for the last
    # multiple of output_every within steps. Both numbers are read from the
    # template here, so a template can change its length without this list
    # changing, and the type is checked against the one named, so the list
    # keeps covering every runner. The values on the line are not matched: the
    # templates' `# check` lines, which the scenario tests evaluate, pin those.
    function(malloy_add_template_check type template)
        set(scenario "${PROJECT_SOURCE_DIR}/scenarios/${template}.scn")
        set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${scenario}")
        file(STRINGS "${scenario}" keys REGEX "^(type|steps|output_every) ")
        set(declared nbody) # the parser's default type
        set(steps "")
        set(every "")
        foreach(key IN LISTS keys)
            if(key MATCHES "^type +([a-z0-9]+)")
                set(declared "${CMAKE_MATCH_1}")
            elseif(key MATCHES "^steps +([0-9]+)")
                set(steps "${CMAKE_MATCH_1}")
            elseif(key MATCHES "^output_every +([0-9]+)")
                set(every "${CMAKE_MATCH_1}")
            endif()
        endforeach()
        if(NOT declared STREQUAL type OR steps STREQUAL "" OR every STREQUAL "")
            message(FATAL_ERROR "scenarios/${template}.scn: an app check needs a "
                                "${type} template that states steps and output_every")
        endif()
        set(last 0) # output_every 0 prints only the step-0 report
        if(every GREATER 0)
            math(EXPR last "${steps} - ${steps} % ${every}")
        endif()
        malloy_add_app_check(malloy_app_runs_${type} "${scenario}" 0 "step +${last}[^0-9]")
    endfunction()

    malloy_add_template_check(nbody three_body_triangle)
    malloy_add_template_check(nbody3d inclined_orbits)
    malloy_add_template_check(particles heat_conduction)
    malloy_add_template_check(particles3d particle_collision3d)
    malloy_add_template_check(rigid tumbling_impact)
    malloy_add_template_check(rigid3d box_sphere_collision)
    malloy_add_template_check(springs spring_chain)
    malloy_add_template_check(springs3d spring_tetrahedron)
    malloy_add_template_check(charges cyclotron)
    malloy_add_template_check(charges3d magnetic_helix)

    # The app's own rigid3d E column includes the applied force's potential
    # (issue #47), which no # check line reads: the scenario tests evaluate the
    # quantities through the library, not the app. Pinned here on the final
    # report of the force template, whose energy is conserved but for the
    # integrator's derived drift; without the potential it would read 6.28e+01.
    malloy_add_app_check(malloy_app_rigid3d_energy_includes_force
        "${PROJECT_SOURCE_DIR}/scenarios/applied_force.scn" 0
        "step +2000 +E +3.02083333e-01")
endif()

# ----------------------------------------------------------------------------
# Fuzzers (off by default; not registered as tests, so the count above is
# unchanged). Coverage-guided robustness for the scenario file, the one
# untrusted input (docs/06, Fuzzing): one target fuzzes the parser, the
# other the whole app behind it (issue #30). They need clang++, and the default
# MSVC presets never build them. The recipe is in fuzz/parse_scenario_fuzzer.cpp.
# ----------------------------------------------------------------------------
if(MALLOY_BUILD_FUZZERS)
    # The GNU-style driver, not clang-cl: CMake links through clang++, which is
    # what makes -fsanitize=fuzzer below pull in libFuzzer and its main().
    if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR
       CMAKE_CXX_COMPILER_FRONTEND_VARIANT STREQUAL "MSVC")
        message(FATAL_ERROR
            "MALLOY_BUILD_FUZZERS needs clang++ (the GNU-style clang driver), "
            "not MSVC or clang-cl. See fuzz/parse_scenario_fuzzer.cpp for the "
            "recipe.")
    endif()

    # Issue #27: instrument EVERY target, not only the fuzzers' entry points.
    # libFuzzer steers by the coverage counters compiled into the code, and
    # the sanitizers see only code compiled with them, so the libraries need
    # the flags as much as the harnesses do. Every target links
    # malloy_project_options, and an INTERFACE target's usage requirements are
    # read when the build is generated, so adding them here, after the targets
    # are defined, still reaches all of them. -fno-sanitize-recover: by default
    # UBSan prints a report and carries on, so libFuzzer would never see it.
    target_compile_options(malloy_project_options INTERFACE
        -g -O1 -fno-omit-frame-pointer
        -fsanitize=fuzzer-no-link,address,undefined
        -fno-sanitize-recover=undefined)
    target_link_options(malloy_project_options INTERFACE
        -fsanitize=address,undefined)

    add_executable(malloy_parse_scenario_fuzzer fuzz/parse_scenario_fuzzer.cpp)
    target_link_libraries(malloy_parse_scenario_fuzzer
        PRIVATE
            malloy::scenario
            malloy::project_warnings
    )
    target_link_options(malloy_parse_scenario_fuzzer PRIVATE -fsanitize=fuzzer)

    # Issue #30: the app's whole path after parsing. The harness compiles
    # apps/nbody_terminal/main.cpp in itself, so it needs what the app needs:
    # the ASCII view, and the scenario loader, which brings every domain.
    add_executable(malloy_app_fuzzer fuzz/app_fuzzer.cpp)
    target_link_libraries(malloy_app_fuzzer
        PRIVATE
            malloy::ascii
            malloy::scenario
            malloy::project_warnings
    )
    target_link_options(malloy_app_fuzzer PRIVATE -fsanitize=fuzzer)

    # Windows (clang++ targeting the MSVC ABI). libFuzzer's prebuilt runtime is
    # compiled against the static release CRT and without ASan's container
    # annotations, and the MSVC STL refuses to link objects that disagree with
    # it on the runtime library, the iterator debug level or the annotations
    # (/failifmismatch). So every target uses the static release CRT, in every
    # build type (which also keeps _DEBUG and its iterator level out of Debug),
    # and the STL annotations are off. The cost: ASan cannot see a vector or
    # string read past size() but inside capacity.
    if(CMAKE_CXX_SIMULATE_ID STREQUAL "MSVC")
        target_compile_definitions(malloy_project_options INTERFACE
            _DISABLE_STRING_ANNOTATION=1
            _DISABLE_VECTOR_ANNOTATION=1
            _DISABLE_OPTIONAL_ANNOTATION=1)
        get_property(malloy_fuzz_build_targets DIRECTORY PROPERTY BUILDSYSTEM_TARGETS)
        foreach(malloy_fuzz_build_target IN LISTS malloy_fuzz_build_targets)
            get_target_property(malloy_fuzz_target_type ${malloy_fuzz_build_target} TYPE)
            if(NOT malloy_fuzz_target_type STREQUAL "INTERFACE_LIBRARY")
                set_property(TARGET ${malloy_fuzz_build_target}
                             PROPERTY MSVC_RUNTIME_LIBRARY MultiThreaded)
            endif()
        endforeach()
    endif()
endif()
