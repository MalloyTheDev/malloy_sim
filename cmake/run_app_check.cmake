# Run the terminal app once on one scenario and check BOTH how it exited and what
# it printed (issue #65). CTest's PASS_REGULAR_EXPRESSION on its own ignores the
# exit code, so an app that printed a refusal yet returned success, the failure
# issue #31 guards against, would still pass.
#
#   cmake -DAPP=<exe> -DSCENARIO=<file> -DEXPECT_EXIT=<code>
#         -DEXPECT_OUTPUT=<regex> -P cmake/run_app_check.cmake
#
# The regex is matched against stdout and stderr together. The check fails, with
# the app's whole output shown, on any other exit code, on a crash (which reports
# no exit code at all, so it can never equal the expected one), or when nothing
# matches.
cmake_minimum_required(VERSION 3.21)

foreach(required IN ITEMS APP SCENARIO EXPECT_EXIT EXPECT_OUTPUT)
    if(NOT DEFINED ${required})
        message(FATAL_ERROR "run_app_check: ${required} is not set")
    endif()
endforeach()

execute_process(
    COMMAND "${APP}" "${SCENARIO}"
    RESULT_VARIABLE exit_code
    OUTPUT_VARIABLE output
    ERROR_VARIABLE output)

set(problem "")
if(NOT "${exit_code}" STREQUAL "${EXPECT_EXIT}")
    set(problem "exited with '${exit_code}', expected ${EXPECT_EXIT}")
elseif(NOT "${output}" MATCHES "${EXPECT_OUTPUT}")
    set(problem "printed nothing matching '${EXPECT_OUTPUT}'")
endif()

if(problem)
    message("${output}")
    message(FATAL_ERROR "run_app_check: ${SCENARIO}: the app ${problem}")
endif()
