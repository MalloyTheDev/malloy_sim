# The viewer's app checks (M41, issue #80), the viewer's sibling of
# cmake/run_app_check.cmake, for the two things that script cannot express:
# running the viewer with flags (or with no argument at all), and comparing its
# report lines with the terminal app's.
#
#   cmake -DMODE=run -DVIEWER=<exe> [-DSCENARIO=<file>] [-DARGS=a|b|c]
#         -DEXPECT_EXIT=<code> -DEXPECT_OUTPUT=<regex> [-DEXPECT_LAST_LINE=<regex>]
#         -P viewer_check.cmake
#
#     Runs the viewer once and requires BOTH the exit code and a match in its
#     stdout and stderr together, as run_app_check does (issue #65), and, when
#     given, a match for the last line it printed. ARGS are extra arguments
#     after the scenario, separated by '|'. Without SCENARIO the viewer gets no
#     scenario argument.
#
#   cmake -DMODE=same -DVIEWER=<exe> -DTERMINAL=<exe> -DSCENARIO=<file>
#         -DFRAME_TIMES=<s>|<s>|... -P viewer_check.cmake
#
#     Runs the terminal app on the scenario, then the viewer with --headless at
#     each frame time, and requires every viewer run to print exactly the
#     terminal's report lines, every one of them, the last tick's included, byte
#     for byte. That is "the physics is untouched by the viewer" and "the frame
#     rate changes nothing", both at tolerance zero.
cmake_minimum_required(VERSION 3.21)

function(viewer_check_fail output problem)
    message("${output}")
    message(FATAL_ERROR "viewer_check: ${problem}")
endfunction()

# The report lines of an output: every line that starts with "step <n> ".
function(viewer_check_report_lines output result)
    string(REPLACE "\r" "" output "${output}")
    string(REGEX MATCHALL "\nstep +[0-9]+ [^\n]*" lines "\n${output}")
    list(TRANSFORM lines REPLACE "^\n" "")
    set(${result} "${lines}" PARENT_SCOPE)
endfunction()

if(NOT DEFINED MODE OR NOT DEFINED VIEWER)
    message(FATAL_ERROR "viewer_check: MODE and VIEWER must be set")
endif()

if(MODE STREQUAL "run")
    foreach(required IN ITEMS EXPECT_EXIT EXPECT_OUTPUT)
        if(NOT DEFINED ${required})
            message(FATAL_ERROR "viewer_check: ${required} is not set")
        endif()
    endforeach()
    set(command "${VIEWER}")
    if(DEFINED SCENARIO)
        list(APPEND command "${SCENARIO}")
    endif()
    if(DEFINED ARGS AND NOT ARGS STREQUAL "")
        string(REPLACE "|" ";" extra "${ARGS}")
        list(APPEND command ${extra})
    endif()
    execute_process(
        COMMAND ${command}
        RESULT_VARIABLE exit_code
        OUTPUT_VARIABLE output
        ERROR_VARIABLE output)
    if(NOT "${exit_code}" STREQUAL "${EXPECT_EXIT}")
        viewer_check_fail("${output}" "the viewer exited with '${exit_code}', expected ${EXPECT_EXIT}")
    endif()
    if(NOT "${output}" MATCHES "${EXPECT_OUTPUT}")
        viewer_check_fail("${output}" "the viewer printed nothing matching '${EXPECT_OUTPUT}'")
    endif()
    # Optionally, what the LAST line printed must match: where a run stopped.
    if(DEFINED EXPECT_LAST_LINE)
        string(REPLACE "\r" "" trimmed "${output}")
        string(STRIP "${trimmed}" trimmed)
        string(FIND "${trimmed}" "\n" newline REVERSE)
        math(EXPR start "${newline} + 1")
        string(SUBSTRING "${trimmed}" ${start} -1 last_line)
        if(NOT "${last_line}" MATCHES "${EXPECT_LAST_LINE}")
            viewer_check_fail("${output}"
                "the viewer's last line '${last_line}' does not match '${EXPECT_LAST_LINE}'")
        endif()
    endif()
elseif(MODE STREQUAL "same")
    foreach(required IN ITEMS TERMINAL SCENARIO FRAME_TIMES)
        if(NOT DEFINED ${required})
            message(FATAL_ERROR "viewer_check: ${required} is not set")
        endif()
    endforeach()
    execute_process(
        COMMAND "${TERMINAL}" "${SCENARIO}"
        RESULT_VARIABLE exit_code
        OUTPUT_VARIABLE output
        ERROR_VARIABLE errors)
    if(NOT exit_code EQUAL 0)
        viewer_check_fail("${output}${errors}" "the terminal app exited with '${exit_code}'")
    endif()
    viewer_check_report_lines("${output}" expected)
    list(LENGTH expected count)
    if(count EQUAL 0)
        viewer_check_fail("${output}" "the terminal app printed no report line")
    endif()
    list(GET expected -1 last)

    string(REPLACE "|" ";" frame_times "${FRAME_TIMES}")
    foreach(frame_time IN LISTS frame_times)
        execute_process(
            COMMAND "${VIEWER}" "${SCENARIO}" --headless --frame-time "${frame_time}"
            RESULT_VARIABLE exit_code
            OUTPUT_VARIABLE output
            ERROR_VARIABLE errors)
        if(NOT exit_code EQUAL 0)
            viewer_check_fail("${output}${errors}"
                "the viewer exited with '${exit_code}' at --frame-time ${frame_time}")
        endif()
        viewer_check_report_lines("${output}" actual)
        if(NOT actual STREQUAL expected)
            list(LENGTH actual actual_count)
            viewer_check_fail("terminal app:\n${expected}\nviewer:\n${actual}"
                "at --frame-time ${frame_time} the viewer printed ${actual_count} report "
                "lines that are not the terminal app's ${count}")
        endif()
    endforeach()
    message("viewer_check: ${count} report lines identical to the terminal app's at "
            "--frame-time ${FRAME_TIMES}; the last: ${last}")
else()
    message(FATAL_ERROR "viewer_check: unknown MODE '${MODE}'")
endif()
