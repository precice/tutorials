add_executable(macro_dumux EXCLUDE_FROM_ALL main.cc)

target_link_libraries(macro_dumux PRIVATE dumux-precice)

add_input_file_links()
add_precice_file_links()

# Dummy test is needed to build porous media solver. The default
# `dumux_add_test` (and the underlying `dune_add_test`) only accept
# one target executable.
dumux_add_test(NAME dummy_test_ff_pm_part_flow_over_square_2d
              TARGET pm_flow_over_square_2d
              LABELS freeflow stokes precice darcy
              TIMEOUT 30
              COMMAND ${CMAKE_SOURCE_DIR}/test/return-test-passed.sh
)

dumux_add_test(NAME test_ff_pm_part_flow_over_square_2d_si_st_first
              TARGET ff_flow_over_square_2d
              TARGET pm_flow_over_square_2d
              LABELS freeflow stokes precice darcy
              TIMEOUT 30
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/test/run-test.py
              CMD_ARGS
              --solver-one ${CMAKE_CURRENT_BINARY_DIR}/ff_flow_over_square_2d
              --solver-two ${CMAKE_CURRENT_BINARY_DIR}/pm_flow_over_square_2d
              --files ${CMAKE_SOURCE_DIR}/test/reference-solutions/flow-over-square-2d/test_stokes_first_si_stokes.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/ff-pm-2d-precice-stokes-00005.vtu
                      ${CMAKE_SOURCE_DIR}/test/reference-solutions/flow-over-square-2d/test_stokes_first_si_darcy.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/ff-pm-2d-precice-darcy-00005.vtu
              --dumux-param-file ${CMAKE_CURRENT_BINARY_DIR}/params-stokesdarcy.input
              --precice-config-file ${CMAKE_CURRENT_BINARY_DIR}/precice-config-si-free-flow-first.xml
              --precice-iteration-files ${CMAKE_SOURCE_DIR}/test/reference-solutions/flow-over-square-2d/test_stokes_first_si_coupliter_stokes.log ${CMAKE_CURRENT_BINARY_DIR}/precice-FreeFlow-iterations.log
              --case-name "flow-over-box-2d-si-stokes-first"
              --relative 5e-6
              --zeroThreshold {"velocity_liq \(m/s\)":1e-14,"p":1e-12}
)
set_tests_properties(test_ff_pm_part_flow_over_square_2d_si_st_first PROPERTIES
    ENVIRONMENT PYTHONPATH=${CMAKE_SOURCE_DIR}/test:$ENV{PYTHONPATH})

dumux_add_test(NAME test_ff_pm_part_flow_over_square_2d_si_st_second
    TARGET pm_flow_over_square_2d
    TARGET ff_flow_over_square_2d
    LABELS freeflow darcy precice darcy
    TIMEOUT 30
    CMAKE_GUARD HAVE_UMFPACK
    COMMAND ${CMAKE_SOURCE_DIR}/test/run-test.py
    CMD_ARGS
            --solver-one ${CMAKE_CURRENT_BINARY_DIR}/ff_flow_over_square_2d
            --solver-two ${CMAKE_CURRENT_BINARY_DIR}/pm_flow_over_square_2d
            --files ${CMAKE_SOURCE_DIR}/test/reference-solutions/flow-over-square-2d/test_stokes_second_si_stokes.vtu
              ${CMAKE_CURRENT_BINARY_DIR}/ff-pm-2d-precice-stokes-00006.vtu
              ${CMAKE_SOURCE_DIR}/test/reference-solutions/flow-over-square-2d/test_stokes_second_si_darcy.vtu
              ${CMAKE_CURRENT_BINARY_DIR}/ff-pm-2d-precice-darcy-00006.vtu
            --dumux-param-file ${CMAKE_CURRENT_BINARY_DIR}/params-stokesdarcy.input
            --precice-config-file ${CMAKE_CURRENT_BINARY_DIR}/precice-config-si-free-flow-second.xml
            --precice-iteration-files ${CMAKE_SOURCE_DIR}/test/reference-solutions/flow-over-square-2d/test_stokes_second_si_coupliter_darcy.log ${CMAKE_CURRENT_BINARY_DIR}/precice-Darcy-iterations.log
    --case-name "flow-over-box-2d-si-stokes-second"
    --relative 5e-6
    --zeroThreshold {"velocity_liq \(m/s\)":1e-14,"p":1e-12}
)
set_tests_properties(test_ff_pm_part_flow_over_square_2d_si_st_second PROPERTIES
  ENVIRONMENT PYTHONPATH=${CMAKE_SOURCE_DIR}/test:$ENV{PYTHONPATH})

dumux_add_test(NAME test_ff_pm_part_flow_over_square_2d_pi_st
              TARGET ff_flow_over_square_2d
              LABELS freeflow stokes precice darcy
              TIMEOUT 30
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/test/run-test.py
              CMD_ARGS
              --solver-one ${CMAKE_CURRENT_BINARY_DIR}/ff_flow_over_square_2d
              --solver-two ${CMAKE_CURRENT_BINARY_DIR}/pm_flow_over_square_2d
              --files ${CMAKE_SOURCE_DIR}/test/reference-solutions/flow-over-square-2d/test_stokes_pi_stokes.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/ff-pm-2d-precice-stokes-00007.vtu
                      ${CMAKE_SOURCE_DIR}/test/reference-solutions/flow-over-square-2d/test_stokes_pi_darcy.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/ff-pm-2d-precice-darcy-00007.vtu
              --dumux-param-file ${CMAKE_CURRENT_BINARY_DIR}/params-stokesdarcy.input
              --precice-config-file ${CMAKE_CURRENT_BINARY_DIR}/precice-config-pi.xml
              --precice-iteration-files ${CMAKE_SOURCE_DIR}/test/reference-solutions/flow-over-square-2d/test_stokes_pi_coupliter_stokes.log
              ${CMAKE_CURRENT_BINARY_DIR}/precice-FreeFlow-iterations.log
              --case-name "flow-over-box-2d-pi"
              --relative 5e-6
              --zeroThreshold {"velocity_liq \(m/s\)":1e-14,"p":1e-12}
)
# Adding DuMuX directory to Python search path such that we can use DuMuX'x fuzzy (inexact) testing of vtu files
set_tests_properties(test_ff_pm_part_flow_over_square_2d_pi_st PROPERTIES
    ENVIRONMENT PYTHONPATH=${CMAKE_SOURCE_DIR}/test:$ENV{PYTHONPATH})

### Setting test dependencies
set_property(TEST test_ff_pm_part_flow_over_square_2d_si_st_second APPEND PROPERTY DEPENDS test_ff_pm_part_flow_over_square_2d_si_st_first)
set_property(TEST test_ff_pm_part_flow_over_square_2d_pi_st APPEND PROPERTY DEPENDS test_ff_pm_part_flow_over_square_2d_si_st_second)

######################
# Navier-Stokes tests
######################

dumux_add_test(NAME test_ff_pm_part_flow_over_square_2d_si_ns_first
              TARGET ff_flow_over_square_2d
              LABELS freeflow navierstokes precice darcy
              TIMEOUT 60
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/test/run-test.py
              CMD_ARGS
              --solver-one ${CMAKE_CURRENT_BINARY_DIR}/ff_flow_over_square_2d
              --solver-two ${CMAKE_CURRENT_BINARY_DIR}/pm_flow_over_square_2d
              --files ${CMAKE_SOURCE_DIR}/test/reference-solutions/flow-over-square-2d/test_navierstokes_first_si_navierstokes.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/ff-pm-2d-precice-navierstokes-00006.vtu
                      ${CMAKE_SOURCE_DIR}/test/reference-solutions/flow-over-square-2d/test_navierstokes_first_si_darcy.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/ff-pm-2d-precice-darcy-00006.vtu
              --dumux-param-file ${CMAKE_CURRENT_BINARY_DIR}/params-navierstokesdarcy.input
              --precice-config-file ${CMAKE_CURRENT_BINARY_DIR}/precice-config-si-free-flow-first.xml
              --precice-iteration-files ${CMAKE_SOURCE_DIR}/test/reference-solutions/flow-over-square-2d/test_navierstokes_first_si_coupliter_navierstokes.log ${CMAKE_CURRENT_BINARY_DIR}/precice-FreeFlow-iterations.log
              --case-name "flow-over-box-2d-si-navierstokes-first"
              --relative 5e-6
              --zeroThreshold {"velocity_liq \(m/s\)":1e-14,"p":1e-12}
)
set_tests_properties(test_ff_pm_part_flow_over_square_2d_si_ns_first PROPERTIES
    ENVIRONMENT PYTHONPATH=${CMAKE_SOURCE_DIR}/test:$ENV{PYTHONPATH})

dumux_add_test(NAME test_ff_pm_part_flow_over_square_2d_si_ns_second
    TARGET ff_flow_over_square_2d
    LABELS freeflow darcy precice darcy
    TIMEOUT 60
    CMAKE_GUARD HAVE_UMFPACK
    COMMAND ${CMAKE_SOURCE_DIR}/test/run-test.py
    CMD_ARGS
    --solver-one ${CMAKE_CURRENT_BINARY_DIR}/ff_flow_over_square_2d
    --solver-two ${CMAKE_CURRENT_BINARY_DIR}/pm_flow_over_square_2d
    --files ${CMAKE_SOURCE_DIR}/test/reference-solutions/flow-over-square-2d/test_navierstokes_second_si_navierstokes.vtu
              ${CMAKE_CURRENT_BINARY_DIR}/ff-pm-2d-precice-navierstokes-00007.vtu
              ${CMAKE_SOURCE_DIR}/test/reference-solutions/flow-over-square-2d/test_navierstokes_second_si_darcy.vtu
              ${CMAKE_CURRENT_BINARY_DIR}/ff-pm-2d-precice-darcy-00007.vtu
            --dumux-param-file ${CMAKE_CURRENT_BINARY_DIR}/params-navierstokesdarcy.input
            --precice-config-file ${CMAKE_CURRENT_BINARY_DIR}/precice-config-si-free-flow-second.xml
            --precice-iteration-files ${CMAKE_SOURCE_DIR}/test/reference-solutions/flow-over-square-2d/test_navierstokes_second_si_coupliter_darcy.log ${CMAKE_CURRENT_BINARY_DIR}/precice-Darcy-iterations.log
    --case-name "flow-over-box-2d-si-navierstokes-second"
    --relative 5e-6
    --zeroThreshold {"velocity_liq \(m/s\)":1e-14,"p":1e-12}
)
set_tests_properties(test_ff_pm_part_flow_over_square_2d_si_ns_second PROPERTIES
 ENVIRONMENT PYTHONPATH=${CMAKE_SOURCE_DIR}/test:$ENV{PYTHONPATH})

dumux_add_test(NAME test_ff_pm_part_flow_over_square_2d_pi_ns
              TARGET ff_flow_over_square_2d
              LABELS freeflow navierstokes precice darcy
              TIMEOUT 60
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/test/run-test.py
              CMD_ARGS
              --solver-one ${CMAKE_CURRENT_BINARY_DIR}/ff_flow_over_square_2d
              --solver-two ${CMAKE_CURRENT_BINARY_DIR}/pm_flow_over_square_2d
              --files ${CMAKE_SOURCE_DIR}/test/reference-solutions/flow-over-square-2d/test_navierstokes_pi_navierstokes.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/ff-pm-2d-precice-navierstokes-00011.vtu
                      ${CMAKE_SOURCE_DIR}/test/reference-solutions/flow-over-square-2d/test_navierstokes_pi_darcy.vtu
                      ${CMAKE_CURRENT_BINARY_DIR}/ff-pm-2d-precice-darcy-00011.vtu
              --dumux-param-file ${CMAKE_CURRENT_BINARY_DIR}/params-navierstokesdarcy.input
              --precice-config-file ${CMAKE_CURRENT_BINARY_DIR}/precice-config-pi.xml
              --precice-iteration-files ${CMAKE_SOURCE_DIR}/test/reference-solutions/flow-over-square-2d/test_navierstokes_pi_coupliter_navierstokes.log
              ${CMAKE_CURRENT_BINARY_DIR}/precice-FreeFlow-iterations.log
              --case-name "flow-over-box-2d-pi"
              --relative 5e-6
              --zeroThreshold {"velocity_liq \(m/s\)":1e-14,"p":1e-12}
)
# Adding DuMuX directory to Python search path such that we can use DuMuX'x fuzzy (inexact) testing of vtu files
set_tests_properties(test_ff_pm_part_flow_over_square_2d_pi_ns PROPERTIES
    ENVIRONMENT PYTHONPATH=${CMAKE_SOURCE_DIR}/test:$ENV{PYTHONPATH})

# Adding DuMuX directory to Python search path such that we can use DuMuX'x fuzzy (inexact) testing of vtu files
set_property(TEST test_ff_pm_part_flow_over_square_2d_si_ns_first APPEND PROPERTY DEPENDS test_ff_pm_part_flow_over_square_2d_pi_st)
set_property(TEST test_ff_pm_part_flow_over_square_2d_si_ns_second APPEND PROPERTY DEPENDS test_ff_pm_part_flow_over_square_2d_si_ns_first)
set_property(TEST test_ff_pm_part_flow_over_square_2d_pi_ns APPEND PROPERTY DEPENDS test_ff_pm_part_flow_over_square_2d_si_ns_second)
