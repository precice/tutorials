add_subdirectory(macro-heat)
add_subdirectory(micro-heat)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/precice-config-full-sim.xml
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/precice-config-with-adaptivity.xml
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})