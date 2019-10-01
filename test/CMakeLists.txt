# Symlink the test data
list( APPEND test_data
    gsw_data_v3_0.nc
)

foreach( FILENAME ${test_data} )
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
                     ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
                     ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach()

# Get a list of test sources
list( APPEND test_srcs
    gsw_poly_check.f90
    gsw_check_functions.f90
)

# Create executable for each test program, and add test
foreach( src ${test_srcs} )

    string( REPLACE ".f90" ".x" exec_name ${src} )
    add_executable( ${exec_name} ${src} )
    add_dependencies( ${exec_name} ${LIB_TARGETS} )
    target_link_libraries( ${exec_name} ${LIB_TARGETS} )
    target_include_directories( ${exec_name}
        PUBLIC $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/${INCLUDE_DIR}> )

    string( REPLACE ".f90" ""   test_name ${src} )
    add_test( ${test_name} COMMAND ${exec_name} )

endforeach()
