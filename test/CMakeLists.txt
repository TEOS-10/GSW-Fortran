# Symlink the test data
list( APPEND gsw_test_data
    gsw_data_v3_0.nc
)

foreach( FILENAME ${gsw_test_data} )
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach( FILENAME )

# Get a list of test programs
list( APPEND gsw_test_files
    gsw_poly_check.f90
    gsw_check_functions.f90
)

# Create test executables and add tests
foreach( test_file ${gsw_test_files} )
    string( REPLACE ".f90" "" test_name ${test_file} )
    add_executable( ${test_name} ${test_file} )
    target_link_libraries( ${test_name} gsw)
    add_test( NAME ${test_name} COMMAND ${test_name} )
endforeach( test_file ${gsw_test_files} )
