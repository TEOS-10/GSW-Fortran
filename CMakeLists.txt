cmake_minimum_required( VERSION 3.12 )

project( gsw VERSION 3.0.5 LANGUAGES Fortran )

## Ecbuild integration
find_package( ecbuild QUIET )
include( ecbuild_system NO_POLICY_SCOPE )
ecbuild_declare_project()
list( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake )
set( CMAKE_DIRECTORY_LABELS ${PROJECT_NAME} )

include( gsw_compiler_flags )

## Sources
add_subdirectory(modules)
add_subdirectory(toolbox)
add_subdirectory(test)

list ( APPEND gsw_src_files
    ${gsw_modules_files}
    ${gsw_toolbox_files}
)

## Library
ecbuild_add_library( TARGET   gsw
                     SOURCES  ${gsw_src_files}
                    )

## Fortran modules
set(MODULE_DIR ${PROJECT_NAME}/module)
set_target_properties(${PROJECT_NAME} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/${MODULE_DIR})
install(DIRECTORY ${CMAKE_BINARY_DIR}/${MODULE_DIR}/ DESTINATION ${CMAKE_INSTALL_LIBDIR}/${MODULE_DIR})
target_include_directories(${PROJECT_NAME} INTERFACE
                            $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/${MODULE_DIR}>
                            $<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}/${MODULE_DIR}>)


## Export package
ecbuild_install_project( NAME ${PROJECT_NAME} )
ecbuild_print_summary()
