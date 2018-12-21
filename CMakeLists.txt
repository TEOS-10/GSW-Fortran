# (C) Copyright 2017 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

################################################################################
# Gibbs-SeaWater (GSW) Oceanographic Toolbox in Fortran
# https://github.com/aerorahul/GSW-Fortran -b jedi-marine
# branch jedi-marine forked from
# https://github.com/TEOS-10/GSW-Fortran
################################################################################

cmake_minimum_required( VERSION 3.3.2 FATAL_ERROR )

project( gsw Fortran )

set( CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH} )

set( ECBUILD_DEFAULT_BUILD_TYPE Release )
set( ENABLE_OS_TESTS           OFF CACHE BOOL "Disable OS tests" FORCE )
set( ENABLE_LARGE_FILE_SUPPORT OFF CACHE BOOL "Disable testing of large file support" FORCE )

include( ecbuild_system NO_POLICY_SCOPE )

ecbuild_requires_macro_version( 2.7 )

################################################################################
# Project
################################################################################

ecbuild_declare_project()

ecbuild_enable_fortran( REQUIRED )

################################################################################
# Dependencies
################################################################################


################################################################################
# Definitions
################################################################################


################################################################################
# Export package info
################################################################################

list( APPEND GSW_TPLS )

set( GSW_INCLUDE_DIRS ${CMAKE_Fortran_MODULE_DIRECTORY})
set( GSW_LIBRARIES gsw )

get_directory_property( GSW_DEFINITIONS COMPILE_DEFINITIONS )

foreach( _tpl ${GSW_TPLS} )
  string( TOUPPER ${_tpl} TPL )
  list( APPEND GSW_EXTRA_DEFINITIONS   ${${TPL}_DEFINITIONS}  ${${TPL}_TPL_DEFINITIONS}  )
  list( APPEND GSW_EXTRA_INCLUDE_DIRS  ${${TPL}_INCLUDE_DIRS} ${${TPL}_TPL_INCLUDE_DIRS} )
  list( APPEND GSW_EXTRA_LIBRARIES     ${${TPL}_LIBRARIES}    ${${TPL}_TPL_LIBRARIES}    )
endforeach()

################################################################################
# Sources
################################################################################

set( GSW_LINKER_LANGUAGE Fortran )

include( gsw_compiler_flags )
include_directories( ${GSW_INCLUDE_DIRS} )

add_subdirectory(modules)
add_subdirectory(toolbox)
add_subdirectory(test)

list ( APPEND gsw_src_files
	${gsw_modules_files}
	${gsw_toolbox_files}
)

################################################################################
# Library
################################################################################

ecbuild_add_library( TARGET   gsw
                     SOURCES  ${gsw_src_files}
                     LINKER_LANGUAGE ${GSW_LINKER_LANGUAGE}
                    )

################################################################################
# Finalise configuration
################################################################################

# prepares a tar.gz of the sources and/or binaries
ecbuild_install_project( NAME gsw )

# print the summary of the configuration
ecbuild_print_summary()
