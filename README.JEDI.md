GSW-Fortran Toolbox Interface for Marine UFO within JEDI framework

(C) Copyright 2018 UCAR.

This software is licensed under the terms of the Apache Licence Version 2.0
which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

--- Requirements ---

See OOPS requirements

--- Building ---

The variables ${SRC_GSW} and ${BUILD} below must be defined for your
environment.
Note: It is good practice to build the code outside of the source tree.

The lines below can be copied into a script or executed manually:

Define environment

    export SRC_GSW=/path/to/source/gsw
    export BUILD=/path/to/build

    export PATH=${PATH}:/path/to/ecbuild/bin

Build GSW-Fortran library

    rm -rf ${BUILD}/gsw; mkdir ${BUILD}/gsw; cd ${BUILD}/gsw
    ecbuild --build=release ${SRC_GSW}
    make -j4
