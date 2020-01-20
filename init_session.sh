#!/bin/bash

# Run: $ source init_session.se

PROJECT_DIR=$PWD
INC_DIR=$PWD/PCMC_lib/include
BUILD_DIR=$PWD/PCMC_lib/build

#################
# Include paths #
#################

if [ -z "$CPLUS_INCLUDE_PATH" ]
then
    CPLUS_INCLUDE_PATH=${INC_DIR};
else
    CPLUS_INCLUDE_PATH+=:${INC_DIR};
fi
CPLUS_INCLUDE_PATH+=:${INC_DIR}/Random;
CPLUS_INCLUDE_PATH+=:${INC_DIR}/Energy;
CPLUS_INCLUDE_PATH+=:${INC_DIR}/Quantum;
#echo $CPLUS_INCLUDE_PATH
export CPLUS_INCLUDE_PATH

#########################
# Build time link paths #
#########################
if [ -z "$LIBRARY_PATH" ]
then
    LIBRARY_PATH=${BUILD_DIR};
else
    LIBRARY_PATH+=:${BUILD_DIR};
fi
export LIBRARY_PATH

printf "Initialization session completed.\n"


