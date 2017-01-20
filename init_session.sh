#!/bin/bash

INC_DIR=$PWD/PCMC_lib/include

#export PROJECT_DIR

if [ -z "$CPLUS_INCLUDE_PATH" ]
then
    CPLUS_INCLUDE_PATH=$INC_DIR*/**/;
else
    CPLUS_INCLUDE_PATH+=:$INC_DIR*/**/;
fi
echo $CPLUS_INCLUDE_PATH
#export CPLUS_INCLUDE_PATH