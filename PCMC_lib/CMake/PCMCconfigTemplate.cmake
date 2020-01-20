FIND_PACKAGE(GSL REQUIRED)
FIND_PACKAGE(MPI REQUIRED)

LINK_DIRECTORIES(${CMAKE_CURRENT_LIST_DIR}/..)
LIST(APPEND PCMC_LIBRARIES pcmc)

INCLUDE_DIRECTORIES(${GSL_INCLUDE_DIRS})
LIST(APPEND PCMC_LIBRARIES ${GSL_LIBRARIES})

INCLUDE_DIRECTORIES(${MPI_CXX_INCLUDE_PATH})
LIST(APPEND PCMC_LIBRARIES ${MPI_CXX_LIBRARIES})

LIST(APPEND CMAKE_CXX_FLAGS ${MPI_CXX_COMPILE_FLAGS})