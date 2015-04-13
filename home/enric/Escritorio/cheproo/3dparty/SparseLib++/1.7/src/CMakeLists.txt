 
cmake_minimum_required(VERSION 2.6) 

project(sparselib) 

FILE(GLOB ALL_HDR "*.h")
FILE(GLOB ALL_SRC "*.cc")
set(ALL_SRC_FILES ${ALL_SRC} ${ALL_HDR})
SOURCE_GROUP("Common" FILES
                ${ALL_HDR}
                ${ALL_SRC}) 
include_directories( ${SPARSELIB_ROOT_DIR}/include
		     ${SPARSELIB_ROOT_DIR}/mv/src
                     ${SPARSELIB_ROOT_DIR}/spblas)


add_library(sparselib STATIC ${ALL_SRC_FILES} )
target_link_libraries( sparselib mv spblas)