 
 
cmake_minimum_required(VERSION 2.6) 

project(spblas) 

FILE(GLOB ALL_HDR "*.h")
FILE(GLOB ALL_SRC "*.cc")
set(ALL_SRC_FILES ${ALL_SRC} ${ALL_HDR})
SOURCE_GROUP("Common" FILES
                ${ALL_HDR}
                ${ALL_SRC}) 

include_directories( ${SPARSELIB_ROOT_DIR}/include) 

add_library(spblas STATIC ${ALL_SRC_FILES} )