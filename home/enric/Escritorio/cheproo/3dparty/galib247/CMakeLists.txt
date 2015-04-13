cmake_minimum_required(VERSION 2.6) 

project(GAlib)

#find source files
file(GLOB GALIB_HDR ga/*.h)
file(GLOB GALIB_SRC ga/*.C)

add_definitions(-D_CRT_SECURE_NO_WARNINGS)

include_directories(${PROJECT_SOURCE_DIR})

list(REMOVE_ITEM GALIB_SRC ${PROJECT_SOURCE_DIR}/ga/GATree.C 
                           ${PROJECT_SOURCE_DIR}/ga/GAList.C
                           ${PROJECT_SOURCE_DIR}/ga/GAListGenome.C
                           ${PROJECT_SOURCE_DIR}/ga/GAAllele.C
                           ${PROJECT_SOURCE_DIR}/ga/GA1DArrayGenome.C
                           ${PROJECT_SOURCE_DIR}/ga/GA2DArrayGenome.C
                           ${PROJECT_SOURCE_DIR}/ga/GA3DArrayGenome.C)
		
add_library(GALib STATIC ${GALIB_HDR} ${GALIB_SRC})