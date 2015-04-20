cmake_minimum_required (VERSION 2.8.6)
project (KratosChemist)



find_package( Qt4 REQUIRED)
set( QT_DONT_USE_QTGUI true )
set(QT_USE_QTXML true)
include( ${QT_USE_FILE} )


SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pthread" )
add_definitions(-DIsMultiThreaded)
message(STATUS "configuring for gcc. You should use version 4.7 or newer")


find_package(Eigen3 REQUIRED )
message (hola2)
if(EIGEN3_FOUND)
    message(STATUS "eigen found")  
    include_directories(${EIGEN3_INCLUDE_DIR})
endif()

FILE(GLOB ALL_HDR "*.h")
FILE(GLOB ALL_SRC "*.cpp")
set(ALL_SRC_FILES ${ALL_SRC} ${ALL_HDR})
SOURCE_GROUP("Common" FILES
                ${ALL_HDR}
                ${ALL_SRC}) 

ADD_LIBRARY (ChemistFranch SHARED ${ALL_SRC_FILES})

target_link_libraries(ChemistFranch ${QT_LIBRARIES})

IF (NOT Chemist)
   message ("Not Found Chemist")
ENDIF ()
SET_TARGET_PROPERTIES(ChemistFranch PROPERTIES LINKER_LANGUAGE C)

INSTALL (FILES ${ChemistFranch} DESTINATION lib)


