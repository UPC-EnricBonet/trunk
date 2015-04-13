cmake_minimum_required(VERSION 2.6)

project(CheprooUnitTest)

set(CHEPROO_ROOT_DIR ${PROJECT_SOURCE_DIR}/..)


file(GLOB ALL_SYSTEST_SRC "${CHEPROO_ROOT_DIR}/testing/UnitTest/*" )

set( ALL_SYSTEST ${ALL_SYSTEST_HDR} ${ALL_SYSTEST_SRC})
source_group("General" FILES
			 ${GENERAL_HDR}
			 ${GENERAL_SRC}
			 ${ALL_SYSTEST}) 



include_directories(${COMMONLIB_SRC_DIR}
					 ${CHEPROO_ROOT_DIR}/lib/cheprooLib
                     ${CHEPROO_ROOT_DIR}/testing/UnitTest
					 ${CHEPROO_ROOT_DIR}/src) 
 
set(CMAKE_MODULE_PATH ${CHEPROO_ROOT_DIR}/cmake)

if (IsQt5)
	find_package( Qt5Core REQUIRED)
else()
	find_package( Qt4 REQUIRED)
	set( QT_DONT_USE_QTGUI true )
	set(QT_USE_QTXML true)
	include( ${QT_USE_FILE} )
endif()


set(EIGEN3_INCLUDE_DIR ${CHEPROO_ROOT_DIR}/../3dparty/eigen-3.1.2)
find_package(Eigen3 REQUIRED)
if(EIGEN3_FOUND)
    include_directories(${EIGEN3_INCLUDE_DIR})
endif()


#find_path(CLAPACK_INCLUDE clapack.h ${clapack_DIR} ${clapack_DIR}/include)
#include_directories(${CLAPACK_INCLUDE})

add_definitions(-D_CRT_SECURE_NO_WARNINGS)

add_executable (CheprooUnitTest ${ALL_SYSTEST}
						      ${ALL_SRC} 
						      ${ALL_HDR})
							  
TARGET_LINK_LIBRARIES(CheprooUnitTest CheprooLib commonLib)
 
 # Use the Widgets module from Qt 5.
 if (IsQt5)
	qt5_use_modules(CheprooUnitTest Core)
endif()
