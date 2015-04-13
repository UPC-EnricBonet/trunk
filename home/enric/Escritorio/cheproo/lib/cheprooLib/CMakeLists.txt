cmake_minimum_required(VERSION 2.6)

project(CheprooLib)

#get eigen
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake) 


find_package( PkgConfig )
pkg_check_modules( EIGEN3 REQUIRED eigen3 )
include_directories( ${EIGEN3_INCLUDE_DIRS} )


					 
FILE(GLOB ALL_FILES_SRC *.h *.cpp)

SOURCE_GROUP("" FILES 
				${ALL_FILES_SRC})

 
 
set(CMAKE_MODULE_PATH ${CHEPROO_ROOT_DIR}/cmake)
if (IsQt5)
	find_package(Qt5Core REQUIRED)
	FIND_PACKAGE( Qt5Xml REQUIRED)
else()
	find_package( Qt4 REQUIRED)
	set( QT_DONT_USE_QTGUI true )
	set(QT_USE_QTXML true)
	include( ${QT_USE_FILE} )
endif()



add_definitions(-D_CRT_SECURE_NO_WARNINGS)

include_directories(${COMMONLIB_SRC_DIR}
					${CHEPROO_ROOT_DIR}/lib/cheprooLib)  


add_library(CheprooLib STATIC 
			    ${ALL_FILES_SRC} )
				
target_link_libraries(CheprooLib commonLib ${QT_LIBRARIES}  )

if (IsQt5)
	 # Use the Widgets module from Qt 5.
	 qt5_use_modules(CheprooLib Core Xml)		 
endif()
