CMAKE_MINIMUM_REQUIRED(VERSION 2.6)  
  
PROJECT (CheprooPlusPlus)

file(GLOB CHEPROO_SRC "*.cpp") #glob looks for all files in current directory. 
 message(STATUS  ${CHEPROO_ROOT_DIR})
include_directories( ${PROJECT_SOURCE_DIR} )

SET( EXECUTABLE_OUTPUT_PATH ${CHEPROO_ROOT_DIR}/bin )
SET( LIBRARY_OUTPUT_PATH ${CHEPROO_ROOT_DIR}/lib )

if (IsQt5)
	FIND_PACKAGE( Qt5Core REQUIRED)
else()
	find_package( Qt4 REQUIRED)
	set( QT_DONT_USE_QTGUI true )
	set(QT_USE_QTXML true)
	include( ${QT_USE_FILE} )
endif()

 
include_directories(${COMMONLIB_SRC_DIR} ${CHEPROO_ROOT_DIR}/lib/proostfiles ${CHEPROO_ROOT_DIR}/lib/cheprooLib) 
 
			
ADD_EXECUTABLE (cheprooPlusPlus ${CHEPROO_SRC} )

TARGET_LINK_LIBRARIES(cheprooPlusPlus commonLib CheprooLib  ${QT_LIBRARIES} )

 
 # Use the Widgets module from Qt 5.
if (IsQt5)
 qt5_use_modules(cheprooPlusPlus Core)		  
endif()
     