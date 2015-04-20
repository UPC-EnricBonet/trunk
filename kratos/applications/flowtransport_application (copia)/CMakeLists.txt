set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosFlowTransportApplication ****")


find_package( Qt4 REQUIRED)
set( QT_DONT_USE_QTGUI true )
set(QT_USE_QTXML true)
set(QT_USE_QTCORE true)
include( ${QT_USE_FILE} )




SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pthread" )
add_definitions(-DIsMultiThreaded)
message(STATUS "configuring for gcc. You should use version 4.7 or newer")

message(STATUS "Viva Kratos!!!!!")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

#include_directories( /usr/include/qt4 )
#include_directories( /usr/include/qt4/QtCore )
#include_directories( /usr/include/qt4/QtXml )

include_directories( /usr/include/eigen3 )

FIND_LIBRARY (CHEMIST_LIB  Chemist ${CMAKE_CURRENT_SOURCE_DIR}/external/Chemist )
FIND_LIBRARY (CHEMIST_LIB_FRANCH  ChemistFranch ${CMAKE_CURRENT_SOURCE_DIR}/external/Chemist )

IF (CHEMIST_LIB_FRANCH)
    message (STATUS "found ChemistFranch")
ENDIF()

IF (CHEMIST_LIB)
    message (STATUS "found Chemist")
ENDIF()

## generate variables with the sources
set( KRATOS_FLOWTRANSPORT_APPLICATION_SOURCES     
	${CMAKE_CURRENT_SOURCE_DIR}//flowtransport_application.cpp 

	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/flowtransport_python_application.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/poisson_2d.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/diffusion_2d.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/darcy.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/pointflow.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/pointsource.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/pointheat.cpp
  
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/flowDarcy.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/transport.cpp

        ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/sinksource.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/sinksourcepressure.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/leakage.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/massflow.cpp
        
)

###############################################################
## define library Kratos which defines the basic python interface
message("------------------------------> found chemist lib =" ${CHEMIST_LIB} )
message("------------------------------> found chemist lib =" ${CHEMIST_LIB_FRANCH} )
add_library(KratosFlowTransportApplication SHARED ${KRATOS_FLOWTRANSPORT_APPLICATION_SOURCES})
target_link_libraries(KratosFlowTransportApplication KratosCore ${CHEMIST_LIB} ${CHEMIST_LIB_FRANCH} ${QT_LIBRARIES} )
set_target_properties(KratosFlowTransportApplication PROPERTIES PREFIX "")
install(TARGETS KratosFlowTransportApplication DESTINATION libs )
install(FILES ${CHEMIST_LIB} DESTINATION libs )
install(FILES ${CHEMIST_LIB_FRANCH} DESTINATION libs )

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# changing the .dll suffix to .pyd 
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosFlowTransportApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/FlowTransportApplication.py" DESTINATION KratosMultiphysics )
