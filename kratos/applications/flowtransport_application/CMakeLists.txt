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

include_directories( ${CMAKE_SOURCE_DIR}/external/externalForLibrary/src)
#Fortran introduction
enable_language (Fortran)


# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # gfortran
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  # ifort (untested)
  set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g")
elseif (Fortran_COMPILER_NAME MATCHES "g77")
  # g77
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -m32")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -m32")
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")


FIND_LIBRARY (CHEMIST_LIB  Chemist ${CMAKE_CURRENT_SOURCE_DIR}/external/Chemist )
FIND_LIBRARY (CHEMIST_LIB_FRANCH  ChemistFranch ${CMAKE_CURRENT_SOURCE_DIR}/external/Chemist )
FIND_LIBRARY (CHEMIST_FOR_LIB ForChemist ${CMAKE_CURRENT_SOURCE_DIR}/external/externalForLibrary)

add_library (A OBJECT ${CMAKE_CURRENT_SOURCE_DIR}/external/externalForLibrary/prove.f)


IF (CHEMIST_LIB_FRANCH)
    message (STATUS "found ChemistFranch")
ENDIF()

IF (CHEMIST_LIB)
    message (STATUS "found Chemist")
ENDIF()

IF (CHEMIST_FOR_LIB)
    message (STATUS "found ForChemist")
ENDIF()

## generate variables with the sources
set( KRATOS_FLOWTRANSPORT_APPLICATION_SOURCES     
	${CMAKE_CURRENT_SOURCE_DIR}/flowtransport_application.cpp 

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
        #${CMAKE_CURRENT_SOURCE_DIR}/external/externalForLibrary/prove.f
        
)

message (${CMAKE_Fortran_COMPILER})
###############################################################
## define library Kratos which defines the basic python interface
message("------------------------------> found chemist lib =" ${CHEMIST_LIB} )
message("------------------------------> found chemist lib =" ${CHEMIST_LIB_FRANCH} )
message("------------------------------> found chemist lib for =" ${CHEMIST_FOR_LIB} )
add_library(KratosFlowTransportApplication SHARED ${KRATOS_FLOWTRANSPORT_APPLICATION_SOURCES} $<TARGET_OBJECTS:A> )
target_link_libraries(KratosFlowTransportApplication KratosCore ${CHEMIST_LIB} ${CHEMIST_LIB_FRANCH} ${QT_LIBRARIES} )
set_target_properties(KratosFlowTransportApplication PROPERTIES PREFIX "")
install(TARGETS KratosFlowTransportApplication DESTINATION libs )
install(FILES ${CHEMIST_LIB} DESTINATION libs )
install(FILES ${CHEMIST_LIB_FRANCH} DESTINATION libs )
install(FILES ${CHEMIST_FOR_LIB} DESTINATION libs )

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
