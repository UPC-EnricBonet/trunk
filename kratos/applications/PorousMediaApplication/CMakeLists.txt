set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosPorousMediaApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
 
## generate variables with the sources
set( KRATOS_POROUS_MEDIA_APPLICATION_SOURCES     
	${CMAKE_CURRENT_SOURCE_DIR}//porous_media_application.cpp 

	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/porous_media_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp

        ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/flow_2d.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/flowtrans_2d.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/flowpressuretrans_2d.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/sinksource.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/sinksourcepressure.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/leakage.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/massflow.cpp 
) 

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosPorousMediaApplication SHARED ${KRATOS_POROUS_MEDIA_APPLICATION_SOURCES})
target_link_libraries(KratosPorousMediaApplication KratosCore )
set_target_properties(KratosPorousMediaApplication PROPERTIES PREFIX "")
install(TARGETS KratosPorousMediaApplication DESTINATION libs )


# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# changing the .dll suffix to .pyd 
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosPorousMediaApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/PorousMediaApplication.py" DESTINATION KratosMultiphysics )
 
