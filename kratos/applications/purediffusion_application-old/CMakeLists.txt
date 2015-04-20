set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosPureDiffusionApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

## generate variables with the sources
set( KRATOS_PUREDIFFUSION_APPLICATION_SOURCES     
	${CMAKE_CURRENT_SOURCE_DIR}//purediffusion_application.cpp 

	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/purediffusion_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/poisson_2d.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/pointsource.cpp
	#${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/your_new_element.cpp
	#${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/your_new_element_2d.h
)

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosPureDiffusionApplication SHARED ${KRATOS_PUREDIFFUSION_APPLICATION_SOURCES})
target_link_libraries(KratosPureDiffusionApplication KratosCore )
set_target_properties(KratosPureDiffusionApplication PROPERTIES PREFIX "")
install(TARGETS KratosPureDiffusionApplication DESTINATION libs )


# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# changing the .dll suffix to .pyd 
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosPureDiffusionApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/PureDiffusionApplication.py" DESTINATION KratosMultiphysics )
