set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosUnit_Benchmark_TestApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

## generate variables with the sources
set( KRATOS_UNIT_BENCHMARK_TEST_APPLICATION_SOURCES     
	${CMAKE_CURRENT_SOURCE_DIR}//unit_benchmark_test_application.cpp 

	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/unit_benchmark_test_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
	#${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/your_new_element.cpp
	#${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/your_new_element_2d.h
)

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosUnit_Benchmark_TestApplication SHARED ${KRATOS_UNIT_BENCHMARK_TEST_APPLICATION_SOURCES})
target_link_libraries(KratosUnit_Benchmark_TestApplication KratosCore )
set_target_properties(KratosUnit_Benchmark_TestApplication PROPERTIES PREFIX "")
install(TARGETS KratosUnit_Benchmark_TestApplication DESTINATION libs )


# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# changing the .dll suffix to .pyd 
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosUnit_Benchmark_TestApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/Unit_Benchmark_TestApplication.py" DESTINATION KratosMultiphysics )
