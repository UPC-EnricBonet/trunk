cmake_minimum_required(VERSION 2.6) 

project(Cheproo++) 


set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake) 
set(CHEPROO_ROOT_DIR ${PROJECT_SOURCE_DIR}) 




#add the tests
if( IS_DIRECTORY ${PROJECT_SOURCE_DIR}/testing)
    add_subdirectory ( ${PROJECT_SOURCE_DIR}/testing )
endif()

#add the executable
if( IS_DIRECTORY ${PROJECT_SOURCE_DIR}/src)
    add_subdirectory ( ${PROJECT_SOURCE_DIR}/src )
endif()


#add the library
if( IS_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/cheprooLib )
	message(STATUS "library found")
    add_subdirectory ( ${PROJECT_SOURCE_DIR}/lib/cheprooLib )
else()
	message(STATUS "library source not found")
endif()
