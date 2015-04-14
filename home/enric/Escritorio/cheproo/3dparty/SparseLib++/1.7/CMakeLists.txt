 
cmake_minimum_required(VERSION 2.6) 

project(allsparseLib)

set(SPARSELIB_ROOT_DIR ${GHSSOFT_PATH}/3dparty/SparseLib++/1.7/)
set( LIBRARY_OUTPUT_PATH ${GHSSOFT_PATH}/lib )
  
add_definitions(-DCOMPLEX=std::complex<double>)
 
if( IS_DIRECTORY ${SPARSELIB_ROOT_DIR}/spblas )
    add_subdirectory ( spblas ) #add_subdirectory runs another cmake file
endif()
if( IS_DIRECTORY ${SPARSELIB_ROOT_DIR}/mv/src )
    add_subdirectory ( mv/src )
endif()
if( IS_DIRECTORY ${SPARSELIB_ROOT_DIR}/src )
    add_subdirectory ( src)
endif()
