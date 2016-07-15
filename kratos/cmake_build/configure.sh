# this is an example file...please DO NOT MODIFY if you don't want to do it for everyone
#to use it, copy it to another file and run it

# additional compiler flags could be added customizing the corresponding var, for example
# -DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -msse3 ". Note that the "defaults are already correctly coded"
#so we should ass here only machine specific stuff

#an effort is made to autodetect all of the libraries needed HOWEVER:
#METIS_APPLICATION needs the var PARMETIS_ROOT_DIR to be specified by the user (not needed if the app is set to OFF)
#TRILINOS_APPLICATION needs the var PARMETIS_ROOT_DIR to be specified by the user (not needed if the app is set to OFF)
#MKL_SOLVERS_APPLICATION needs the var MKLSOLVER_INCLUDE_DIR and MKLSOLVER_LIB_DIR to be specified by the user (not needed if the app is set to OFF)
#note that the MKLSOLVER_LIB_DIR should include /lib/em64t. This is needed as intel is changing location of mkl at every update of the compiler!!

#the user should also note that the symbol "\" marks that the command continues on the next line. IT SHOULD ONLY BE FOLLOWED
#BY the "ENTER" and NOT by any space!!

clear

#you may want to decomment this the first time you compile
rm CMakeCache.txt
rm *.cmake
rm -rf CMakeFiles\

cmake ..  									                            \
-DCMAKE_C_COMPILER=/usr/bin/gcc							                \
-DCMAKE_INSTALL_RPATH="/home/enric/kratos/libs"                      \
-DCMAKE_MODULE_PATH="/usr/share/cmake-2.8/Modules/FindPythonLibs.cmake" \
-DCMAKE_INSTALL_RPATH_USE_LINK_PATH=TRUE                                \
-DCMAKE_CXX_COMPILER=/usr/bin/g++ 						                \
-DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -msse3 -std=c++11"                   \
-DBOOST_ROOT="/home/enric/boost_1_57_0"                                    \
-DPYTHON_LIBRARY="/usr/lib/python3.4/config-3.4m-x86_64-linux-gnu/libpython3.4m.so"        \
-DPYTHON_INCLUDE_DIR="/usr/include/python3.4"                                       \
-DCMAKE_C_FLAGS="${CMAKE_C_FLAGS} -msse3 " 					            \
-DCMAKE_BUILD_TYPE=Debug                                                                    \
-DINCOMPRESSIBLE_FLUID_APPLICATION=OFF  						            \
-DMESHING_APPLICATION=OFF 							                    \
-DEXTERNAL_SOLVERS_APPLICATION=OFF						                \
-DPFEM_APPLICATION=OFF 								                    \
-DSTRUCTURAL_APPLICATION=OFF 							                \
-DCONVECTION_DIFFUSION_APPLICATION=OFF						            \
-DFLUID_DYNAMICS_APPLICATION=OFF 						                \
-DALE_APPLICATION=OFF 								                    \
-DFSI_APPLICATION=OFF 								                    \
-DOPENCL_APPLICATION=OFF							                    \
-DMIXED_ELEMENT_APPLICATION=OFF							                \
-DSTRUCTURAL_MECHANICS_APPLICATION=OFF                                                            \
-DMKL_SOLVERS_APPLICATION=OFF							                \
-DMKLSOLVER_INCLUDE_DIR="/opt/intel/Compiler/11.1/072/mkl/include"	    \
-DMKLSOLVER_LIB_DIR="/opt/intel/Compiler/11.1/072/mkl/lib/em64t"	    \
-DMETIS_APPLICATION=OFF								                    \
-DPARMETIS_ROOT_DIR="/home/enric/compiled_libraries/ParMetis-3.1.1"  \
-DTRILINOS_APPLICATION=OFF							                    \
-DTRILINOS_ROOT="/home/enric/compiled_libraries/trilinos-10.2.0"		\
-DDEM_APPLICATION=OFF								                    \
-DSWIMMING_DEM_APPLICATION=OFF                                          \
-DSOLID_MECHANICS_APPLICATION=OFF   					                    \
-DPFEM_SOLID_MECHANICS_APPLICATION=OFF 					                \
-DINSTALL_EMBEDDED_PYTHON=OFF                                           \
-DPFEM_2_APPLICATION=OFF                                           \
-DSOLUTETRANSPORT_APPLICATION=OFF                                          \
-DPUREDIFFUSION_APPLICATION=OFF                                           \
-DPUREDIFFUSIONTEMPLATES_APPLICATION=OFF                                           \
-DFLOWTRANSPORT_APPLICATION=OFF                                           \
-DFLOW_TRANSPORT_APPLICATION=OFF                                           \
-DPOROUS_MEDIA_APPLICATION=OFF                                           \
-DCO2_CELIA_APPLICATION=OFF                                           \
-DMULTISCALE_APPLICATION=OFF                                         \
-DFLUIDS_FINITE_ELEM_APPLICATION=OFF                                         \
-DUNIT_BENCHMARK_TEST_APPLICATION=OFF                                     \
-DMANAGER_APPLICATION=ON                                     \
-DHENRY_APPLICATION=ON                                     \
#decomment this to have it verbose
# make VERBOSE=1 -j4
make -j2
make install
