# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.3

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/enric/kratos

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/enric/kratos/cmake_build

# Include any dependencies generated for this target.
include applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/depend.make

# Include the progress variables for this target.
include applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/progress.make

# Include the compile flags for this target's objects.
include applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/flags.make

applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/manager_application.cpp.o: applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/flags.make
applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/manager_application.cpp.o: ../applications/ManagerApplication/manager_application.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/enric/kratos/cmake_build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/manager_application.cpp.o"
	cd /home/enric/kratos/cmake_build/applications/ManagerApplication && /usr/bin/g++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/KratosManagerApplication.dir/manager_application.cpp.o -c /home/enric/kratos/applications/ManagerApplication/manager_application.cpp

applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/manager_application.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/KratosManagerApplication.dir/manager_application.cpp.i"
	cd /home/enric/kratos/cmake_build/applications/ManagerApplication && /usr/bin/g++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/enric/kratos/applications/ManagerApplication/manager_application.cpp > CMakeFiles/KratosManagerApplication.dir/manager_application.cpp.i

applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/manager_application.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/KratosManagerApplication.dir/manager_application.cpp.s"
	cd /home/enric/kratos/cmake_build/applications/ManagerApplication && /usr/bin/g++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/enric/kratos/applications/ManagerApplication/manager_application.cpp -o CMakeFiles/KratosManagerApplication.dir/manager_application.cpp.s

applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/manager_application.cpp.o.requires:

.PHONY : applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/manager_application.cpp.o.requires

applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/manager_application.cpp.o.provides: applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/manager_application.cpp.o.requires
	$(MAKE) -f applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/build.make applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/manager_application.cpp.o.provides.build
.PHONY : applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/manager_application.cpp.o.provides

applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/manager_application.cpp.o.provides.build: applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/manager_application.cpp.o


applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/custom_python/manager_python_application.cpp.o: applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/flags.make
applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/custom_python/manager_python_application.cpp.o: ../applications/ManagerApplication/custom_python/manager_python_application.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/enric/kratos/cmake_build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/custom_python/manager_python_application.cpp.o"
	cd /home/enric/kratos/cmake_build/applications/ManagerApplication && /usr/bin/g++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/KratosManagerApplication.dir/custom_python/manager_python_application.cpp.o -c /home/enric/kratos/applications/ManagerApplication/custom_python/manager_python_application.cpp

applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/custom_python/manager_python_application.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/KratosManagerApplication.dir/custom_python/manager_python_application.cpp.i"
	cd /home/enric/kratos/cmake_build/applications/ManagerApplication && /usr/bin/g++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/enric/kratos/applications/ManagerApplication/custom_python/manager_python_application.cpp > CMakeFiles/KratosManagerApplication.dir/custom_python/manager_python_application.cpp.i

applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/custom_python/manager_python_application.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/KratosManagerApplication.dir/custom_python/manager_python_application.cpp.s"
	cd /home/enric/kratos/cmake_build/applications/ManagerApplication && /usr/bin/g++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/enric/kratos/applications/ManagerApplication/custom_python/manager_python_application.cpp -o CMakeFiles/KratosManagerApplication.dir/custom_python/manager_python_application.cpp.s

applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/custom_python/manager_python_application.cpp.o.requires:

.PHONY : applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/custom_python/manager_python_application.cpp.o.requires

applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/custom_python/manager_python_application.cpp.o.provides: applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/custom_python/manager_python_application.cpp.o.requires
	$(MAKE) -f applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/build.make applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/custom_python/manager_python_application.cpp.o.provides.build
.PHONY : applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/custom_python/manager_python_application.cpp.o.provides

applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/custom_python/manager_python_application.cpp.o.provides.build: applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/custom_python/manager_python_application.cpp.o


applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/custom_python/add_custom_utilities_to_python.cpp.o: applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/flags.make
applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/custom_python/add_custom_utilities_to_python.cpp.o: ../applications/ManagerApplication/custom_python/add_custom_utilities_to_python.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/enric/kratos/cmake_build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/custom_python/add_custom_utilities_to_python.cpp.o"
	cd /home/enric/kratos/cmake_build/applications/ManagerApplication && /usr/bin/g++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/KratosManagerApplication.dir/custom_python/add_custom_utilities_to_python.cpp.o -c /home/enric/kratos/applications/ManagerApplication/custom_python/add_custom_utilities_to_python.cpp

applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/custom_python/add_custom_utilities_to_python.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/KratosManagerApplication.dir/custom_python/add_custom_utilities_to_python.cpp.i"
	cd /home/enric/kratos/cmake_build/applications/ManagerApplication && /usr/bin/g++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/enric/kratos/applications/ManagerApplication/custom_python/add_custom_utilities_to_python.cpp > CMakeFiles/KratosManagerApplication.dir/custom_python/add_custom_utilities_to_python.cpp.i

applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/custom_python/add_custom_utilities_to_python.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/KratosManagerApplication.dir/custom_python/add_custom_utilities_to_python.cpp.s"
	cd /home/enric/kratos/cmake_build/applications/ManagerApplication && /usr/bin/g++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/enric/kratos/applications/ManagerApplication/custom_python/add_custom_utilities_to_python.cpp -o CMakeFiles/KratosManagerApplication.dir/custom_python/add_custom_utilities_to_python.cpp.s

applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/custom_python/add_custom_utilities_to_python.cpp.o.requires:

.PHONY : applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/custom_python/add_custom_utilities_to_python.cpp.o.requires

applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/custom_python/add_custom_utilities_to_python.cpp.o.provides: applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/custom_python/add_custom_utilities_to_python.cpp.o.requires
	$(MAKE) -f applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/build.make applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/custom_python/add_custom_utilities_to_python.cpp.o.provides.build
.PHONY : applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/custom_python/add_custom_utilities_to_python.cpp.o.provides

applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/custom_python/add_custom_utilities_to_python.cpp.o.provides.build: applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/custom_python/add_custom_utilities_to_python.cpp.o


# Object files for target KratosManagerApplication
KratosManagerApplication_OBJECTS = \
"CMakeFiles/KratosManagerApplication.dir/manager_application.cpp.o" \
"CMakeFiles/KratosManagerApplication.dir/custom_python/manager_python_application.cpp.o" \
"CMakeFiles/KratosManagerApplication.dir/custom_python/add_custom_utilities_to_python.cpp.o"

# External object files for target KratosManagerApplication
KratosManagerApplication_EXTERNAL_OBJECTS =

applications/ManagerApplication/KratosManagerApplication.so: applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/manager_application.cpp.o
applications/ManagerApplication/KratosManagerApplication.so: applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/custom_python/manager_python_application.cpp.o
applications/ManagerApplication/KratosManagerApplication.so: applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/custom_python/add_custom_utilities_to_python.cpp.o
applications/ManagerApplication/KratosManagerApplication.so: applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/build.make
applications/ManagerApplication/KratosManagerApplication.so: kratos/libKratosCore.so
applications/ManagerApplication/KratosManagerApplication.so: /home/enric/boost_1_57_0/stage/lib/libboost_python.so.1.57.0
applications/ManagerApplication/KratosManagerApplication.so: /usr/lib/python3.4/config-3.4m-x86_64-linux-gnu/libpython3.4m.so
applications/ManagerApplication/KratosManagerApplication.so: external_libraries/gidpost/libgidpost.a
applications/ManagerApplication/KratosManagerApplication.so: external_libraries/zlib/libz.so.1.2.8
applications/ManagerApplication/KratosManagerApplication.so: applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/enric/kratos/cmake_build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX shared library KratosManagerApplication.so"
	cd /home/enric/kratos/cmake_build/applications/ManagerApplication && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/KratosManagerApplication.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/build: applications/ManagerApplication/KratosManagerApplication.so

.PHONY : applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/build

applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/requires: applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/manager_application.cpp.o.requires
applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/requires: applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/custom_python/manager_python_application.cpp.o.requires
applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/requires: applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/custom_python/add_custom_utilities_to_python.cpp.o.requires

.PHONY : applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/requires

applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/clean:
	cd /home/enric/kratos/cmake_build/applications/ManagerApplication && $(CMAKE_COMMAND) -P CMakeFiles/KratosManagerApplication.dir/cmake_clean.cmake
.PHONY : applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/clean

applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/depend:
	cd /home/enric/kratos/cmake_build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/enric/kratos /home/enric/kratos/applications/ManagerApplication /home/enric/kratos/cmake_build /home/enric/kratos/cmake_build/applications/ManagerApplication /home/enric/kratos/cmake_build/applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : applications/ManagerApplication/CMakeFiles/KratosManagerApplication.dir/depend

