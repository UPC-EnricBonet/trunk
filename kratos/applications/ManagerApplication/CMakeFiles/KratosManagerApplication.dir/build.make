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
CMAKE_BINARY_DIR = /home/enric/kratos

# Include any dependencies generated for this target.
include applications/manager_application/CMakeFiles/KratosManagerApplication.dir/depend.make

# Include the progress variables for this target.
include applications/manager_application/CMakeFiles/KratosManagerApplication.dir/progress.make

# Include the compile flags for this target's objects.
include applications/manager_application/CMakeFiles/KratosManagerApplication.dir/flags.make

applications/manager_application/CMakeFiles/KratosManagerApplication.dir/manager_application.cpp.o: applications/manager_application/CMakeFiles/KratosManagerApplication.dir/flags.make
applications/manager_application/CMakeFiles/KratosManagerApplication.dir/manager_application.cpp.o: applications/manager_application/manager_application.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/enric/kratos/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object applications/manager_application/CMakeFiles/KratosManagerApplication.dir/manager_application.cpp.o"
	cd /home/enric/kratos/applications/manager_application && /usr/bin/g++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/KratosManagerApplication.dir/manager_application.cpp.o -c /home/enric/kratos/applications/manager_application/manager_application.cpp

applications/manager_application/CMakeFiles/KratosManagerApplication.dir/manager_application.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/KratosManagerApplication.dir/manager_application.cpp.i"
	cd /home/enric/kratos/applications/manager_application && /usr/bin/g++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/enric/kratos/applications/manager_application/manager_application.cpp > CMakeFiles/KratosManagerApplication.dir/manager_application.cpp.i

applications/manager_application/CMakeFiles/KratosManagerApplication.dir/manager_application.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/KratosManagerApplication.dir/manager_application.cpp.s"
	cd /home/enric/kratos/applications/manager_application && /usr/bin/g++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/enric/kratos/applications/manager_application/manager_application.cpp -o CMakeFiles/KratosManagerApplication.dir/manager_application.cpp.s

applications/manager_application/CMakeFiles/KratosManagerApplication.dir/manager_application.cpp.o.requires:

.PHONY : applications/manager_application/CMakeFiles/KratosManagerApplication.dir/manager_application.cpp.o.requires

applications/manager_application/CMakeFiles/KratosManagerApplication.dir/manager_application.cpp.o.provides: applications/manager_application/CMakeFiles/KratosManagerApplication.dir/manager_application.cpp.o.requires
	$(MAKE) -f applications/manager_application/CMakeFiles/KratosManagerApplication.dir/build.make applications/manager_application/CMakeFiles/KratosManagerApplication.dir/manager_application.cpp.o.provides.build
.PHONY : applications/manager_application/CMakeFiles/KratosManagerApplication.dir/manager_application.cpp.o.provides

applications/manager_application/CMakeFiles/KratosManagerApplication.dir/manager_application.cpp.o.provides.build: applications/manager_application/CMakeFiles/KratosManagerApplication.dir/manager_application.cpp.o


applications/manager_application/CMakeFiles/KratosManagerApplication.dir/custom_python/manager_python_application.cpp.o: applications/manager_application/CMakeFiles/KratosManagerApplication.dir/flags.make
applications/manager_application/CMakeFiles/KratosManagerApplication.dir/custom_python/manager_python_application.cpp.o: applications/manager_application/custom_python/manager_python_application.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/enric/kratos/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object applications/manager_application/CMakeFiles/KratosManagerApplication.dir/custom_python/manager_python_application.cpp.o"
	cd /home/enric/kratos/applications/manager_application && /usr/bin/g++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/KratosManagerApplication.dir/custom_python/manager_python_application.cpp.o -c /home/enric/kratos/applications/manager_application/custom_python/manager_python_application.cpp

applications/manager_application/CMakeFiles/KratosManagerApplication.dir/custom_python/manager_python_application.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/KratosManagerApplication.dir/custom_python/manager_python_application.cpp.i"
	cd /home/enric/kratos/applications/manager_application && /usr/bin/g++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/enric/kratos/applications/manager_application/custom_python/manager_python_application.cpp > CMakeFiles/KratosManagerApplication.dir/custom_python/manager_python_application.cpp.i

applications/manager_application/CMakeFiles/KratosManagerApplication.dir/custom_python/manager_python_application.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/KratosManagerApplication.dir/custom_python/manager_python_application.cpp.s"
	cd /home/enric/kratos/applications/manager_application && /usr/bin/g++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/enric/kratos/applications/manager_application/custom_python/manager_python_application.cpp -o CMakeFiles/KratosManagerApplication.dir/custom_python/manager_python_application.cpp.s

applications/manager_application/CMakeFiles/KratosManagerApplication.dir/custom_python/manager_python_application.cpp.o.requires:

.PHONY : applications/manager_application/CMakeFiles/KratosManagerApplication.dir/custom_python/manager_python_application.cpp.o.requires

applications/manager_application/CMakeFiles/KratosManagerApplication.dir/custom_python/manager_python_application.cpp.o.provides: applications/manager_application/CMakeFiles/KratosManagerApplication.dir/custom_python/manager_python_application.cpp.o.requires
	$(MAKE) -f applications/manager_application/CMakeFiles/KratosManagerApplication.dir/build.make applications/manager_application/CMakeFiles/KratosManagerApplication.dir/custom_python/manager_python_application.cpp.o.provides.build
.PHONY : applications/manager_application/CMakeFiles/KratosManagerApplication.dir/custom_python/manager_python_application.cpp.o.provides

applications/manager_application/CMakeFiles/KratosManagerApplication.dir/custom_python/manager_python_application.cpp.o.provides.build: applications/manager_application/CMakeFiles/KratosManagerApplication.dir/custom_python/manager_python_application.cpp.o


applications/manager_application/CMakeFiles/KratosManagerApplication.dir/custom_python/add_custom_utilities_to_python.cpp.o: applications/manager_application/CMakeFiles/KratosManagerApplication.dir/flags.make
applications/manager_application/CMakeFiles/KratosManagerApplication.dir/custom_python/add_custom_utilities_to_python.cpp.o: applications/manager_application/custom_python/add_custom_utilities_to_python.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/enric/kratos/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object applications/manager_application/CMakeFiles/KratosManagerApplication.dir/custom_python/add_custom_utilities_to_python.cpp.o"
	cd /home/enric/kratos/applications/manager_application && /usr/bin/g++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/KratosManagerApplication.dir/custom_python/add_custom_utilities_to_python.cpp.o -c /home/enric/kratos/applications/manager_application/custom_python/add_custom_utilities_to_python.cpp

applications/manager_application/CMakeFiles/KratosManagerApplication.dir/custom_python/add_custom_utilities_to_python.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/KratosManagerApplication.dir/custom_python/add_custom_utilities_to_python.cpp.i"
	cd /home/enric/kratos/applications/manager_application && /usr/bin/g++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/enric/kratos/applications/manager_application/custom_python/add_custom_utilities_to_python.cpp > CMakeFiles/KratosManagerApplication.dir/custom_python/add_custom_utilities_to_python.cpp.i

applications/manager_application/CMakeFiles/KratosManagerApplication.dir/custom_python/add_custom_utilities_to_python.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/KratosManagerApplication.dir/custom_python/add_custom_utilities_to_python.cpp.s"
	cd /home/enric/kratos/applications/manager_application && /usr/bin/g++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/enric/kratos/applications/manager_application/custom_python/add_custom_utilities_to_python.cpp -o CMakeFiles/KratosManagerApplication.dir/custom_python/add_custom_utilities_to_python.cpp.s

applications/manager_application/CMakeFiles/KratosManagerApplication.dir/custom_python/add_custom_utilities_to_python.cpp.o.requires:

.PHONY : applications/manager_application/CMakeFiles/KratosManagerApplication.dir/custom_python/add_custom_utilities_to_python.cpp.o.requires

applications/manager_application/CMakeFiles/KratosManagerApplication.dir/custom_python/add_custom_utilities_to_python.cpp.o.provides: applications/manager_application/CMakeFiles/KratosManagerApplication.dir/custom_python/add_custom_utilities_to_python.cpp.o.requires
	$(MAKE) -f applications/manager_application/CMakeFiles/KratosManagerApplication.dir/build.make applications/manager_application/CMakeFiles/KratosManagerApplication.dir/custom_python/add_custom_utilities_to_python.cpp.o.provides.build
.PHONY : applications/manager_application/CMakeFiles/KratosManagerApplication.dir/custom_python/add_custom_utilities_to_python.cpp.o.provides

applications/manager_application/CMakeFiles/KratosManagerApplication.dir/custom_python/add_custom_utilities_to_python.cpp.o.provides.build: applications/manager_application/CMakeFiles/KratosManagerApplication.dir/custom_python/add_custom_utilities_to_python.cpp.o


# Object files for target KratosManagerApplication
KratosManagerApplication_OBJECTS = \
"CMakeFiles/KratosManagerApplication.dir/manager_application.cpp.o" \
"CMakeFiles/KratosManagerApplication.dir/custom_python/manager_python_application.cpp.o" \
"CMakeFiles/KratosManagerApplication.dir/custom_python/add_custom_utilities_to_python.cpp.o"

# External object files for target KratosManagerApplication
KratosManagerApplication_EXTERNAL_OBJECTS =

applications/manager_application/KratosManagerApplication.so: applications/manager_application/CMakeFiles/KratosManagerApplication.dir/manager_application.cpp.o
applications/manager_application/KratosManagerApplication.so: applications/manager_application/CMakeFiles/KratosManagerApplication.dir/custom_python/manager_python_application.cpp.o
applications/manager_application/KratosManagerApplication.so: applications/manager_application/CMakeFiles/KratosManagerApplication.dir/custom_python/add_custom_utilities_to_python.cpp.o
applications/manager_application/KratosManagerApplication.so: applications/manager_application/CMakeFiles/KratosManagerApplication.dir/build.make
applications/manager_application/KratosManagerApplication.so: kratos/libKratosCore.so
applications/manager_application/KratosManagerApplication.so: /usr/lib/x86_64-linux-gnu/libboost_python-py27.so.1.54.0
applications/manager_application/KratosManagerApplication.so: /usr/lib/python3.4/config-3.4m-x86_64-linux-gnu/libpython3.4m.so
applications/manager_application/KratosManagerApplication.so: external_libraries/gidpost/libgidpost.a
applications/manager_application/KratosManagerApplication.so: external_libraries/zlib/libz.so.1.2.8
applications/manager_application/KratosManagerApplication.so: applications/manager_application/CMakeFiles/KratosManagerApplication.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/enric/kratos/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX shared library KratosManagerApplication.so"
	cd /home/enric/kratos/applications/manager_application && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/KratosManagerApplication.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
applications/manager_application/CMakeFiles/KratosManagerApplication.dir/build: applications/manager_application/KratosManagerApplication.so

.PHONY : applications/manager_application/CMakeFiles/KratosManagerApplication.dir/build

applications/manager_application/CMakeFiles/KratosManagerApplication.dir/requires: applications/manager_application/CMakeFiles/KratosManagerApplication.dir/manager_application.cpp.o.requires
applications/manager_application/CMakeFiles/KratosManagerApplication.dir/requires: applications/manager_application/CMakeFiles/KratosManagerApplication.dir/custom_python/manager_python_application.cpp.o.requires
applications/manager_application/CMakeFiles/KratosManagerApplication.dir/requires: applications/manager_application/CMakeFiles/KratosManagerApplication.dir/custom_python/add_custom_utilities_to_python.cpp.o.requires

.PHONY : applications/manager_application/CMakeFiles/KratosManagerApplication.dir/requires

applications/manager_application/CMakeFiles/KratosManagerApplication.dir/clean:
	cd /home/enric/kratos/applications/manager_application && $(CMAKE_COMMAND) -P CMakeFiles/KratosManagerApplication.dir/cmake_clean.cmake
.PHONY : applications/manager_application/CMakeFiles/KratosManagerApplication.dir/clean

applications/manager_application/CMakeFiles/KratosManagerApplication.dir/depend:
	cd /home/enric/kratos && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/enric/kratos /home/enric/kratos/applications/manager_application /home/enric/kratos /home/enric/kratos/applications/manager_application /home/enric/kratos/applications/manager_application/CMakeFiles/KratosManagerApplication.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : applications/manager_application/CMakeFiles/KratosManagerApplication.dir/depend

