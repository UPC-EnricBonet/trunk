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
include external_libraries/zlib/CMakeFiles/example64.dir/depend.make

# Include the progress variables for this target.
include external_libraries/zlib/CMakeFiles/example64.dir/progress.make

# Include the compile flags for this target's objects.
include external_libraries/zlib/CMakeFiles/example64.dir/flags.make

external_libraries/zlib/CMakeFiles/example64.dir/test/example.o: external_libraries/zlib/CMakeFiles/example64.dir/flags.make
external_libraries/zlib/CMakeFiles/example64.dir/test/example.o: external_libraries/zlib/test/example.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/enric/kratos/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object external_libraries/zlib/CMakeFiles/example64.dir/test/example.o"
	cd /home/enric/kratos/external_libraries/zlib && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/example64.dir/test/example.o   -c /home/enric/kratos/external_libraries/zlib/test/example.c

external_libraries/zlib/CMakeFiles/example64.dir/test/example.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/example64.dir/test/example.i"
	cd /home/enric/kratos/external_libraries/zlib && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -E /home/enric/kratos/external_libraries/zlib/test/example.c > CMakeFiles/example64.dir/test/example.i

external_libraries/zlib/CMakeFiles/example64.dir/test/example.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/example64.dir/test/example.s"
	cd /home/enric/kratos/external_libraries/zlib && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -S /home/enric/kratos/external_libraries/zlib/test/example.c -o CMakeFiles/example64.dir/test/example.s

external_libraries/zlib/CMakeFiles/example64.dir/test/example.o.requires:

.PHONY : external_libraries/zlib/CMakeFiles/example64.dir/test/example.o.requires

external_libraries/zlib/CMakeFiles/example64.dir/test/example.o.provides: external_libraries/zlib/CMakeFiles/example64.dir/test/example.o.requires
	$(MAKE) -f external_libraries/zlib/CMakeFiles/example64.dir/build.make external_libraries/zlib/CMakeFiles/example64.dir/test/example.o.provides.build
.PHONY : external_libraries/zlib/CMakeFiles/example64.dir/test/example.o.provides

external_libraries/zlib/CMakeFiles/example64.dir/test/example.o.provides.build: external_libraries/zlib/CMakeFiles/example64.dir/test/example.o


# Object files for target example64
example64_OBJECTS = \
"CMakeFiles/example64.dir/test/example.o"

# External object files for target example64
example64_EXTERNAL_OBJECTS =

external_libraries/zlib/example64: external_libraries/zlib/CMakeFiles/example64.dir/test/example.o
external_libraries/zlib/example64: external_libraries/zlib/CMakeFiles/example64.dir/build.make
external_libraries/zlib/example64: external_libraries/zlib/libz.so.1.2.8
external_libraries/zlib/example64: external_libraries/zlib/CMakeFiles/example64.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/enric/kratos/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable example64"
	cd /home/enric/kratos/external_libraries/zlib && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/example64.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
external_libraries/zlib/CMakeFiles/example64.dir/build: external_libraries/zlib/example64

.PHONY : external_libraries/zlib/CMakeFiles/example64.dir/build

external_libraries/zlib/CMakeFiles/example64.dir/requires: external_libraries/zlib/CMakeFiles/example64.dir/test/example.o.requires

.PHONY : external_libraries/zlib/CMakeFiles/example64.dir/requires

external_libraries/zlib/CMakeFiles/example64.dir/clean:
	cd /home/enric/kratos/external_libraries/zlib && $(CMAKE_COMMAND) -P CMakeFiles/example64.dir/cmake_clean.cmake
.PHONY : external_libraries/zlib/CMakeFiles/example64.dir/clean

external_libraries/zlib/CMakeFiles/example64.dir/depend:
	cd /home/enric/kratos && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/enric/kratos /home/enric/kratos/external_libraries/zlib /home/enric/kratos /home/enric/kratos/external_libraries/zlib /home/enric/kratos/external_libraries/zlib/CMakeFiles/example64.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : external_libraries/zlib/CMakeFiles/example64.dir/depend
