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
include external_libraries/zlib/CMakeFiles/zlib.dir/depend.make

# Include the progress variables for this target.
include external_libraries/zlib/CMakeFiles/zlib.dir/progress.make

# Include the compile flags for this target's objects.
include external_libraries/zlib/CMakeFiles/zlib.dir/flags.make

external_libraries/zlib/CMakeFiles/zlib.dir/adler32.o: external_libraries/zlib/CMakeFiles/zlib.dir/flags.make
external_libraries/zlib/CMakeFiles/zlib.dir/adler32.o: external_libraries/zlib/adler32.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/enric/kratos/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object external_libraries/zlib/CMakeFiles/zlib.dir/adler32.o"
	cd /home/enric/kratos/external_libraries/zlib && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/zlib.dir/adler32.o   -c /home/enric/kratos/external_libraries/zlib/adler32.c

external_libraries/zlib/CMakeFiles/zlib.dir/adler32.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlib.dir/adler32.i"
	cd /home/enric/kratos/external_libraries/zlib && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -E /home/enric/kratos/external_libraries/zlib/adler32.c > CMakeFiles/zlib.dir/adler32.i

external_libraries/zlib/CMakeFiles/zlib.dir/adler32.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlib.dir/adler32.s"
	cd /home/enric/kratos/external_libraries/zlib && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -S /home/enric/kratos/external_libraries/zlib/adler32.c -o CMakeFiles/zlib.dir/adler32.s

external_libraries/zlib/CMakeFiles/zlib.dir/adler32.o.requires:

.PHONY : external_libraries/zlib/CMakeFiles/zlib.dir/adler32.o.requires

external_libraries/zlib/CMakeFiles/zlib.dir/adler32.o.provides: external_libraries/zlib/CMakeFiles/zlib.dir/adler32.o.requires
	$(MAKE) -f external_libraries/zlib/CMakeFiles/zlib.dir/build.make external_libraries/zlib/CMakeFiles/zlib.dir/adler32.o.provides.build
.PHONY : external_libraries/zlib/CMakeFiles/zlib.dir/adler32.o.provides

external_libraries/zlib/CMakeFiles/zlib.dir/adler32.o.provides.build: external_libraries/zlib/CMakeFiles/zlib.dir/adler32.o


external_libraries/zlib/CMakeFiles/zlib.dir/compress.o: external_libraries/zlib/CMakeFiles/zlib.dir/flags.make
external_libraries/zlib/CMakeFiles/zlib.dir/compress.o: external_libraries/zlib/compress.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/enric/kratos/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object external_libraries/zlib/CMakeFiles/zlib.dir/compress.o"
	cd /home/enric/kratos/external_libraries/zlib && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/zlib.dir/compress.o   -c /home/enric/kratos/external_libraries/zlib/compress.c

external_libraries/zlib/CMakeFiles/zlib.dir/compress.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlib.dir/compress.i"
	cd /home/enric/kratos/external_libraries/zlib && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -E /home/enric/kratos/external_libraries/zlib/compress.c > CMakeFiles/zlib.dir/compress.i

external_libraries/zlib/CMakeFiles/zlib.dir/compress.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlib.dir/compress.s"
	cd /home/enric/kratos/external_libraries/zlib && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -S /home/enric/kratos/external_libraries/zlib/compress.c -o CMakeFiles/zlib.dir/compress.s

external_libraries/zlib/CMakeFiles/zlib.dir/compress.o.requires:

.PHONY : external_libraries/zlib/CMakeFiles/zlib.dir/compress.o.requires

external_libraries/zlib/CMakeFiles/zlib.dir/compress.o.provides: external_libraries/zlib/CMakeFiles/zlib.dir/compress.o.requires
	$(MAKE) -f external_libraries/zlib/CMakeFiles/zlib.dir/build.make external_libraries/zlib/CMakeFiles/zlib.dir/compress.o.provides.build
.PHONY : external_libraries/zlib/CMakeFiles/zlib.dir/compress.o.provides

external_libraries/zlib/CMakeFiles/zlib.dir/compress.o.provides.build: external_libraries/zlib/CMakeFiles/zlib.dir/compress.o


external_libraries/zlib/CMakeFiles/zlib.dir/crc32.o: external_libraries/zlib/CMakeFiles/zlib.dir/flags.make
external_libraries/zlib/CMakeFiles/zlib.dir/crc32.o: external_libraries/zlib/crc32.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/enric/kratos/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object external_libraries/zlib/CMakeFiles/zlib.dir/crc32.o"
	cd /home/enric/kratos/external_libraries/zlib && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/zlib.dir/crc32.o   -c /home/enric/kratos/external_libraries/zlib/crc32.c

external_libraries/zlib/CMakeFiles/zlib.dir/crc32.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlib.dir/crc32.i"
	cd /home/enric/kratos/external_libraries/zlib && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -E /home/enric/kratos/external_libraries/zlib/crc32.c > CMakeFiles/zlib.dir/crc32.i

external_libraries/zlib/CMakeFiles/zlib.dir/crc32.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlib.dir/crc32.s"
	cd /home/enric/kratos/external_libraries/zlib && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -S /home/enric/kratos/external_libraries/zlib/crc32.c -o CMakeFiles/zlib.dir/crc32.s

external_libraries/zlib/CMakeFiles/zlib.dir/crc32.o.requires:

.PHONY : external_libraries/zlib/CMakeFiles/zlib.dir/crc32.o.requires

external_libraries/zlib/CMakeFiles/zlib.dir/crc32.o.provides: external_libraries/zlib/CMakeFiles/zlib.dir/crc32.o.requires
	$(MAKE) -f external_libraries/zlib/CMakeFiles/zlib.dir/build.make external_libraries/zlib/CMakeFiles/zlib.dir/crc32.o.provides.build
.PHONY : external_libraries/zlib/CMakeFiles/zlib.dir/crc32.o.provides

external_libraries/zlib/CMakeFiles/zlib.dir/crc32.o.provides.build: external_libraries/zlib/CMakeFiles/zlib.dir/crc32.o


external_libraries/zlib/CMakeFiles/zlib.dir/deflate.o: external_libraries/zlib/CMakeFiles/zlib.dir/flags.make
external_libraries/zlib/CMakeFiles/zlib.dir/deflate.o: external_libraries/zlib/deflate.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/enric/kratos/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building C object external_libraries/zlib/CMakeFiles/zlib.dir/deflate.o"
	cd /home/enric/kratos/external_libraries/zlib && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/zlib.dir/deflate.o   -c /home/enric/kratos/external_libraries/zlib/deflate.c

external_libraries/zlib/CMakeFiles/zlib.dir/deflate.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlib.dir/deflate.i"
	cd /home/enric/kratos/external_libraries/zlib && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -E /home/enric/kratos/external_libraries/zlib/deflate.c > CMakeFiles/zlib.dir/deflate.i

external_libraries/zlib/CMakeFiles/zlib.dir/deflate.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlib.dir/deflate.s"
	cd /home/enric/kratos/external_libraries/zlib && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -S /home/enric/kratos/external_libraries/zlib/deflate.c -o CMakeFiles/zlib.dir/deflate.s

external_libraries/zlib/CMakeFiles/zlib.dir/deflate.o.requires:

.PHONY : external_libraries/zlib/CMakeFiles/zlib.dir/deflate.o.requires

external_libraries/zlib/CMakeFiles/zlib.dir/deflate.o.provides: external_libraries/zlib/CMakeFiles/zlib.dir/deflate.o.requires
	$(MAKE) -f external_libraries/zlib/CMakeFiles/zlib.dir/build.make external_libraries/zlib/CMakeFiles/zlib.dir/deflate.o.provides.build
.PHONY : external_libraries/zlib/CMakeFiles/zlib.dir/deflate.o.provides

external_libraries/zlib/CMakeFiles/zlib.dir/deflate.o.provides.build: external_libraries/zlib/CMakeFiles/zlib.dir/deflate.o


external_libraries/zlib/CMakeFiles/zlib.dir/gzclose.o: external_libraries/zlib/CMakeFiles/zlib.dir/flags.make
external_libraries/zlib/CMakeFiles/zlib.dir/gzclose.o: external_libraries/zlib/gzclose.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/enric/kratos/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building C object external_libraries/zlib/CMakeFiles/zlib.dir/gzclose.o"
	cd /home/enric/kratos/external_libraries/zlib && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/zlib.dir/gzclose.o   -c /home/enric/kratos/external_libraries/zlib/gzclose.c

external_libraries/zlib/CMakeFiles/zlib.dir/gzclose.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlib.dir/gzclose.i"
	cd /home/enric/kratos/external_libraries/zlib && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -E /home/enric/kratos/external_libraries/zlib/gzclose.c > CMakeFiles/zlib.dir/gzclose.i

external_libraries/zlib/CMakeFiles/zlib.dir/gzclose.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlib.dir/gzclose.s"
	cd /home/enric/kratos/external_libraries/zlib && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -S /home/enric/kratos/external_libraries/zlib/gzclose.c -o CMakeFiles/zlib.dir/gzclose.s

external_libraries/zlib/CMakeFiles/zlib.dir/gzclose.o.requires:

.PHONY : external_libraries/zlib/CMakeFiles/zlib.dir/gzclose.o.requires

external_libraries/zlib/CMakeFiles/zlib.dir/gzclose.o.provides: external_libraries/zlib/CMakeFiles/zlib.dir/gzclose.o.requires
	$(MAKE) -f external_libraries/zlib/CMakeFiles/zlib.dir/build.make external_libraries/zlib/CMakeFiles/zlib.dir/gzclose.o.provides.build
.PHONY : external_libraries/zlib/CMakeFiles/zlib.dir/gzclose.o.provides

external_libraries/zlib/CMakeFiles/zlib.dir/gzclose.o.provides.build: external_libraries/zlib/CMakeFiles/zlib.dir/gzclose.o


external_libraries/zlib/CMakeFiles/zlib.dir/gzlib.o: external_libraries/zlib/CMakeFiles/zlib.dir/flags.make
external_libraries/zlib/CMakeFiles/zlib.dir/gzlib.o: external_libraries/zlib/gzlib.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/enric/kratos/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building C object external_libraries/zlib/CMakeFiles/zlib.dir/gzlib.o"
	cd /home/enric/kratos/external_libraries/zlib && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/zlib.dir/gzlib.o   -c /home/enric/kratos/external_libraries/zlib/gzlib.c

external_libraries/zlib/CMakeFiles/zlib.dir/gzlib.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlib.dir/gzlib.i"
	cd /home/enric/kratos/external_libraries/zlib && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -E /home/enric/kratos/external_libraries/zlib/gzlib.c > CMakeFiles/zlib.dir/gzlib.i

external_libraries/zlib/CMakeFiles/zlib.dir/gzlib.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlib.dir/gzlib.s"
	cd /home/enric/kratos/external_libraries/zlib && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -S /home/enric/kratos/external_libraries/zlib/gzlib.c -o CMakeFiles/zlib.dir/gzlib.s

external_libraries/zlib/CMakeFiles/zlib.dir/gzlib.o.requires:

.PHONY : external_libraries/zlib/CMakeFiles/zlib.dir/gzlib.o.requires

external_libraries/zlib/CMakeFiles/zlib.dir/gzlib.o.provides: external_libraries/zlib/CMakeFiles/zlib.dir/gzlib.o.requires
	$(MAKE) -f external_libraries/zlib/CMakeFiles/zlib.dir/build.make external_libraries/zlib/CMakeFiles/zlib.dir/gzlib.o.provides.build
.PHONY : external_libraries/zlib/CMakeFiles/zlib.dir/gzlib.o.provides

external_libraries/zlib/CMakeFiles/zlib.dir/gzlib.o.provides.build: external_libraries/zlib/CMakeFiles/zlib.dir/gzlib.o


external_libraries/zlib/CMakeFiles/zlib.dir/gzread.o: external_libraries/zlib/CMakeFiles/zlib.dir/flags.make
external_libraries/zlib/CMakeFiles/zlib.dir/gzread.o: external_libraries/zlib/gzread.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/enric/kratos/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building C object external_libraries/zlib/CMakeFiles/zlib.dir/gzread.o"
	cd /home/enric/kratos/external_libraries/zlib && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/zlib.dir/gzread.o   -c /home/enric/kratos/external_libraries/zlib/gzread.c

external_libraries/zlib/CMakeFiles/zlib.dir/gzread.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlib.dir/gzread.i"
	cd /home/enric/kratos/external_libraries/zlib && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -E /home/enric/kratos/external_libraries/zlib/gzread.c > CMakeFiles/zlib.dir/gzread.i

external_libraries/zlib/CMakeFiles/zlib.dir/gzread.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlib.dir/gzread.s"
	cd /home/enric/kratos/external_libraries/zlib && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -S /home/enric/kratos/external_libraries/zlib/gzread.c -o CMakeFiles/zlib.dir/gzread.s

external_libraries/zlib/CMakeFiles/zlib.dir/gzread.o.requires:

.PHONY : external_libraries/zlib/CMakeFiles/zlib.dir/gzread.o.requires

external_libraries/zlib/CMakeFiles/zlib.dir/gzread.o.provides: external_libraries/zlib/CMakeFiles/zlib.dir/gzread.o.requires
	$(MAKE) -f external_libraries/zlib/CMakeFiles/zlib.dir/build.make external_libraries/zlib/CMakeFiles/zlib.dir/gzread.o.provides.build
.PHONY : external_libraries/zlib/CMakeFiles/zlib.dir/gzread.o.provides

external_libraries/zlib/CMakeFiles/zlib.dir/gzread.o.provides.build: external_libraries/zlib/CMakeFiles/zlib.dir/gzread.o


external_libraries/zlib/CMakeFiles/zlib.dir/gzwrite.o: external_libraries/zlib/CMakeFiles/zlib.dir/flags.make
external_libraries/zlib/CMakeFiles/zlib.dir/gzwrite.o: external_libraries/zlib/gzwrite.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/enric/kratos/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building C object external_libraries/zlib/CMakeFiles/zlib.dir/gzwrite.o"
	cd /home/enric/kratos/external_libraries/zlib && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/zlib.dir/gzwrite.o   -c /home/enric/kratos/external_libraries/zlib/gzwrite.c

external_libraries/zlib/CMakeFiles/zlib.dir/gzwrite.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlib.dir/gzwrite.i"
	cd /home/enric/kratos/external_libraries/zlib && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -E /home/enric/kratos/external_libraries/zlib/gzwrite.c > CMakeFiles/zlib.dir/gzwrite.i

external_libraries/zlib/CMakeFiles/zlib.dir/gzwrite.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlib.dir/gzwrite.s"
	cd /home/enric/kratos/external_libraries/zlib && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -S /home/enric/kratos/external_libraries/zlib/gzwrite.c -o CMakeFiles/zlib.dir/gzwrite.s

external_libraries/zlib/CMakeFiles/zlib.dir/gzwrite.o.requires:

.PHONY : external_libraries/zlib/CMakeFiles/zlib.dir/gzwrite.o.requires

external_libraries/zlib/CMakeFiles/zlib.dir/gzwrite.o.provides: external_libraries/zlib/CMakeFiles/zlib.dir/gzwrite.o.requires
	$(MAKE) -f external_libraries/zlib/CMakeFiles/zlib.dir/build.make external_libraries/zlib/CMakeFiles/zlib.dir/gzwrite.o.provides.build
.PHONY : external_libraries/zlib/CMakeFiles/zlib.dir/gzwrite.o.provides

external_libraries/zlib/CMakeFiles/zlib.dir/gzwrite.o.provides.build: external_libraries/zlib/CMakeFiles/zlib.dir/gzwrite.o


external_libraries/zlib/CMakeFiles/zlib.dir/inflate.o: external_libraries/zlib/CMakeFiles/zlib.dir/flags.make
external_libraries/zlib/CMakeFiles/zlib.dir/inflate.o: external_libraries/zlib/inflate.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/enric/kratos/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building C object external_libraries/zlib/CMakeFiles/zlib.dir/inflate.o"
	cd /home/enric/kratos/external_libraries/zlib && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/zlib.dir/inflate.o   -c /home/enric/kratos/external_libraries/zlib/inflate.c

external_libraries/zlib/CMakeFiles/zlib.dir/inflate.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlib.dir/inflate.i"
	cd /home/enric/kratos/external_libraries/zlib && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -E /home/enric/kratos/external_libraries/zlib/inflate.c > CMakeFiles/zlib.dir/inflate.i

external_libraries/zlib/CMakeFiles/zlib.dir/inflate.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlib.dir/inflate.s"
	cd /home/enric/kratos/external_libraries/zlib && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -S /home/enric/kratos/external_libraries/zlib/inflate.c -o CMakeFiles/zlib.dir/inflate.s

external_libraries/zlib/CMakeFiles/zlib.dir/inflate.o.requires:

.PHONY : external_libraries/zlib/CMakeFiles/zlib.dir/inflate.o.requires

external_libraries/zlib/CMakeFiles/zlib.dir/inflate.o.provides: external_libraries/zlib/CMakeFiles/zlib.dir/inflate.o.requires
	$(MAKE) -f external_libraries/zlib/CMakeFiles/zlib.dir/build.make external_libraries/zlib/CMakeFiles/zlib.dir/inflate.o.provides.build
.PHONY : external_libraries/zlib/CMakeFiles/zlib.dir/inflate.o.provides

external_libraries/zlib/CMakeFiles/zlib.dir/inflate.o.provides.build: external_libraries/zlib/CMakeFiles/zlib.dir/inflate.o


external_libraries/zlib/CMakeFiles/zlib.dir/infback.o: external_libraries/zlib/CMakeFiles/zlib.dir/flags.make
external_libraries/zlib/CMakeFiles/zlib.dir/infback.o: external_libraries/zlib/infback.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/enric/kratos/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building C object external_libraries/zlib/CMakeFiles/zlib.dir/infback.o"
	cd /home/enric/kratos/external_libraries/zlib && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/zlib.dir/infback.o   -c /home/enric/kratos/external_libraries/zlib/infback.c

external_libraries/zlib/CMakeFiles/zlib.dir/infback.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlib.dir/infback.i"
	cd /home/enric/kratos/external_libraries/zlib && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -E /home/enric/kratos/external_libraries/zlib/infback.c > CMakeFiles/zlib.dir/infback.i

external_libraries/zlib/CMakeFiles/zlib.dir/infback.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlib.dir/infback.s"
	cd /home/enric/kratos/external_libraries/zlib && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -S /home/enric/kratos/external_libraries/zlib/infback.c -o CMakeFiles/zlib.dir/infback.s

external_libraries/zlib/CMakeFiles/zlib.dir/infback.o.requires:

.PHONY : external_libraries/zlib/CMakeFiles/zlib.dir/infback.o.requires

external_libraries/zlib/CMakeFiles/zlib.dir/infback.o.provides: external_libraries/zlib/CMakeFiles/zlib.dir/infback.o.requires
	$(MAKE) -f external_libraries/zlib/CMakeFiles/zlib.dir/build.make external_libraries/zlib/CMakeFiles/zlib.dir/infback.o.provides.build
.PHONY : external_libraries/zlib/CMakeFiles/zlib.dir/infback.o.provides

external_libraries/zlib/CMakeFiles/zlib.dir/infback.o.provides.build: external_libraries/zlib/CMakeFiles/zlib.dir/infback.o


external_libraries/zlib/CMakeFiles/zlib.dir/inftrees.o: external_libraries/zlib/CMakeFiles/zlib.dir/flags.make
external_libraries/zlib/CMakeFiles/zlib.dir/inftrees.o: external_libraries/zlib/inftrees.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/enric/kratos/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building C object external_libraries/zlib/CMakeFiles/zlib.dir/inftrees.o"
	cd /home/enric/kratos/external_libraries/zlib && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/zlib.dir/inftrees.o   -c /home/enric/kratos/external_libraries/zlib/inftrees.c

external_libraries/zlib/CMakeFiles/zlib.dir/inftrees.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlib.dir/inftrees.i"
	cd /home/enric/kratos/external_libraries/zlib && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -E /home/enric/kratos/external_libraries/zlib/inftrees.c > CMakeFiles/zlib.dir/inftrees.i

external_libraries/zlib/CMakeFiles/zlib.dir/inftrees.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlib.dir/inftrees.s"
	cd /home/enric/kratos/external_libraries/zlib && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -S /home/enric/kratos/external_libraries/zlib/inftrees.c -o CMakeFiles/zlib.dir/inftrees.s

external_libraries/zlib/CMakeFiles/zlib.dir/inftrees.o.requires:

.PHONY : external_libraries/zlib/CMakeFiles/zlib.dir/inftrees.o.requires

external_libraries/zlib/CMakeFiles/zlib.dir/inftrees.o.provides: external_libraries/zlib/CMakeFiles/zlib.dir/inftrees.o.requires
	$(MAKE) -f external_libraries/zlib/CMakeFiles/zlib.dir/build.make external_libraries/zlib/CMakeFiles/zlib.dir/inftrees.o.provides.build
.PHONY : external_libraries/zlib/CMakeFiles/zlib.dir/inftrees.o.provides

external_libraries/zlib/CMakeFiles/zlib.dir/inftrees.o.provides.build: external_libraries/zlib/CMakeFiles/zlib.dir/inftrees.o


external_libraries/zlib/CMakeFiles/zlib.dir/inffast.o: external_libraries/zlib/CMakeFiles/zlib.dir/flags.make
external_libraries/zlib/CMakeFiles/zlib.dir/inffast.o: external_libraries/zlib/inffast.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/enric/kratos/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Building C object external_libraries/zlib/CMakeFiles/zlib.dir/inffast.o"
	cd /home/enric/kratos/external_libraries/zlib && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/zlib.dir/inffast.o   -c /home/enric/kratos/external_libraries/zlib/inffast.c

external_libraries/zlib/CMakeFiles/zlib.dir/inffast.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlib.dir/inffast.i"
	cd /home/enric/kratos/external_libraries/zlib && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -E /home/enric/kratos/external_libraries/zlib/inffast.c > CMakeFiles/zlib.dir/inffast.i

external_libraries/zlib/CMakeFiles/zlib.dir/inffast.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlib.dir/inffast.s"
	cd /home/enric/kratos/external_libraries/zlib && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -S /home/enric/kratos/external_libraries/zlib/inffast.c -o CMakeFiles/zlib.dir/inffast.s

external_libraries/zlib/CMakeFiles/zlib.dir/inffast.o.requires:

.PHONY : external_libraries/zlib/CMakeFiles/zlib.dir/inffast.o.requires

external_libraries/zlib/CMakeFiles/zlib.dir/inffast.o.provides: external_libraries/zlib/CMakeFiles/zlib.dir/inffast.o.requires
	$(MAKE) -f external_libraries/zlib/CMakeFiles/zlib.dir/build.make external_libraries/zlib/CMakeFiles/zlib.dir/inffast.o.provides.build
.PHONY : external_libraries/zlib/CMakeFiles/zlib.dir/inffast.o.provides

external_libraries/zlib/CMakeFiles/zlib.dir/inffast.o.provides.build: external_libraries/zlib/CMakeFiles/zlib.dir/inffast.o


external_libraries/zlib/CMakeFiles/zlib.dir/trees.o: external_libraries/zlib/CMakeFiles/zlib.dir/flags.make
external_libraries/zlib/CMakeFiles/zlib.dir/trees.o: external_libraries/zlib/trees.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/enric/kratos/CMakeFiles --progress-num=$(CMAKE_PROGRESS_13) "Building C object external_libraries/zlib/CMakeFiles/zlib.dir/trees.o"
	cd /home/enric/kratos/external_libraries/zlib && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/zlib.dir/trees.o   -c /home/enric/kratos/external_libraries/zlib/trees.c

external_libraries/zlib/CMakeFiles/zlib.dir/trees.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlib.dir/trees.i"
	cd /home/enric/kratos/external_libraries/zlib && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -E /home/enric/kratos/external_libraries/zlib/trees.c > CMakeFiles/zlib.dir/trees.i

external_libraries/zlib/CMakeFiles/zlib.dir/trees.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlib.dir/trees.s"
	cd /home/enric/kratos/external_libraries/zlib && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -S /home/enric/kratos/external_libraries/zlib/trees.c -o CMakeFiles/zlib.dir/trees.s

external_libraries/zlib/CMakeFiles/zlib.dir/trees.o.requires:

.PHONY : external_libraries/zlib/CMakeFiles/zlib.dir/trees.o.requires

external_libraries/zlib/CMakeFiles/zlib.dir/trees.o.provides: external_libraries/zlib/CMakeFiles/zlib.dir/trees.o.requires
	$(MAKE) -f external_libraries/zlib/CMakeFiles/zlib.dir/build.make external_libraries/zlib/CMakeFiles/zlib.dir/trees.o.provides.build
.PHONY : external_libraries/zlib/CMakeFiles/zlib.dir/trees.o.provides

external_libraries/zlib/CMakeFiles/zlib.dir/trees.o.provides.build: external_libraries/zlib/CMakeFiles/zlib.dir/trees.o


external_libraries/zlib/CMakeFiles/zlib.dir/uncompr.o: external_libraries/zlib/CMakeFiles/zlib.dir/flags.make
external_libraries/zlib/CMakeFiles/zlib.dir/uncompr.o: external_libraries/zlib/uncompr.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/enric/kratos/CMakeFiles --progress-num=$(CMAKE_PROGRESS_14) "Building C object external_libraries/zlib/CMakeFiles/zlib.dir/uncompr.o"
	cd /home/enric/kratos/external_libraries/zlib && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/zlib.dir/uncompr.o   -c /home/enric/kratos/external_libraries/zlib/uncompr.c

external_libraries/zlib/CMakeFiles/zlib.dir/uncompr.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlib.dir/uncompr.i"
	cd /home/enric/kratos/external_libraries/zlib && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -E /home/enric/kratos/external_libraries/zlib/uncompr.c > CMakeFiles/zlib.dir/uncompr.i

external_libraries/zlib/CMakeFiles/zlib.dir/uncompr.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlib.dir/uncompr.s"
	cd /home/enric/kratos/external_libraries/zlib && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -S /home/enric/kratos/external_libraries/zlib/uncompr.c -o CMakeFiles/zlib.dir/uncompr.s

external_libraries/zlib/CMakeFiles/zlib.dir/uncompr.o.requires:

.PHONY : external_libraries/zlib/CMakeFiles/zlib.dir/uncompr.o.requires

external_libraries/zlib/CMakeFiles/zlib.dir/uncompr.o.provides: external_libraries/zlib/CMakeFiles/zlib.dir/uncompr.o.requires
	$(MAKE) -f external_libraries/zlib/CMakeFiles/zlib.dir/build.make external_libraries/zlib/CMakeFiles/zlib.dir/uncompr.o.provides.build
.PHONY : external_libraries/zlib/CMakeFiles/zlib.dir/uncompr.o.provides

external_libraries/zlib/CMakeFiles/zlib.dir/uncompr.o.provides.build: external_libraries/zlib/CMakeFiles/zlib.dir/uncompr.o


external_libraries/zlib/CMakeFiles/zlib.dir/zutil.o: external_libraries/zlib/CMakeFiles/zlib.dir/flags.make
external_libraries/zlib/CMakeFiles/zlib.dir/zutil.o: external_libraries/zlib/zutil.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/enric/kratos/CMakeFiles --progress-num=$(CMAKE_PROGRESS_15) "Building C object external_libraries/zlib/CMakeFiles/zlib.dir/zutil.o"
	cd /home/enric/kratos/external_libraries/zlib && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/zlib.dir/zutil.o   -c /home/enric/kratos/external_libraries/zlib/zutil.c

external_libraries/zlib/CMakeFiles/zlib.dir/zutil.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/zlib.dir/zutil.i"
	cd /home/enric/kratos/external_libraries/zlib && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -E /home/enric/kratos/external_libraries/zlib/zutil.c > CMakeFiles/zlib.dir/zutil.i

external_libraries/zlib/CMakeFiles/zlib.dir/zutil.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/zlib.dir/zutil.s"
	cd /home/enric/kratos/external_libraries/zlib && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -S /home/enric/kratos/external_libraries/zlib/zutil.c -o CMakeFiles/zlib.dir/zutil.s

external_libraries/zlib/CMakeFiles/zlib.dir/zutil.o.requires:

.PHONY : external_libraries/zlib/CMakeFiles/zlib.dir/zutil.o.requires

external_libraries/zlib/CMakeFiles/zlib.dir/zutil.o.provides: external_libraries/zlib/CMakeFiles/zlib.dir/zutil.o.requires
	$(MAKE) -f external_libraries/zlib/CMakeFiles/zlib.dir/build.make external_libraries/zlib/CMakeFiles/zlib.dir/zutil.o.provides.build
.PHONY : external_libraries/zlib/CMakeFiles/zlib.dir/zutil.o.provides

external_libraries/zlib/CMakeFiles/zlib.dir/zutil.o.provides.build: external_libraries/zlib/CMakeFiles/zlib.dir/zutil.o


# Object files for target zlib
zlib_OBJECTS = \
"CMakeFiles/zlib.dir/adler32.o" \
"CMakeFiles/zlib.dir/compress.o" \
"CMakeFiles/zlib.dir/crc32.o" \
"CMakeFiles/zlib.dir/deflate.o" \
"CMakeFiles/zlib.dir/gzclose.o" \
"CMakeFiles/zlib.dir/gzlib.o" \
"CMakeFiles/zlib.dir/gzread.o" \
"CMakeFiles/zlib.dir/gzwrite.o" \
"CMakeFiles/zlib.dir/inflate.o" \
"CMakeFiles/zlib.dir/infback.o" \
"CMakeFiles/zlib.dir/inftrees.o" \
"CMakeFiles/zlib.dir/inffast.o" \
"CMakeFiles/zlib.dir/trees.o" \
"CMakeFiles/zlib.dir/uncompr.o" \
"CMakeFiles/zlib.dir/zutil.o"

# External object files for target zlib
zlib_EXTERNAL_OBJECTS =

external_libraries/zlib/libz.so.1.2.8: external_libraries/zlib/CMakeFiles/zlib.dir/adler32.o
external_libraries/zlib/libz.so.1.2.8: external_libraries/zlib/CMakeFiles/zlib.dir/compress.o
external_libraries/zlib/libz.so.1.2.8: external_libraries/zlib/CMakeFiles/zlib.dir/crc32.o
external_libraries/zlib/libz.so.1.2.8: external_libraries/zlib/CMakeFiles/zlib.dir/deflate.o
external_libraries/zlib/libz.so.1.2.8: external_libraries/zlib/CMakeFiles/zlib.dir/gzclose.o
external_libraries/zlib/libz.so.1.2.8: external_libraries/zlib/CMakeFiles/zlib.dir/gzlib.o
external_libraries/zlib/libz.so.1.2.8: external_libraries/zlib/CMakeFiles/zlib.dir/gzread.o
external_libraries/zlib/libz.so.1.2.8: external_libraries/zlib/CMakeFiles/zlib.dir/gzwrite.o
external_libraries/zlib/libz.so.1.2.8: external_libraries/zlib/CMakeFiles/zlib.dir/inflate.o
external_libraries/zlib/libz.so.1.2.8: external_libraries/zlib/CMakeFiles/zlib.dir/infback.o
external_libraries/zlib/libz.so.1.2.8: external_libraries/zlib/CMakeFiles/zlib.dir/inftrees.o
external_libraries/zlib/libz.so.1.2.8: external_libraries/zlib/CMakeFiles/zlib.dir/inffast.o
external_libraries/zlib/libz.so.1.2.8: external_libraries/zlib/CMakeFiles/zlib.dir/trees.o
external_libraries/zlib/libz.so.1.2.8: external_libraries/zlib/CMakeFiles/zlib.dir/uncompr.o
external_libraries/zlib/libz.so.1.2.8: external_libraries/zlib/CMakeFiles/zlib.dir/zutil.o
external_libraries/zlib/libz.so.1.2.8: external_libraries/zlib/CMakeFiles/zlib.dir/build.make
external_libraries/zlib/libz.so.1.2.8: external_libraries/zlib/CMakeFiles/zlib.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/enric/kratos/CMakeFiles --progress-num=$(CMAKE_PROGRESS_16) "Linking C shared library libz.so"
	cd /home/enric/kratos/external_libraries/zlib && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/zlib.dir/link.txt --verbose=$(VERBOSE)
	cd /home/enric/kratos/external_libraries/zlib && $(CMAKE_COMMAND) -E cmake_symlink_library libz.so.1.2.8 libz.so.1 libz.so

external_libraries/zlib/libz.so.1: external_libraries/zlib/libz.so.1.2.8
	@$(CMAKE_COMMAND) -E touch_nocreate external_libraries/zlib/libz.so.1

external_libraries/zlib/libz.so: external_libraries/zlib/libz.so.1.2.8
	@$(CMAKE_COMMAND) -E touch_nocreate external_libraries/zlib/libz.so

# Rule to build all files generated by this target.
external_libraries/zlib/CMakeFiles/zlib.dir/build: external_libraries/zlib/libz.so

.PHONY : external_libraries/zlib/CMakeFiles/zlib.dir/build

external_libraries/zlib/CMakeFiles/zlib.dir/requires: external_libraries/zlib/CMakeFiles/zlib.dir/adler32.o.requires
external_libraries/zlib/CMakeFiles/zlib.dir/requires: external_libraries/zlib/CMakeFiles/zlib.dir/compress.o.requires
external_libraries/zlib/CMakeFiles/zlib.dir/requires: external_libraries/zlib/CMakeFiles/zlib.dir/crc32.o.requires
external_libraries/zlib/CMakeFiles/zlib.dir/requires: external_libraries/zlib/CMakeFiles/zlib.dir/deflate.o.requires
external_libraries/zlib/CMakeFiles/zlib.dir/requires: external_libraries/zlib/CMakeFiles/zlib.dir/gzclose.o.requires
external_libraries/zlib/CMakeFiles/zlib.dir/requires: external_libraries/zlib/CMakeFiles/zlib.dir/gzlib.o.requires
external_libraries/zlib/CMakeFiles/zlib.dir/requires: external_libraries/zlib/CMakeFiles/zlib.dir/gzread.o.requires
external_libraries/zlib/CMakeFiles/zlib.dir/requires: external_libraries/zlib/CMakeFiles/zlib.dir/gzwrite.o.requires
external_libraries/zlib/CMakeFiles/zlib.dir/requires: external_libraries/zlib/CMakeFiles/zlib.dir/inflate.o.requires
external_libraries/zlib/CMakeFiles/zlib.dir/requires: external_libraries/zlib/CMakeFiles/zlib.dir/infback.o.requires
external_libraries/zlib/CMakeFiles/zlib.dir/requires: external_libraries/zlib/CMakeFiles/zlib.dir/inftrees.o.requires
external_libraries/zlib/CMakeFiles/zlib.dir/requires: external_libraries/zlib/CMakeFiles/zlib.dir/inffast.o.requires
external_libraries/zlib/CMakeFiles/zlib.dir/requires: external_libraries/zlib/CMakeFiles/zlib.dir/trees.o.requires
external_libraries/zlib/CMakeFiles/zlib.dir/requires: external_libraries/zlib/CMakeFiles/zlib.dir/uncompr.o.requires
external_libraries/zlib/CMakeFiles/zlib.dir/requires: external_libraries/zlib/CMakeFiles/zlib.dir/zutil.o.requires

.PHONY : external_libraries/zlib/CMakeFiles/zlib.dir/requires

external_libraries/zlib/CMakeFiles/zlib.dir/clean:
	cd /home/enric/kratos/external_libraries/zlib && $(CMAKE_COMMAND) -P CMakeFiles/zlib.dir/cmake_clean.cmake
.PHONY : external_libraries/zlib/CMakeFiles/zlib.dir/clean

external_libraries/zlib/CMakeFiles/zlib.dir/depend:
	cd /home/enric/kratos && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/enric/kratos /home/enric/kratos/external_libraries/zlib /home/enric/kratos /home/enric/kratos/external_libraries/zlib /home/enric/kratos/external_libraries/zlib/CMakeFiles/zlib.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : external_libraries/zlib/CMakeFiles/zlib.dir/depend
