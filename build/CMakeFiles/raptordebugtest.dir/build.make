# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/toanhoi/Desktop/RaptorDebugTest

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/toanhoi/Desktop/RaptorDebugTest/build

# Include any dependencies generated for this target.
include CMakeFiles/raptordebugtest.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/raptordebugtest.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/raptordebugtest.dir/flags.make

CMakeFiles/raptordebugtest.dir/main.cpp.o: CMakeFiles/raptordebugtest.dir/flags.make
CMakeFiles/raptordebugtest.dir/main.cpp.o: ../main.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/toanhoi/Desktop/RaptorDebugTest/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/raptordebugtest.dir/main.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/raptordebugtest.dir/main.cpp.o -c /home/toanhoi/Desktop/RaptorDebugTest/main.cpp

CMakeFiles/raptordebugtest.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/raptordebugtest.dir/main.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/toanhoi/Desktop/RaptorDebugTest/main.cpp > CMakeFiles/raptordebugtest.dir/main.cpp.i

CMakeFiles/raptordebugtest.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/raptordebugtest.dir/main.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/toanhoi/Desktop/RaptorDebugTest/main.cpp -o CMakeFiles/raptordebugtest.dir/main.cpp.s

CMakeFiles/raptordebugtest.dir/main.cpp.o.requires:
.PHONY : CMakeFiles/raptordebugtest.dir/main.cpp.o.requires

CMakeFiles/raptordebugtest.dir/main.cpp.o.provides: CMakeFiles/raptordebugtest.dir/main.cpp.o.requires
	$(MAKE) -f CMakeFiles/raptordebugtest.dir/build.make CMakeFiles/raptordebugtest.dir/main.cpp.o.provides.build
.PHONY : CMakeFiles/raptordebugtest.dir/main.cpp.o.provides

CMakeFiles/raptordebugtest.dir/main.cpp.o.provides.build: CMakeFiles/raptordebugtest.dir/main.cpp.o

# Object files for target raptordebugtest
raptordebugtest_OBJECTS = \
"CMakeFiles/raptordebugtest.dir/main.cpp.o"

# External object files for target raptordebugtest
raptordebugtest_EXTERNAL_OBJECTS =

raptordebugtest: CMakeFiles/raptordebugtest.dir/main.cpp.o
raptordebugtest: CMakeFiles/raptordebugtest.dir/build.make
raptordebugtest: CMakeFiles/raptordebugtest.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable raptordebugtest"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/raptordebugtest.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/raptordebugtest.dir/build: raptordebugtest
.PHONY : CMakeFiles/raptordebugtest.dir/build

CMakeFiles/raptordebugtest.dir/requires: CMakeFiles/raptordebugtest.dir/main.cpp.o.requires
.PHONY : CMakeFiles/raptordebugtest.dir/requires

CMakeFiles/raptordebugtest.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/raptordebugtest.dir/cmake_clean.cmake
.PHONY : CMakeFiles/raptordebugtest.dir/clean

CMakeFiles/raptordebugtest.dir/depend:
	cd /home/toanhoi/Desktop/RaptorDebugTest/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/toanhoi/Desktop/RaptorDebugTest /home/toanhoi/Desktop/RaptorDebugTest /home/toanhoi/Desktop/RaptorDebugTest/build /home/toanhoi/Desktop/RaptorDebugTest/build /home/toanhoi/Desktop/RaptorDebugTest/build/CMakeFiles/raptordebugtest.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/raptordebugtest.dir/depend

