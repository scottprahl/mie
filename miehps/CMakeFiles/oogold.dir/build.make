# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.7

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
CMAKE_COMMAND = /opt/local/bin/cmake

# The command to remove a file.
RM = /opt/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/prahl/Desktop/miehps-1.0.1

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/prahl/Desktop/miehps-1.0.1

# Include any dependencies generated for this target.
include CMakeFiles/oogold.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/oogold.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/oogold.dir/flags.make

CMakeFiles/oogold.dir/oogold.cpp.o: CMakeFiles/oogold.dir/flags.make
CMakeFiles/oogold.dir/oogold.cpp.o: oogold.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/prahl/Desktop/miehps-1.0.1/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/oogold.dir/oogold.cpp.o"
	g++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/oogold.dir/oogold.cpp.o -c /Users/prahl/Desktop/miehps-1.0.1/oogold.cpp

CMakeFiles/oogold.dir/oogold.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/oogold.dir/oogold.cpp.i"
	g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/prahl/Desktop/miehps-1.0.1/oogold.cpp > CMakeFiles/oogold.dir/oogold.cpp.i

CMakeFiles/oogold.dir/oogold.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/oogold.dir/oogold.cpp.s"
	g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/prahl/Desktop/miehps-1.0.1/oogold.cpp -o CMakeFiles/oogold.dir/oogold.cpp.s

CMakeFiles/oogold.dir/oogold.cpp.o.requires:

.PHONY : CMakeFiles/oogold.dir/oogold.cpp.o.requires

CMakeFiles/oogold.dir/oogold.cpp.o.provides: CMakeFiles/oogold.dir/oogold.cpp.o.requires
	$(MAKE) -f CMakeFiles/oogold.dir/build.make CMakeFiles/oogold.dir/oogold.cpp.o.provides.build
.PHONY : CMakeFiles/oogold.dir/oogold.cpp.o.provides

CMakeFiles/oogold.dir/oogold.cpp.o.provides.build: CMakeFiles/oogold.dir/oogold.cpp.o


CMakeFiles/oogold.dir/miehps.cpp.o: CMakeFiles/oogold.dir/flags.make
CMakeFiles/oogold.dir/miehps.cpp.o: miehps.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/prahl/Desktop/miehps-1.0.1/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/oogold.dir/miehps.cpp.o"
	g++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/oogold.dir/miehps.cpp.o -c /Users/prahl/Desktop/miehps-1.0.1/miehps.cpp

CMakeFiles/oogold.dir/miehps.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/oogold.dir/miehps.cpp.i"
	g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/prahl/Desktop/miehps-1.0.1/miehps.cpp > CMakeFiles/oogold.dir/miehps.cpp.i

CMakeFiles/oogold.dir/miehps.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/oogold.dir/miehps.cpp.s"
	g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/prahl/Desktop/miehps-1.0.1/miehps.cpp -o CMakeFiles/oogold.dir/miehps.cpp.s

CMakeFiles/oogold.dir/miehps.cpp.o.requires:

.PHONY : CMakeFiles/oogold.dir/miehps.cpp.o.requires

CMakeFiles/oogold.dir/miehps.cpp.o.provides: CMakeFiles/oogold.dir/miehps.cpp.o.requires
	$(MAKE) -f CMakeFiles/oogold.dir/build.make CMakeFiles/oogold.dir/miehps.cpp.o.provides.build
.PHONY : CMakeFiles/oogold.dir/miehps.cpp.o.provides

CMakeFiles/oogold.dir/miehps.cpp.o.provides.build: CMakeFiles/oogold.dir/miehps.cpp.o


# Object files for target oogold
oogold_OBJECTS = \
"CMakeFiles/oogold.dir/oogold.cpp.o" \
"CMakeFiles/oogold.dir/miehps.cpp.o"

# External object files for target oogold
oogold_EXTERNAL_OBJECTS =

oogold: CMakeFiles/oogold.dir/oogold.cpp.o
oogold: CMakeFiles/oogold.dir/miehps.cpp.o
oogold: CMakeFiles/oogold.dir/build.make
oogold: CMakeFiles/oogold.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/prahl/Desktop/miehps-1.0.1/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable oogold"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/oogold.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/oogold.dir/build: oogold

.PHONY : CMakeFiles/oogold.dir/build

CMakeFiles/oogold.dir/requires: CMakeFiles/oogold.dir/oogold.cpp.o.requires
CMakeFiles/oogold.dir/requires: CMakeFiles/oogold.dir/miehps.cpp.o.requires

.PHONY : CMakeFiles/oogold.dir/requires

CMakeFiles/oogold.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/oogold.dir/cmake_clean.cmake
.PHONY : CMakeFiles/oogold.dir/clean

CMakeFiles/oogold.dir/depend:
	cd /Users/prahl/Desktop/miehps-1.0.1 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/prahl/Desktop/miehps-1.0.1 /Users/prahl/Desktop/miehps-1.0.1 /Users/prahl/Desktop/miehps-1.0.1 /Users/prahl/Desktop/miehps-1.0.1 /Users/prahl/Desktop/miehps-1.0.1/CMakeFiles/oogold.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/oogold.dir/depend
