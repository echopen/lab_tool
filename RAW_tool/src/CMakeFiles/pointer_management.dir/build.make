# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/lecoucl/Projects/EchOpen/lab_tool/RAW_tool

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/lecoucl/Projects/EchOpen/lab_tool/RAW_tool

# Include any dependencies generated for this target.
include src/CMakeFiles/pointer_management.dir/depend.make

# Include the progress variables for this target.
include src/CMakeFiles/pointer_management.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/pointer_management.dir/flags.make

src/CMakeFiles/pointer_management.dir/pointer_management.cpp.o: src/CMakeFiles/pointer_management.dir/flags.make
src/CMakeFiles/pointer_management.dir/pointer_management.cpp.o: src/pointer_management.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/lecoucl/Projects/EchOpen/lab_tool/RAW_tool/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/CMakeFiles/pointer_management.dir/pointer_management.cpp.o"
	cd /home/lecoucl/Projects/EchOpen/lab_tool/RAW_tool/src && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/pointer_management.dir/pointer_management.cpp.o -c /home/lecoucl/Projects/EchOpen/lab_tool/RAW_tool/src/pointer_management.cpp

src/CMakeFiles/pointer_management.dir/pointer_management.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/pointer_management.dir/pointer_management.cpp.i"
	cd /home/lecoucl/Projects/EchOpen/lab_tool/RAW_tool/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/lecoucl/Projects/EchOpen/lab_tool/RAW_tool/src/pointer_management.cpp > CMakeFiles/pointer_management.dir/pointer_management.cpp.i

src/CMakeFiles/pointer_management.dir/pointer_management.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/pointer_management.dir/pointer_management.cpp.s"
	cd /home/lecoucl/Projects/EchOpen/lab_tool/RAW_tool/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/lecoucl/Projects/EchOpen/lab_tool/RAW_tool/src/pointer_management.cpp -o CMakeFiles/pointer_management.dir/pointer_management.cpp.s

src/CMakeFiles/pointer_management.dir/pointer_management.cpp.o.requires:

.PHONY : src/CMakeFiles/pointer_management.dir/pointer_management.cpp.o.requires

src/CMakeFiles/pointer_management.dir/pointer_management.cpp.o.provides: src/CMakeFiles/pointer_management.dir/pointer_management.cpp.o.requires
	$(MAKE) -f src/CMakeFiles/pointer_management.dir/build.make src/CMakeFiles/pointer_management.dir/pointer_management.cpp.o.provides.build
.PHONY : src/CMakeFiles/pointer_management.dir/pointer_management.cpp.o.provides

src/CMakeFiles/pointer_management.dir/pointer_management.cpp.o.provides.build: src/CMakeFiles/pointer_management.dir/pointer_management.cpp.o


# Object files for target pointer_management
pointer_management_OBJECTS = \
"CMakeFiles/pointer_management.dir/pointer_management.cpp.o"

# External object files for target pointer_management
pointer_management_EXTERNAL_OBJECTS =

src/libpointer_management.a: src/CMakeFiles/pointer_management.dir/pointer_management.cpp.o
src/libpointer_management.a: src/CMakeFiles/pointer_management.dir/build.make
src/libpointer_management.a: src/CMakeFiles/pointer_management.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/lecoucl/Projects/EchOpen/lab_tool/RAW_tool/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libpointer_management.a"
	cd /home/lecoucl/Projects/EchOpen/lab_tool/RAW_tool/src && $(CMAKE_COMMAND) -P CMakeFiles/pointer_management.dir/cmake_clean_target.cmake
	cd /home/lecoucl/Projects/EchOpen/lab_tool/RAW_tool/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/pointer_management.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/pointer_management.dir/build: src/libpointer_management.a

.PHONY : src/CMakeFiles/pointer_management.dir/build

src/CMakeFiles/pointer_management.dir/requires: src/CMakeFiles/pointer_management.dir/pointer_management.cpp.o.requires

.PHONY : src/CMakeFiles/pointer_management.dir/requires

src/CMakeFiles/pointer_management.dir/clean:
	cd /home/lecoucl/Projects/EchOpen/lab_tool/RAW_tool/src && $(CMAKE_COMMAND) -P CMakeFiles/pointer_management.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/pointer_management.dir/clean

src/CMakeFiles/pointer_management.dir/depend:
	cd /home/lecoucl/Projects/EchOpen/lab_tool/RAW_tool && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/lecoucl/Projects/EchOpen/lab_tool/RAW_tool /home/lecoucl/Projects/EchOpen/lab_tool/RAW_tool/src /home/lecoucl/Projects/EchOpen/lab_tool/RAW_tool /home/lecoucl/Projects/EchOpen/lab_tool/RAW_tool/src /home/lecoucl/Projects/EchOpen/lab_tool/RAW_tool/src/CMakeFiles/pointer_management.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/pointer_management.dir/depend

