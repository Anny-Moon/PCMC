# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.8

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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.8.1/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.8.1/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/annasinelnikova/Documents/git_projects/PCMC/PCMC_lib/test

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/annasinelnikova/Documents/git_projects/PCMC/PCMC_lib/test

# Include any dependencies generated for this target.
include CMakeFiles/pcmcTest.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/pcmcTest.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/pcmcTest.dir/flags.make

CMakeFiles/pcmcTest.dir/source/main.cpp.o: CMakeFiles/pcmcTest.dir/flags.make
CMakeFiles/pcmcTest.dir/source/main.cpp.o: source/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/annasinelnikova/Documents/git_projects/PCMC/PCMC_lib/test/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/pcmcTest.dir/source/main.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/pcmcTest.dir/source/main.cpp.o -c /Users/annasinelnikova/Documents/git_projects/PCMC/PCMC_lib/test/source/main.cpp

CMakeFiles/pcmcTest.dir/source/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/pcmcTest.dir/source/main.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/annasinelnikova/Documents/git_projects/PCMC/PCMC_lib/test/source/main.cpp > CMakeFiles/pcmcTest.dir/source/main.cpp.i

CMakeFiles/pcmcTest.dir/source/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/pcmcTest.dir/source/main.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/annasinelnikova/Documents/git_projects/PCMC/PCMC_lib/test/source/main.cpp -o CMakeFiles/pcmcTest.dir/source/main.cpp.s

CMakeFiles/pcmcTest.dir/source/main.cpp.o.requires:

.PHONY : CMakeFiles/pcmcTest.dir/source/main.cpp.o.requires

CMakeFiles/pcmcTest.dir/source/main.cpp.o.provides: CMakeFiles/pcmcTest.dir/source/main.cpp.o.requires
	$(MAKE) -f CMakeFiles/pcmcTest.dir/build.make CMakeFiles/pcmcTest.dir/source/main.cpp.o.provides.build
.PHONY : CMakeFiles/pcmcTest.dir/source/main.cpp.o.provides

CMakeFiles/pcmcTest.dir/source/main.cpp.o.provides.build: CMakeFiles/pcmcTest.dir/source/main.cpp.o


# Object files for target pcmcTest
pcmcTest_OBJECTS = \
"CMakeFiles/pcmcTest.dir/source/main.cpp.o"

# External object files for target pcmcTest
pcmcTest_EXTERNAL_OBJECTS =

pcmcTest: CMakeFiles/pcmcTest.dir/source/main.cpp.o
pcmcTest: CMakeFiles/pcmcTest.dir/build.make
pcmcTest: /usr/local/lib/libgtest.a
pcmcTest: /usr/local/lib/libgtest_main.a
pcmcTest: /usr/local/lib/libgsl.dylib
pcmcTest: /usr/local/lib/libgslcblas.dylib
pcmcTest: /usr/local/lib/libmpi.dylib
pcmcTest: CMakeFiles/pcmcTest.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/annasinelnikova/Documents/git_projects/PCMC/PCMC_lib/test/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable pcmcTest"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/pcmcTest.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/pcmcTest.dir/build: pcmcTest

.PHONY : CMakeFiles/pcmcTest.dir/build

CMakeFiles/pcmcTest.dir/requires: CMakeFiles/pcmcTest.dir/source/main.cpp.o.requires

.PHONY : CMakeFiles/pcmcTest.dir/requires

CMakeFiles/pcmcTest.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/pcmcTest.dir/cmake_clean.cmake
.PHONY : CMakeFiles/pcmcTest.dir/clean

CMakeFiles/pcmcTest.dir/depend:
	cd /Users/annasinelnikova/Documents/git_projects/PCMC/PCMC_lib/test && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/annasinelnikova/Documents/git_projects/PCMC/PCMC_lib/test /Users/annasinelnikova/Documents/git_projects/PCMC/PCMC_lib/test /Users/annasinelnikova/Documents/git_projects/PCMC/PCMC_lib/test /Users/annasinelnikova/Documents/git_projects/PCMC/PCMC_lib/test /Users/annasinelnikova/Documents/git_projects/PCMC/PCMC_lib/test/CMakeFiles/pcmcTest.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/pcmcTest.dir/depend

