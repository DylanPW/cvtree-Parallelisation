# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.12

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

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "C:\Program Files\JetBrains\CLion 2018.2.4\bin\cmake\win\bin\cmake.exe"

# The command to remove a file.
RM = "C:\Program Files\JetBrains\CLion 2018.2.4\bin\cmake\win\bin\cmake.exe" -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = C:\Projects\cvtree

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = C:\Projects\cvtree\cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/cvtree.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/cvtree.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/cvtree.dir/flags.make

CMakeFiles/cvtree.dir/parallel.cpp.obj: CMakeFiles/cvtree.dir/flags.make
CMakeFiles/cvtree.dir/parallel.cpp.obj: ../parallel.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Projects\cvtree\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/cvtree.dir/parallel.cpp.obj"
	C:\mingw-w64\x86_64-8.1.0-posix-seh-rt_v6-rev0\mingw64\bin\g++.exe  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\cvtree.dir\parallel.cpp.obj -c C:\Projects\cvtree\parallel.cpp

CMakeFiles/cvtree.dir/parallel.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/cvtree.dir/parallel.cpp.i"
	C:\mingw-w64\x86_64-8.1.0-posix-seh-rt_v6-rev0\mingw64\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\Projects\cvtree\parallel.cpp > CMakeFiles\cvtree.dir\parallel.cpp.i

CMakeFiles/cvtree.dir/parallel.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/cvtree.dir/parallel.cpp.s"
	C:\mingw-w64\x86_64-8.1.0-posix-seh-rt_v6-rev0\mingw64\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S C:\Projects\cvtree\parallel.cpp -o CMakeFiles\cvtree.dir\parallel.cpp.s

# Object files for target cvtree
cvtree_OBJECTS = \
"CMakeFiles/cvtree.dir/parallel.cpp.obj"

# External object files for target cvtree
cvtree_EXTERNAL_OBJECTS =

cvtree.exe: CMakeFiles/cvtree.dir/parallel.cpp.obj
cvtree.exe: CMakeFiles/cvtree.dir/build.make
cvtree.exe: CMakeFiles/cvtree.dir/linklibs.rsp
cvtree.exe: CMakeFiles/cvtree.dir/objects1.rsp
cvtree.exe: CMakeFiles/cvtree.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=C:\Projects\cvtree\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable cvtree.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\cvtree.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/cvtree.dir/build: cvtree.exe

.PHONY : CMakeFiles/cvtree.dir/build

CMakeFiles/cvtree.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\cvtree.dir\cmake_clean.cmake
.PHONY : CMakeFiles/cvtree.dir/clean

CMakeFiles/cvtree.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" C:\Projects\cvtree C:\Projects\cvtree C:\Projects\cvtree\cmake-build-debug C:\Projects\cvtree\cmake-build-debug C:\Projects\cvtree\cmake-build-debug\CMakeFiles\cvtree.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/cvtree.dir/depend

