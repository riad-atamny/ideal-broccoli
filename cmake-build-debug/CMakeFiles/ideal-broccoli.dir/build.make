# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.15

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
CMAKE_COMMAND = "C:\Program Files\JetBrains\CLion 2019.2.5\bin\cmake\win\bin\cmake.exe"

# The command to remove a file.
RM = "C:\Program Files\JetBrains\CLion 2019.2.5\bin\cmake\win\bin\cmake.exe" -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = C:\Users\ASUS\CLionProjects\ideal-broccoli

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = C:\Users\ASUS\CLionProjects\ideal-broccoli\cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/ideal-broccoli.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/ideal-broccoli.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/ideal-broccoli.dir/flags.make

CMakeFiles/ideal-broccoli.dir/main.cpp.obj: CMakeFiles/ideal-broccoli.dir/flags.make
CMakeFiles/ideal-broccoli.dir/main.cpp.obj: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\ASUS\CLionProjects\ideal-broccoli\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/ideal-broccoli.dir/main.cpp.obj"
	C:\TDM-GCC-64new\bin\g++.exe  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\ideal-broccoli.dir\main.cpp.obj -c C:\Users\ASUS\CLionProjects\ideal-broccoli\main.cpp

CMakeFiles/ideal-broccoli.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ideal-broccoli.dir/main.cpp.i"
	C:\TDM-GCC-64new\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\Users\ASUS\CLionProjects\ideal-broccoli\main.cpp > CMakeFiles\ideal-broccoli.dir\main.cpp.i

CMakeFiles/ideal-broccoli.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ideal-broccoli.dir/main.cpp.s"
	C:\TDM-GCC-64new\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S C:\Users\ASUS\CLionProjects\ideal-broccoli\main.cpp -o CMakeFiles\ideal-broccoli.dir\main.cpp.s

# Object files for target ideal-broccoli
ideal__broccoli_OBJECTS = \
"CMakeFiles/ideal-broccoli.dir/main.cpp.obj"

# External object files for target ideal-broccoli
ideal__broccoli_EXTERNAL_OBJECTS =

ideal-broccoli.exe: CMakeFiles/ideal-broccoli.dir/main.cpp.obj
ideal-broccoli.exe: CMakeFiles/ideal-broccoli.dir/build.make
ideal-broccoli.exe: CMakeFiles/ideal-broccoli.dir/linklibs.rsp
ideal-broccoli.exe: CMakeFiles/ideal-broccoli.dir/objects1.rsp
ideal-broccoli.exe: CMakeFiles/ideal-broccoli.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=C:\Users\ASUS\CLionProjects\ideal-broccoli\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ideal-broccoli.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\ideal-broccoli.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/ideal-broccoli.dir/build: ideal-broccoli.exe

.PHONY : CMakeFiles/ideal-broccoli.dir/build

CMakeFiles/ideal-broccoli.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\ideal-broccoli.dir\cmake_clean.cmake
.PHONY : CMakeFiles/ideal-broccoli.dir/clean

CMakeFiles/ideal-broccoli.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" C:\Users\ASUS\CLionProjects\ideal-broccoli C:\Users\ASUS\CLionProjects\ideal-broccoli C:\Users\ASUS\CLionProjects\ideal-broccoli\cmake-build-debug C:\Users\ASUS\CLionProjects\ideal-broccoli\cmake-build-debug C:\Users\ASUS\CLionProjects\ideal-broccoli\cmake-build-debug\CMakeFiles\ideal-broccoli.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/ideal-broccoli.dir/depend

