# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.15

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
CMAKE_COMMAND = /home/sentinela/clion-2019.3.3/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /home/sentinela/clion-2019.3.3/bin/cmake/linux/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/sentinela/Documentos/4ano/CG/TP/fase2/engine

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/sentinela/Documentos/4ano/CG/TP/fase2/engine/cmake-build-debug

# Include any dependencies generated for this target.
include tinyxml2/CMakeFiles/tinyxml2.dir/depend.make

# Include the progress variables for this target.
include tinyxml2/CMakeFiles/tinyxml2.dir/progress.make

# Include the compile flags for this target's objects.
include tinyxml2/CMakeFiles/tinyxml2.dir/flags.make

tinyxml2/CMakeFiles/tinyxml2.dir/tinyxml2.cpp.o: tinyxml2/CMakeFiles/tinyxml2.dir/flags.make
tinyxml2/CMakeFiles/tinyxml2.dir/tinyxml2.cpp.o: ../tinyxml2/tinyxml2.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/sentinela/Documentos/4ano/CG/TP/fase2/engine/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object tinyxml2/CMakeFiles/tinyxml2.dir/tinyxml2.cpp.o"
	cd /home/sentinela/Documentos/4ano/CG/TP/fase2/engine/cmake-build-debug/tinyxml2 && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/tinyxml2.dir/tinyxml2.cpp.o -c /home/sentinela/Documentos/4ano/CG/TP/fase2/engine/tinyxml2/tinyxml2.cpp

tinyxml2/CMakeFiles/tinyxml2.dir/tinyxml2.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/tinyxml2.dir/tinyxml2.cpp.i"
	cd /home/sentinela/Documentos/4ano/CG/TP/fase2/engine/cmake-build-debug/tinyxml2 && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/sentinela/Documentos/4ano/CG/TP/fase2/engine/tinyxml2/tinyxml2.cpp > CMakeFiles/tinyxml2.dir/tinyxml2.cpp.i

tinyxml2/CMakeFiles/tinyxml2.dir/tinyxml2.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/tinyxml2.dir/tinyxml2.cpp.s"
	cd /home/sentinela/Documentos/4ano/CG/TP/fase2/engine/cmake-build-debug/tinyxml2 && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/sentinela/Documentos/4ano/CG/TP/fase2/engine/tinyxml2/tinyxml2.cpp -o CMakeFiles/tinyxml2.dir/tinyxml2.cpp.s

# Object files for target tinyxml2
tinyxml2_OBJECTS = \
"CMakeFiles/tinyxml2.dir/tinyxml2.cpp.o"

# External object files for target tinyxml2
tinyxml2_EXTERNAL_OBJECTS =

tinyxml2/libtinyxml2.a: tinyxml2/CMakeFiles/tinyxml2.dir/tinyxml2.cpp.o
tinyxml2/libtinyxml2.a: tinyxml2/CMakeFiles/tinyxml2.dir/build.make
tinyxml2/libtinyxml2.a: tinyxml2/CMakeFiles/tinyxml2.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/sentinela/Documentos/4ano/CG/TP/fase2/engine/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libtinyxml2.a"
	cd /home/sentinela/Documentos/4ano/CG/TP/fase2/engine/cmake-build-debug/tinyxml2 && $(CMAKE_COMMAND) -P CMakeFiles/tinyxml2.dir/cmake_clean_target.cmake
	cd /home/sentinela/Documentos/4ano/CG/TP/fase2/engine/cmake-build-debug/tinyxml2 && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/tinyxml2.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tinyxml2/CMakeFiles/tinyxml2.dir/build: tinyxml2/libtinyxml2.a

.PHONY : tinyxml2/CMakeFiles/tinyxml2.dir/build

tinyxml2/CMakeFiles/tinyxml2.dir/clean:
	cd /home/sentinela/Documentos/4ano/CG/TP/fase2/engine/cmake-build-debug/tinyxml2 && $(CMAKE_COMMAND) -P CMakeFiles/tinyxml2.dir/cmake_clean.cmake
.PHONY : tinyxml2/CMakeFiles/tinyxml2.dir/clean

tinyxml2/CMakeFiles/tinyxml2.dir/depend:
	cd /home/sentinela/Documentos/4ano/CG/TP/fase2/engine/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/sentinela/Documentos/4ano/CG/TP/fase2/engine /home/sentinela/Documentos/4ano/CG/TP/fase2/engine/tinyxml2 /home/sentinela/Documentos/4ano/CG/TP/fase2/engine/cmake-build-debug /home/sentinela/Documentos/4ano/CG/TP/fase2/engine/cmake-build-debug/tinyxml2 /home/sentinela/Documentos/4ano/CG/TP/fase2/engine/cmake-build-debug/tinyxml2/CMakeFiles/tinyxml2.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tinyxml2/CMakeFiles/tinyxml2.dir/depend

