# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
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
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/adam/OpenCCO

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/adam/OpenCCO/build

# Include any dependencies generated for this target.
include tools/CMakeFiles/graph2render.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include tools/CMakeFiles/graph2render.dir/compiler_depend.make

# Include the progress variables for this target.
include tools/CMakeFiles/graph2render.dir/progress.make

# Include the compile flags for this target's objects.
include tools/CMakeFiles/graph2render.dir/flags.make

tools/CMakeFiles/graph2render.dir/graph2render.cpp.o: tools/CMakeFiles/graph2render.dir/flags.make
tools/CMakeFiles/graph2render.dir/graph2render.cpp.o: ../tools/graph2render.cpp
tools/CMakeFiles/graph2render.dir/graph2render.cpp.o: tools/CMakeFiles/graph2render.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/adam/OpenCCO/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object tools/CMakeFiles/graph2render.dir/graph2render.cpp.o"
	cd /home/adam/OpenCCO/build/tools && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT tools/CMakeFiles/graph2render.dir/graph2render.cpp.o -MF CMakeFiles/graph2render.dir/graph2render.cpp.o.d -o CMakeFiles/graph2render.dir/graph2render.cpp.o -c /home/adam/OpenCCO/tools/graph2render.cpp

tools/CMakeFiles/graph2render.dir/graph2render.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/graph2render.dir/graph2render.cpp.i"
	cd /home/adam/OpenCCO/build/tools && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/adam/OpenCCO/tools/graph2render.cpp > CMakeFiles/graph2render.dir/graph2render.cpp.i

tools/CMakeFiles/graph2render.dir/graph2render.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/graph2render.dir/graph2render.cpp.s"
	cd /home/adam/OpenCCO/build/tools && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/adam/OpenCCO/tools/graph2render.cpp -o CMakeFiles/graph2render.dir/graph2render.cpp.s

tools/CMakeFiles/graph2render.dir/TreeImageRenderer.cpp.o: tools/CMakeFiles/graph2render.dir/flags.make
tools/CMakeFiles/graph2render.dir/TreeImageRenderer.cpp.o: ../tools/TreeImageRenderer.cpp
tools/CMakeFiles/graph2render.dir/TreeImageRenderer.cpp.o: tools/CMakeFiles/graph2render.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/adam/OpenCCO/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object tools/CMakeFiles/graph2render.dir/TreeImageRenderer.cpp.o"
	cd /home/adam/OpenCCO/build/tools && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT tools/CMakeFiles/graph2render.dir/TreeImageRenderer.cpp.o -MF CMakeFiles/graph2render.dir/TreeImageRenderer.cpp.o.d -o CMakeFiles/graph2render.dir/TreeImageRenderer.cpp.o -c /home/adam/OpenCCO/tools/TreeImageRenderer.cpp

tools/CMakeFiles/graph2render.dir/TreeImageRenderer.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/graph2render.dir/TreeImageRenderer.cpp.i"
	cd /home/adam/OpenCCO/build/tools && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/adam/OpenCCO/tools/TreeImageRenderer.cpp > CMakeFiles/graph2render.dir/TreeImageRenderer.cpp.i

tools/CMakeFiles/graph2render.dir/TreeImageRenderer.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/graph2render.dir/TreeImageRenderer.cpp.s"
	cd /home/adam/OpenCCO/build/tools && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/adam/OpenCCO/tools/TreeImageRenderer.cpp -o CMakeFiles/graph2render.dir/TreeImageRenderer.cpp.s

# Object files for target graph2render
graph2render_OBJECTS = \
"CMakeFiles/graph2render.dir/graph2render.cpp.o" \
"CMakeFiles/graph2render.dir/TreeImageRenderer.cpp.o"

# External object files for target graph2render
graph2render_EXTERNAL_OBJECTS =

tools/graph2render: tools/CMakeFiles/graph2render.dir/graph2render.cpp.o
tools/graph2render: tools/CMakeFiles/graph2render.dir/TreeImageRenderer.cpp.o
tools/graph2render: tools/CMakeFiles/graph2render.dir/build.make
tools/graph2render: /usr/local/lib/libDGtal.a
tools/graph2render: /usr/lib/x86_64-linux-gnu/libz.so
tools/graph2render: /usr/lib/x86_64-linux-gnu/libgmpxx.so
tools/graph2render: /usr/lib/x86_64-linux-gnu/libgmp.so
tools/graph2render: /usr/local/lib/libITKLabelMap-5.3.a
tools/graph2render: /usr/local/lib/libITKFastMarching-5.3.a
tools/graph2render: /usr/local/lib/libITKConvolution-5.3.a
tools/graph2render: /usr/local/lib/libITKPolynomials-5.3.a
tools/graph2render: /usr/local/lib/libITKBiasCorrection-5.3.a
tools/graph2render: /usr/local/lib/libITKColormap-5.3.a
tools/graph2render: /usr/local/lib/libITKDICOMParser-5.3.a
tools/graph2render: /usr/local/lib/libITKDeformableMesh-5.3.a
tools/graph2render: /usr/local/lib/libITKDenoising-5.3.a
tools/graph2render: /usr/local/lib/libITKDiffusionTensorImage-5.3.a
tools/graph2render: /usr/local/lib/libITKPDEDeformableRegistration-5.3.a
tools/graph2render: /usr/local/lib/libITKIOBioRad-5.3.a
tools/graph2render: /usr/local/lib/libITKIOBruker-5.3.a
tools/graph2render: /usr/local/lib/libITKIOCSV-5.3.a
tools/graph2render: /usr/local/lib/libITKIOGE-5.3.a
tools/graph2render: /usr/local/lib/libITKIOHDF5-5.3.a
tools/graph2render: /usr/local/lib/libITKIOJPEG2000-5.3.a
tools/graph2render: /usr/local/lib/libitkopenjpeg-5.3.a
tools/graph2render: /usr/local/lib/libITKIOLSM-5.3.a
tools/graph2render: /usr/local/lib/libITKIOMINC-5.3.a
tools/graph2render: /usr/local/lib/libitkminc2-5.3.a
tools/graph2render: /usr/local/lib/libITKIOMRC-5.3.a
tools/graph2render: /usr/local/lib/libITKIOSiemens-5.3.a
tools/graph2render: /usr/local/lib/libITKIOIPL-5.3.a
tools/graph2render: /usr/local/lib/libITKIOSpatialObjects-5.3.a
tools/graph2render: /usr/local/lib/libITKIOXML-5.3.a
tools/graph2render: /usr/local/lib/libITKIOStimulate-5.3.a
tools/graph2render: /usr/local/lib/libITKIOTransformHDF5-5.3.a
tools/graph2render: /usr/local/lib/libitkhdf5_cpp-static-5.3.a
tools/graph2render: /usr/local/lib/libitkhdf5_hl-static-5.3.a
tools/graph2render: /usr/local/lib/libitkhdf5-static-5.3.a
tools/graph2render: /usr/local/lib/libITKIOTransformInsightLegacy-5.3.a
tools/graph2render: /usr/local/lib/libITKIOTransformMatlab-5.3.a
tools/graph2render: /usr/local/lib/libITKIOTransformBase-5.3.a
tools/graph2render: /usr/local/lib/libITKTransformFactory-5.3.a
tools/graph2render: /usr/local/lib/libITKKLMRegionGrowing-5.3.a
tools/graph2render: /usr/local/lib/libITKMarkovRandomFieldsClassifiers-5.3.a
tools/graph2render: /usr/local/lib/libITKQuadEdgeMeshFiltering-5.3.a
tools/graph2render: /usr/local/lib/libITKRegionGrowing-5.3.a
tools/graph2render: /usr/local/lib/libITKRegistrationMethodsv4-5.3.a
tools/graph2render: /usr/local/lib/libITKImageFeature-5.3.a
tools/graph2render: /usr/local/lib/libITKOptimizersv4-5.3.a
tools/graph2render: /usr/local/lib/libITKOptimizers-5.3.a
tools/graph2render: /usr/local/lib/libitklbfgs-5.3.a
tools/graph2render: /usr/local/lib/libITKTestKernel-5.3.a
tools/graph2render: /usr/local/lib/libITKFFT-5.3.a
tools/graph2render: /usr/local/lib/libITKIOBMP-5.3.a
tools/graph2render: /usr/local/lib/libITKIOGDCM-5.3.a
tools/graph2render: /usr/local/lib/libitkgdcmMSFF-5.3.a
tools/graph2render: /usr/local/lib/libitkgdcmDICT-5.3.a
tools/graph2render: /usr/local/lib/libitkgdcmIOD-5.3.a
tools/graph2render: /usr/local/lib/libitkgdcmDSED-5.3.a
tools/graph2render: /usr/local/lib/libitkgdcmCommon-5.3.a
tools/graph2render: /usr/local/lib/libitkgdcmjpeg8-5.3.a
tools/graph2render: /usr/local/lib/libitkgdcmjpeg12-5.3.a
tools/graph2render: /usr/local/lib/libitkgdcmjpeg16-5.3.a
tools/graph2render: /usr/local/lib/libitkgdcmopenjp2-5.3.a
tools/graph2render: /usr/local/lib/libitkgdcmcharls-5.3.a
tools/graph2render: /usr/local/lib/libitkgdcmuuid-5.3.a
tools/graph2render: /usr/local/lib/libITKIOGIPL-5.3.a
tools/graph2render: /usr/local/lib/libITKIOJPEG-5.3.a
tools/graph2render: /usr/local/lib/libITKIOTIFF-5.3.a
tools/graph2render: /usr/local/lib/libitktiff-5.3.a
tools/graph2render: /usr/local/lib/libitkjpeg-5.3.a
tools/graph2render: /usr/local/lib/libITKIOMeshBYU-5.3.a
tools/graph2render: /usr/local/lib/libITKIOMeshFreeSurfer-5.3.a
tools/graph2render: /usr/local/lib/libITKIOMeshGifti-5.3.a
tools/graph2render: /usr/local/lib/libITKgiftiio-5.3.a
tools/graph2render: /usr/local/lib/libITKEXPAT-5.3.a
tools/graph2render: /usr/local/lib/libITKIOMeshOBJ-5.3.a
tools/graph2render: /usr/local/lib/libITKIOMeshOFF-5.3.a
tools/graph2render: /usr/local/lib/libITKIOMeshVTK-5.3.a
tools/graph2render: /usr/local/lib/libITKIOMeshBase-5.3.a
tools/graph2render: /usr/local/lib/libITKQuadEdgeMesh-5.3.a
tools/graph2render: /usr/local/lib/libITKIOMeta-5.3.a
tools/graph2render: /usr/local/lib/libITKMetaIO-5.3.a
tools/graph2render: /usr/local/lib/libITKIONIFTI-5.3.a
tools/graph2render: /usr/local/lib/libITKniftiio-5.3.a
tools/graph2render: /usr/local/lib/libITKznz-5.3.a
tools/graph2render: /usr/lib/x86_64-linux-gnu/libm.so
tools/graph2render: /usr/local/lib/libITKIONRRD-5.3.a
tools/graph2render: /usr/local/lib/libITKNrrdIO-5.3.a
tools/graph2render: /usr/local/lib/libITKIOPNG-5.3.a
tools/graph2render: /usr/local/lib/libitkpng-5.3.a
tools/graph2render: /usr/local/lib/libitkzlib-5.3.a
tools/graph2render: /usr/local/lib/libITKIOVTK-5.3.a
tools/graph2render: /usr/local/lib/libITKVTK-5.3.a
tools/graph2render: /usr/local/lib/libITKVideoIO-5.3.a
tools/graph2render: /usr/local/lib/libITKIOImageBase-5.3.a
tools/graph2render: /usr/local/lib/libITKVideoCore-5.3.a
tools/graph2render: /usr/local/lib/libITKWatersheds-5.3.a
tools/graph2render: /usr/local/lib/libITKMathematicalMorphology-5.3.a
tools/graph2render: /usr/local/lib/libITKStatistics-5.3.a
tools/graph2render: /usr/local/lib/libitkNetlibSlatec-5.3.a
tools/graph2render: /usr/local/lib/libITKSpatialObjects-5.3.a
tools/graph2render: /usr/local/lib/libITKMesh-5.3.a
tools/graph2render: /usr/local/lib/libITKTransform-5.3.a
tools/graph2render: /usr/local/lib/libITKPath-5.3.a
tools/graph2render: /usr/local/lib/libITKCommon-5.3.a
tools/graph2render: /usr/local/lib/libitkdouble-conversion-5.3.a
tools/graph2render: /usr/local/lib/libitksys-5.3.a
tools/graph2render: /usr/local/lib/libITKVNLInstantiation-5.3.a
tools/graph2render: /usr/local/lib/libitkvnl_algo-5.3.a
tools/graph2render: /usr/local/lib/libitkvnl-5.3.a
tools/graph2render: /usr/local/lib/libitkv3p_netlib-5.3.a
tools/graph2render: /usr/local/lib/libitkvcl-5.3.a
tools/graph2render: /usr/local/lib/libITKSmoothing-5.3.a
tools/graph2render: /usr/lib/x86_64-linux-gnu/libcairo.so
tools/graph2render: /usr/lib/x86_64-linux-gnu/libQt5OpenGL.so.5.15.3
tools/graph2render: /usr/lib/x86_64-linux-gnu/libQt5Widgets.so.5.15.3
tools/graph2render: /usr/lib/x86_64-linux-gnu/libQt5Gui.so.5.15.3
tools/graph2render: /usr/lib/x86_64-linux-gnu/libQt5Xml.so.5.15.3
tools/graph2render: /usr/lib/x86_64-linux-gnu/libQt5Core.so.5.15.3
tools/graph2render: /usr/local/lib/libQGLViewer-qt5.so
tools/graph2render: /usr/lib/x86_64-linux-gnu/libOpenGL.so
tools/graph2render: /usr/lib/x86_64-linux-gnu/libGLX.so
tools/graph2render: /usr/lib/x86_64-linux-gnu/libGLU.so
tools/graph2render: tools/CMakeFiles/graph2render.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/adam/OpenCCO/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable graph2render"
	cd /home/adam/OpenCCO/build/tools && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/graph2render.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tools/CMakeFiles/graph2render.dir/build: tools/graph2render
.PHONY : tools/CMakeFiles/graph2render.dir/build

tools/CMakeFiles/graph2render.dir/clean:
	cd /home/adam/OpenCCO/build/tools && $(CMAKE_COMMAND) -P CMakeFiles/graph2render.dir/cmake_clean.cmake
.PHONY : tools/CMakeFiles/graph2render.dir/clean

tools/CMakeFiles/graph2render.dir/depend:
	cd /home/adam/OpenCCO/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/adam/OpenCCO /home/adam/OpenCCO/tools /home/adam/OpenCCO/build /home/adam/OpenCCO/build/tools /home/adam/OpenCCO/build/tools/CMakeFiles/graph2render.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tools/CMakeFiles/graph2render.dir/depend

