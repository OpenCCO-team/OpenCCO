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
include tests/CMakeFiles/testGeom.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include tests/CMakeFiles/testGeom.dir/compiler_depend.make

# Include the progress variables for this target.
include tests/CMakeFiles/testGeom.dir/progress.make

# Include the compile flags for this target's objects.
include tests/CMakeFiles/testGeom.dir/flags.make

tests/CMakeFiles/testGeom.dir/testGeom.cpp.o: tests/CMakeFiles/testGeom.dir/flags.make
tests/CMakeFiles/testGeom.dir/testGeom.cpp.o: ../tests/testGeom.cpp
tests/CMakeFiles/testGeom.dir/testGeom.cpp.o: tests/CMakeFiles/testGeom.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/adam/OpenCCO/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object tests/CMakeFiles/testGeom.dir/testGeom.cpp.o"
	cd /home/adam/OpenCCO/build/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT tests/CMakeFiles/testGeom.dir/testGeom.cpp.o -MF CMakeFiles/testGeom.dir/testGeom.cpp.o.d -o CMakeFiles/testGeom.dir/testGeom.cpp.o -c /home/adam/OpenCCO/tests/testGeom.cpp

tests/CMakeFiles/testGeom.dir/testGeom.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/testGeom.dir/testGeom.cpp.i"
	cd /home/adam/OpenCCO/build/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/adam/OpenCCO/tests/testGeom.cpp > CMakeFiles/testGeom.dir/testGeom.cpp.i

tests/CMakeFiles/testGeom.dir/testGeom.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/testGeom.dir/testGeom.cpp.s"
	cd /home/adam/OpenCCO/build/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/adam/OpenCCO/tests/testGeom.cpp -o CMakeFiles/testGeom.dir/testGeom.cpp.s

# Object files for target testGeom
testGeom_OBJECTS = \
"CMakeFiles/testGeom.dir/testGeom.cpp.o"

# External object files for target testGeom
testGeom_EXTERNAL_OBJECTS =

tests/testGeom: tests/CMakeFiles/testGeom.dir/testGeom.cpp.o
tests/testGeom: tests/CMakeFiles/testGeom.dir/build.make
tests/testGeom: /usr/local/lib/libceres.a
tests/testGeom: /usr/local/lib/libDGtal.a
tests/testGeom: /usr/lib/x86_64-linux-gnu/libglog.so.0.4.0
tests/testGeom: /usr/lib/x86_64-linux-gnu/libunwind.so
tests/testGeom: /usr/lib/x86_64-linux-gnu/libgflags.so.2.2.2
tests/testGeom: /usr/lib/x86_64-linux-gnu/libspqr.so
tests/testGeom: /usr/lib/x86_64-linux-gnu/libcholmod.so
tests/testGeom: /usr/lib/x86_64-linux-gnu/libamd.so
tests/testGeom: /usr/lib/x86_64-linux-gnu/libcamd.so
tests/testGeom: /usr/lib/x86_64-linux-gnu/libccolamd.so
tests/testGeom: /usr/lib/x86_64-linux-gnu/libcolamd.so
tests/testGeom: /usr/lib/x86_64-linux-gnu/libsuitesparseconfig.so
tests/testGeom: /usr/lib/x86_64-linux-gnu/libtbb.so.12.5
tests/testGeom: /usr/lib/x86_64-linux-gnu/libcxsparse.so
tests/testGeom: /usr/lib/x86_64-linux-gnu/liblapack.so
tests/testGeom: /usr/lib/x86_64-linux-gnu/libblas.so
tests/testGeom: /usr/lib/x86_64-linux-gnu/libz.so
tests/testGeom: /usr/lib/x86_64-linux-gnu/libgmpxx.so
tests/testGeom: /usr/lib/x86_64-linux-gnu/libgmp.so
tests/testGeom: /usr/local/lib/libITKLabelMap-5.3.a
tests/testGeom: /usr/local/lib/libITKFastMarching-5.3.a
tests/testGeom: /usr/local/lib/libITKConvolution-5.3.a
tests/testGeom: /usr/local/lib/libITKPolynomials-5.3.a
tests/testGeom: /usr/local/lib/libITKBiasCorrection-5.3.a
tests/testGeom: /usr/local/lib/libITKColormap-5.3.a
tests/testGeom: /usr/local/lib/libITKDICOMParser-5.3.a
tests/testGeom: /usr/local/lib/libITKDeformableMesh-5.3.a
tests/testGeom: /usr/local/lib/libITKDenoising-5.3.a
tests/testGeom: /usr/local/lib/libITKDiffusionTensorImage-5.3.a
tests/testGeom: /usr/local/lib/libITKPDEDeformableRegistration-5.3.a
tests/testGeom: /usr/local/lib/libITKIOBioRad-5.3.a
tests/testGeom: /usr/local/lib/libITKIOBruker-5.3.a
tests/testGeom: /usr/local/lib/libITKIOCSV-5.3.a
tests/testGeom: /usr/local/lib/libITKIOGE-5.3.a
tests/testGeom: /usr/local/lib/libITKIOHDF5-5.3.a
tests/testGeom: /usr/local/lib/libITKIOJPEG2000-5.3.a
tests/testGeom: /usr/local/lib/libitkopenjpeg-5.3.a
tests/testGeom: /usr/local/lib/libITKIOLSM-5.3.a
tests/testGeom: /usr/local/lib/libITKIOMINC-5.3.a
tests/testGeom: /usr/local/lib/libitkminc2-5.3.a
tests/testGeom: /usr/local/lib/libITKIOMRC-5.3.a
tests/testGeom: /usr/local/lib/libITKIOSiemens-5.3.a
tests/testGeom: /usr/local/lib/libITKIOIPL-5.3.a
tests/testGeom: /usr/local/lib/libITKIOSpatialObjects-5.3.a
tests/testGeom: /usr/local/lib/libITKIOXML-5.3.a
tests/testGeom: /usr/local/lib/libITKIOStimulate-5.3.a
tests/testGeom: /usr/local/lib/libITKIOTransformHDF5-5.3.a
tests/testGeom: /usr/local/lib/libitkhdf5_cpp-static-5.3.a
tests/testGeom: /usr/local/lib/libitkhdf5_hl-static-5.3.a
tests/testGeom: /usr/local/lib/libitkhdf5-static-5.3.a
tests/testGeom: /usr/local/lib/libITKIOTransformInsightLegacy-5.3.a
tests/testGeom: /usr/local/lib/libITKIOTransformMatlab-5.3.a
tests/testGeom: /usr/local/lib/libITKIOTransformBase-5.3.a
tests/testGeom: /usr/local/lib/libITKTransformFactory-5.3.a
tests/testGeom: /usr/local/lib/libITKKLMRegionGrowing-5.3.a
tests/testGeom: /usr/local/lib/libITKMarkovRandomFieldsClassifiers-5.3.a
tests/testGeom: /usr/local/lib/libITKQuadEdgeMeshFiltering-5.3.a
tests/testGeom: /usr/local/lib/libITKRegionGrowing-5.3.a
tests/testGeom: /usr/local/lib/libITKRegistrationMethodsv4-5.3.a
tests/testGeom: /usr/local/lib/libITKImageFeature-5.3.a
tests/testGeom: /usr/local/lib/libITKOptimizersv4-5.3.a
tests/testGeom: /usr/local/lib/libITKOptimizers-5.3.a
tests/testGeom: /usr/local/lib/libitklbfgs-5.3.a
tests/testGeom: /usr/local/lib/libITKTestKernel-5.3.a
tests/testGeom: /usr/local/lib/libITKFFT-5.3.a
tests/testGeom: /usr/local/lib/libITKIOBMP-5.3.a
tests/testGeom: /usr/local/lib/libITKIOGDCM-5.3.a
tests/testGeom: /usr/local/lib/libitkgdcmMSFF-5.3.a
tests/testGeom: /usr/local/lib/libitkgdcmDICT-5.3.a
tests/testGeom: /usr/local/lib/libitkgdcmIOD-5.3.a
tests/testGeom: /usr/local/lib/libitkgdcmDSED-5.3.a
tests/testGeom: /usr/local/lib/libitkgdcmCommon-5.3.a
tests/testGeom: /usr/local/lib/libitkgdcmjpeg8-5.3.a
tests/testGeom: /usr/local/lib/libitkgdcmjpeg12-5.3.a
tests/testGeom: /usr/local/lib/libitkgdcmjpeg16-5.3.a
tests/testGeom: /usr/local/lib/libitkgdcmopenjp2-5.3.a
tests/testGeom: /usr/local/lib/libitkgdcmcharls-5.3.a
tests/testGeom: /usr/local/lib/libitkgdcmuuid-5.3.a
tests/testGeom: /usr/local/lib/libITKIOGIPL-5.3.a
tests/testGeom: /usr/local/lib/libITKIOJPEG-5.3.a
tests/testGeom: /usr/local/lib/libITKIOTIFF-5.3.a
tests/testGeom: /usr/local/lib/libitktiff-5.3.a
tests/testGeom: /usr/local/lib/libitkjpeg-5.3.a
tests/testGeom: /usr/local/lib/libITKIOMeshBYU-5.3.a
tests/testGeom: /usr/local/lib/libITKIOMeshFreeSurfer-5.3.a
tests/testGeom: /usr/local/lib/libITKIOMeshGifti-5.3.a
tests/testGeom: /usr/local/lib/libITKgiftiio-5.3.a
tests/testGeom: /usr/local/lib/libITKEXPAT-5.3.a
tests/testGeom: /usr/local/lib/libITKIOMeshOBJ-5.3.a
tests/testGeom: /usr/local/lib/libITKIOMeshOFF-5.3.a
tests/testGeom: /usr/local/lib/libITKIOMeshVTK-5.3.a
tests/testGeom: /usr/local/lib/libITKIOMeshBase-5.3.a
tests/testGeom: /usr/local/lib/libITKQuadEdgeMesh-5.3.a
tests/testGeom: /usr/local/lib/libITKIOMeta-5.3.a
tests/testGeom: /usr/local/lib/libITKMetaIO-5.3.a
tests/testGeom: /usr/local/lib/libITKIONIFTI-5.3.a
tests/testGeom: /usr/local/lib/libITKniftiio-5.3.a
tests/testGeom: /usr/local/lib/libITKznz-5.3.a
tests/testGeom: /usr/lib/x86_64-linux-gnu/libm.so
tests/testGeom: /usr/local/lib/libITKIONRRD-5.3.a
tests/testGeom: /usr/local/lib/libITKNrrdIO-5.3.a
tests/testGeom: /usr/local/lib/libITKIOPNG-5.3.a
tests/testGeom: /usr/local/lib/libitkpng-5.3.a
tests/testGeom: /usr/local/lib/libitkzlib-5.3.a
tests/testGeom: /usr/local/lib/libITKIOVTK-5.3.a
tests/testGeom: /usr/local/lib/libITKVTK-5.3.a
tests/testGeom: /usr/local/lib/libITKVideoIO-5.3.a
tests/testGeom: /usr/local/lib/libITKIOImageBase-5.3.a
tests/testGeom: /usr/local/lib/libITKVideoCore-5.3.a
tests/testGeom: /usr/local/lib/libITKWatersheds-5.3.a
tests/testGeom: /usr/local/lib/libITKMathematicalMorphology-5.3.a
tests/testGeom: /usr/local/lib/libITKStatistics-5.3.a
tests/testGeom: /usr/local/lib/libitkNetlibSlatec-5.3.a
tests/testGeom: /usr/local/lib/libITKSpatialObjects-5.3.a
tests/testGeom: /usr/local/lib/libITKMesh-5.3.a
tests/testGeom: /usr/local/lib/libITKTransform-5.3.a
tests/testGeom: /usr/local/lib/libITKPath-5.3.a
tests/testGeom: /usr/local/lib/libITKCommon-5.3.a
tests/testGeom: /usr/local/lib/libitkdouble-conversion-5.3.a
tests/testGeom: /usr/local/lib/libitksys-5.3.a
tests/testGeom: /usr/local/lib/libITKVNLInstantiation-5.3.a
tests/testGeom: /usr/local/lib/libitkvnl_algo-5.3.a
tests/testGeom: /usr/local/lib/libitkvnl-5.3.a
tests/testGeom: /usr/local/lib/libitkv3p_netlib-5.3.a
tests/testGeom: /usr/local/lib/libitkvcl-5.3.a
tests/testGeom: /usr/local/lib/libITKSmoothing-5.3.a
tests/testGeom: /usr/lib/x86_64-linux-gnu/libcairo.so
tests/testGeom: /usr/lib/x86_64-linux-gnu/libQt5OpenGL.so.5.15.3
tests/testGeom: /usr/lib/x86_64-linux-gnu/libQt5Widgets.so.5.15.3
tests/testGeom: /usr/lib/x86_64-linux-gnu/libQt5Gui.so.5.15.3
tests/testGeom: /usr/lib/x86_64-linux-gnu/libQt5Xml.so.5.15.3
tests/testGeom: /usr/lib/x86_64-linux-gnu/libQt5Core.so.5.15.3
tests/testGeom: /usr/local/lib/libQGLViewer-qt5.so
tests/testGeom: /usr/lib/x86_64-linux-gnu/libOpenGL.so
tests/testGeom: /usr/lib/x86_64-linux-gnu/libGLX.so
tests/testGeom: /usr/lib/x86_64-linux-gnu/libGLU.so
tests/testGeom: tests/CMakeFiles/testGeom.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/adam/OpenCCO/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable testGeom"
	cd /home/adam/OpenCCO/build/tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/testGeom.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tests/CMakeFiles/testGeom.dir/build: tests/testGeom
.PHONY : tests/CMakeFiles/testGeom.dir/build

tests/CMakeFiles/testGeom.dir/clean:
	cd /home/adam/OpenCCO/build/tests && $(CMAKE_COMMAND) -P CMakeFiles/testGeom.dir/cmake_clean.cmake
.PHONY : tests/CMakeFiles/testGeom.dir/clean

tests/CMakeFiles/testGeom.dir/depend:
	cd /home/adam/OpenCCO/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/adam/OpenCCO /home/adam/OpenCCO/tests /home/adam/OpenCCO/build /home/adam/OpenCCO/build/tests /home/adam/OpenCCO/build/tests/CMakeFiles/testGeom.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/CMakeFiles/testGeom.dir/depend

