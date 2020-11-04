# RatRace

A method for indirect quad mesh generation.

# Prerequisites

This software requires third party libraries:
- OpenMesh 8.1 https://www.openmesh.org
- Glog https://github.com/google/glog

These libraries must be installed first using CMake.

Additionally, the source code for Blossom V is required: https://pub.ist.ac.at/~vnk/software.html
The code must be stored in "./src/Blossom5/".

# Setup
- This software uses CUDA. Set the according architecture in CMakeLists.txt, Line 30.
A list of architectures can be found here: https://en.wikipedia.org/wiki/CUDA

- This software uses the C++17 standard

## Windows
- create a folder "build" in the root directory of this project.
- open cmake-gui and select the root folder of this project as source and the newly created build folder as build directory.
- if you changed the install directory of the third party libraries, click "Add Entry" and enter "CMAKE_INSTALL_PREFIX". Enter the path to the installation of the third party libraries.
- Click "Configure", then "Generate".
- open "build/RatRace.sln" with Visual Studio and set "RatRace" as Startup Project.

## Linux (not tested)
- create a folder "build" in the root directory of this project.
- "cd" your way into this folder in a terminal and call "cmake ..".
- if you have changed the install directory of the third party libraries, call "cmake .. -DCMAKE_INSTALL_PREFIX=your_install_directory" instead.
- to build the project call make.

# Run
An example is provided in the main function.