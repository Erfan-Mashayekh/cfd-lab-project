![](FluidchenLogo.png)

Fluidchen is a CFD Solver developed for the CFD Lab taught at TUM Informatics, Chair of Scientific Computing in Computer Science.

After forking, use this `README.md` however you want: feel free to remove anything you don't need,
or add any additional details we should know to run the code.

## Working with fluidchen

You will extend this code step-by-step starting from a pure framework to a parallel CFD solver. Please follow these [instructions for work with git and submitting the assignments](docs/first-steps.md).

## Software Requirements

* VTK 7 or higher
* GCC 9 (optional)
  
Detailed information is given below.

## Installing

```shell
git clone https://gitlab.lrz.de/oguzziya/GroupX_CFDLab.git
cd GroupX_CFDLab
mkdir build && cd build
cmake ..
make
make install # optional
```

These commands will create the executable `fluidchen` and copy it to the default directory `/usr/local/bin` . If you want to install to a custom path, execute the cmake command as

```shell
cmake -DCMAKE_INSTALL_PREFIX=/path/to/directory ..
```

After `make && make install` **fluidchen** will be installed to `/path/to/directory/bin` . Note that you may need to update your `PATH` environment variable.

By default, **fluidchen** is installed in `DEBUG` mode. To obtain full performance, you can execute cmake as

```shell
cmake -DCMAKE_BUILD_TYPE=RELEASE ..
```

or

```shell
cmake -DCMAKE_CXX_FLAGS="-O3" ..
```

A good idea would be that you setup your computers as runners for [GitLab CI](https://docs.gitlab.com/ee/ci/)
(see the file `.gitlab-ci.yml` here) to check the code building automatically every time you push.

## Running

In order to run **Fluidchen**, the case file should be given as input parameter. Some default case files are located in the `example_cases` directory. If you installed **Fluidchen**, you can execute them from anywhere you want as
For Serial:

```shell
fluidchen /path/to/fluidchen/example_cases/LidDrivenCavity/LidDrivenCavity.dat
```

This will run the case file and create the output folder `/path/to/case/case_name_Output` which holds the `.vtk` files of the solution. The output folder is created in the same location as your case file. Note that this may require write permissions in the given directory.

It is essential to provide a geometry file. If input file does not contain a geometry file it will issue an error and would exit. Even for the default LidDrivenCavity example, a geometry file has been added. 

### GCC version

You can get you current version of GCC by running:

```shell
g++ -v
```

### Defining your GCC version

If you have GCC 9 or newer, you can set in the `CMakeLists.txt` file:

```cmake
set(gpp9 True)
```

If you have a version lower than 9, then you don't have to modify the `CMakeLists.txt` file.

This will affect how we are using the C++ filesystem library, which is available already in GCC 7 as an experimental feature.

### Setup of VTK and GCC 9 (Ubuntu **20.04**)

```shell
apt-get update &&
apt-get upgrade -y &&
apt-get install -y build-essential cmake libvtk7-dev libfmt-dev openmpi-bin libopenmpi-dev
```

### Setup of VTK and GCC 9 (Ubuntu **18.04**)

If you want, you can upgrade your compiler version to have access to more recent C++ features.
This is, however, optional.

```shell
apt-get update &&
apt-get install -y software-properties-common &&
add-apt-repository -y ppa:ubuntu-toolchain-r/test &&
apt-get upgrade -y &&
apt-get install -y build-essential cmake libvtk7-dev libfmt-dev gcc-9 g++-9 
apt-get install -y gcc-9 g++-9 openmpi-bin libopenmpi-dev
```

### Dependencies for macOS

In macOS, you can use default `clang` compiler. Do not install `gcc` compiler since it might cause problems with the standard library and VTK. Other dependencies can be installed by using `homebrew` package manager as

```shell
brew install cmake
brew install open-mpi
brew install vtk
```

**macOS Troubleshooting**
- In macOS, the default `g++` command is linked to `clang++` command, which means, `g++` command does not run the GCC compiler but the Clang compiler. 
- Setup of GCC compiler is experienced to be cumbersome and clashes with lots of other dependencies, therefore please do not use GCC compiler on this project.
- If CMake cannot find the correct C++ binary, you can set it by
```
export CXX=`which clang++``
export CMAKE_CXX_COMPILER=`which clang++``
```
which is going to set the corresponding environment variables to the path of Clang compiler. Please note that if you run these commands on a terminal session, they are only going to be valid on that terminal session. In order to make these changes permanent, you can add these lines to your `~/.zshrc` file.
- Although installation of VTK looks simple, sometimes it is possible that CMake cannot find some necessary libraries for VTK, most famous one being Qt5. If you face an error something like:
```
CMake Error at /usr/local/lib/cmake/vtk-9.0/VTK-vtk-module-find-packages.cmake:115 (find_package):
 By not providing "FindQt5.cmake" in CMAKE_MODULE_PATH this project has
 asked CMake to find a package configuration file provided by "Qt5", but
 CMake did not find one.

 Could not find a package configuration file provided by "Qt5" (requested
 version 5.15) with any of the following names:

   Qt5Config.cmake
   qt5-config.cmake

 Add the installation prefix of "Qt5" to CMAKE_PREFIX_PATH or set "Qt5_DIR"
 to a directory containing one of the above files.  If "Qt5" provides a
 separate development package or SDK, be sure it has been installed.

```
which means that CMake could not find Qt5. Solution is simple fortunately. First, make sure that you have Qt5 installed:
```
brew install qt5
```
Then extend `CMAKE_PREFIX_PATH`, which are the locations where CMake tries to find packages, by adding following lines to your `.zshrc` file
```
export CMAKE_PREFIX_PATH="/usr/local/opt/qt5:$CMAKE_PREFIX_PATH"
```
Please not that your installation might be in another location. The most possible another location is `/usr/local/Cellar/qt@5/5.**.*/`, which depends on the Qt5 version. 

## Using CMake

CMake is a C++ build system generator, which simplifies the building process compared e.g. to a system-specific Makefile. The CMake configuration is defined in the `CMakeList.txt` file.

In order to build your code with CMake, you can follow this (quite common) procedure:

1. Create a build directory: `mkdir build`
2. Get inside it: `cd build`
3. Configure and generate the build system: `cmake ..` (Note the two dots, this means that the `CmakeLists.txt` File is in the folder above)
4. Build your code: `make` (build the executable)

### Troubleshooting: VTK not found

You might run into a problem where the VTK library is not found. To fix this, you can try the following steps:

1. Find the installation path of your VTK library 
2. Define this path as an environment variable, as e.g. `export VTK_DIR=".../lib/cmake/vtk-8.2"`
3. Start in a clean build folder
4. Run `cmake ..` again

### Set a different GCC version

If you have multiple compiler versions installed you can set the GCC version which should be used by `cmake` like this:

```shell
export CXX=`which g++-7`
```

Make sure to use a backtick (\`) to get the `which` command executed. Afterwards, you can run `cmake ..`.
