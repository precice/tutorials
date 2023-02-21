# Build instruction

To build the Oscillator FMU you need [CMake](https://cmake.org/) &GreaterEqual; 3.16 and a supported [build tool](https://cmake.org/cmake/help/latest/manual/cmake-generators.7.html) e.g. Visual Studio &GreaterEqual; 2013 , Xcode or make.

Generate the build files with the following commands:

```bash
mkdir build
cd build
cmake -DFMI_TYPE=CS -DFMI_VERSION=3 ..

```

Then run `make` or your preferred build tool to create the FMU. The model will be in the parent folder.
