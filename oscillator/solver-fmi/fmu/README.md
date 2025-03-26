# Build instruction

To build the Oscillator FMU you need [CMake](https://cmake.org/) &GreaterEqual; 3.16 and a supported [build tool](https://cmake.org/cmake/help/latest/manual/cmake-generators.7.html) e.g. Visual Studio &GreaterEqual; 2013 , Xcode or make.

Generate the build files with the following commands:

```bash
mkdir build
cd build
cmake -DFMI_TYPE=CS -DFMI_VERSION=3 ..
make
cp ./Oscillator.fmu ../..
```

Instead of `make`, you can also use your preferred build tool to create the FMU.

## License

The FMU model `Oscillator.fmu` used for this tutorial is based on the [Reference-FMUs](https://github.com/modelica/Reference-FMUs) from the Modelica Association Project "FMI", which are released under the [2-Clause BSD License](https://github.com/precice/tutorials/blob/master/oscillator/fmi/fmu/thirdparty/LICENSE.txt).
