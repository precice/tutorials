rm build -rf
mkdir build
cd build
cmake .. -DFEBio_SDK="~/FEBioStudio/sdk"
make
