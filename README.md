Repository combines the following components:
  - A library which provides a tool to fit data (RelationHolder and ParameterKeeper).
  - Implementation of physics-methods (K-matrix, Production process, Unitarisation, Deck projection)
  - Program to fit the compass PW-data. It works based on configuration file.

Requirements:
  - `cmake` >= 3.5 (likely it works even with older version but has not been checked)
  - `gcc` >= 4.7 to support c++11 (mainly lambda functions)
  - `boost` >= 1.53 (ublas is used for matrix calculations)
  - `libconfig` 
  - `ROOT`

Installation:

1. Install required packages. 
One can check version with the following commands (centos 7).
  * `yum list installed | grep boost`
  * `cmake --version`
  * `gcc --version`
  * `root-config --version`


2. Install wavefitter
  * git clone https://github.com/mmikhasenko/wavesfitter.git
  * cd wavesfitter
  * mkdir build
  * cmake ..
  * make


## Details on Ubuntu 18 installation:
```bash
# 
sudo apt-get install cmake
sudo apt-get install libconfig-dev libconfig++-dev
# 
# install root with mathmore using sources
# 
sudo apt-get install git dpkg-dev cmake g++ gcc binutils libx11-dev libxpm-dev libxft-dev libxext-dev # root dependencies
# 
cd Documents 
git clone https://github.com/root-project/root.git root-source
mkdir root-build
cd root-build
cmake ../root-source -Dmathmore=ON
cmake --build . -- -j3
source bin/thisroot.sh
# 
#  build
git clone https://github.com/mmikhasenko/wavesfitter.git
cd wavesfitter
mkdir build
cmake .. -DCMAKE_CXX_STANDARD=14 
make wavesfitter
```