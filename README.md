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
[] `cmake --version`
[] `yum list installed | grep boost`
[] `gcc --version`
[] `root-config --version`

2. Install wavefitter
[] git clone https://Misha_Mikhasenko@bitbucket.org/Misha_Mikhasenko/wavesfitter.git wavesfitter
[] cd wavesfitter
[] mkdir build
[] cmake ..
[] make