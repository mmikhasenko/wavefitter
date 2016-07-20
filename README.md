# README #

This README would normally document whatever steps are necessary to get your application up and running.

### What is this repository for? ###

  - A library which provides a tool to fit data (RelationHolder and ParameterKeeper).
  - Implementation of physics-methods (K-matrix, Production process, Unitarisation, Deck projection)
  - Program to fit the compass PW-data. It works based on confiburation file.

### Dependencies ###
1. cmake >= 3.5 (likely it works even with older version but has not been checke* d)
1. gcc >= 4.7 to support c++11 (mainly lambda funct* ions)
1. boost >= 1.53 (ublas is used for matrix calcu* lations)
1. libconfig
1. ROOT


### How do I get set up? ###

* Install required packages.
One can check version with the following commands (centos 7).

```
#!bash

cmake --version
yum list installed | grep boost
gcc --version
root-config --version

```

* Install wavefitter

```
#!bash

git clone https://Misha_Mikhasenko@bitbucket.org/Misha_Mikhasenko/wavesfitter.git wavesfitter
cd wavesfitter
mkdir build
cmake ..
make

```

* Configuration

* How to run tests
  One can write tests to check out several simple examples. One finds tests in kmatrix/test_*.cc

* Deployment instructions

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* Repo owner or admin
* Other community or team contact