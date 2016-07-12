# README #

This README would normally document whatever steps are necessary to get your application up and running.

### What is this repository for? ###

* Quick summary:
Repository combines the followin components:
  - A library which provides a tool to fit data (RelationHolder and ParameterKeeper).
  - Implementation of physics-methods (K-matrix, Production process, Unitarisation, Deck projection)
  - Program to fit the compass PW-data. It works based on confiburation file.

### How do I get set up? ###

* Summary of set up
1. Install required packages.
One can check version with the following commands (centos 7).
[] cmake --version
[] yum list installed | grep boost
[] gcc --version
[] root-config --version

2. Install wavefitter
[] git clone https://Misha_Mikhasenko@bitbucket.org/Misha_Mikhasenko/wavesfitter.git wavesfitter
[] cd wavesfitter
[] mkdir build
[] cmake ..
[] make

* Configuration
* Dependencies
  - cmake >= 3.5 (likely it works even with older version but has not been checked)
  - gcc >= 4.7 to support c++11 (mainly lambda functions)
  - boost >= 1.53 (ublas is used for matrix calculations)
  - libconfig
  - ROOT
* Database configuration
  Configuration 

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