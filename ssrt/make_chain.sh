#!/bin/bash

result_path="${HOME}/results/ssrt/simpleBW"

set_of_s1=(0.5 -1.0 -20.0 -100.0)
nAttempt=1000
fit_fout="/tmp/test.root"

count=0
while [ "x${set_of_s1[count]}" != "x" ] 
do
    echo "---------------------------------------------------------------------"
    echo "-----------------------------COMPILETION-----------------------------"
    echo "---------------------------------------------------------------------"
    g++ -o test plot_second_sheet.cc $(root-config --libs) -I. -I$(root-config --incdir) -DNPAR_FIT=3 -DS1=${set_of_s1[count]}
    echo "---------------------------------------------------------------------"
    echo "----------------------------RUNNING FIT------------------------------"
    echo "---------------------------------------------------------------------"
    ./test FIT $nAttempt $fit_fout
    echo "---------------------------------------------------------------------"
    echo "----------------------------PLOT RESULT------------------------------"
    echo "---------------------------------------------------------------------"
    ./test PLOT 1 $fit_fout
    root -l -q -b "~/Documents/root-scripts/browse_minimization_result.C++(\"${fit_fout}\")"
    echo "---------------------------------------------------------------------"
    echo "----------------------------FIND A POLE------------------------------"
    echo "---------------------------------------------------------------------"
    ./test FIND 100 $fit_fout
    echo "---------------------------------------------------------------------"
    echo "------------------------MOVING EVERYTHING----------------------------"
    echo "---------------------------------------------------------------------"
    mv $fit_fout ${result_path}/fit_${set_of_s1[count]}.root
    mv c2.png ${result_path}/sheets_${set_of_s1[count]}.png
    mv c3.png ${result_path}/fit_data_${set_of_s1[count]}.png
    mv /tmp/find_poles.root ${result_path}/poles_${set_of_s1[count]}.root

    count=$(( $count + 1 ))
done
