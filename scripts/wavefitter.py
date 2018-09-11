#!/usr/bin/env python

# Misha Mikhasenko  [03.01.2018]

# parse parameters
# this script suppose to fill gaps in config file
# and create new combined config file
#
# There should be possibility to leave placeholders in the *.cfg files
# then, the placeholders are filled by arguments of the file

import argparse
import random
import os

parser = argparse.ArgumentParser(
                    description="""
Create customized config file and run wavefitter with these settings.
The files \"model.cfg\", \"fit.cgf\", \"continuation.cgf\", \"plot.cgf\" should be present in the folder.
The use <0>, <1>, and so on for the placeholders in the config files. It will be filled by the provided args
(see options below).
                                 """)
# model arguments
parser.add_argument('--model_args', action='append',
                    default=[],
                    help='list values for placeholders in the model.cfg in the right order',
                    )
# fit option and arguments
parser.add_argument('--fit', action='store_true', default=False,
                    help='Set a add \"fit\" section to config file')
parser.add_argument('--fit_args', action='append',
                    default=[],
                    help='list values for placeholders in the fit.cfg in the right order',
                    )
# continue option and arguments
parser.add_argument('--continuation', action='store_true', default=False,
                    help='Set a add \"continuation\" section to config file')
parser.add_argument('--continuation_args', action='append', default=[],
                    help='list values for placeholders in the continuation.cfg in the right order',
                    )
# plot option and arguments
parser.add_argument('--plot', action='store_true', default=False,
                    help='Set a add \"plot\" section to config file')
parser.add_argument('--plot_args', action='append', default=[],
                    help='list values for placeholders in the plot.cfg in the right order',
                    )
parser.add_argument('--version', action='version', version='%(prog)s 1.0')

results = parser.parse_args();
# print results
print 'model_args        =', results.model_args
print 'fit               =', results.fit
print 'fit_args          =', results.fit_args
print 'continuation      =', results.continuation
print 'continuation_args =', results.continuation_args
print 'plot              =', results.plot
print 'plot_args         =', results.plot_args

# files are assumed in the directory
# might be it is worth to parse them later
model_name = "model.cfg"
fit_name = "fit.cfg"
continuation_name = "continuation.cfg"
plot_name = "plot.cfg"

# list of replacements to get around {} in the config file
reps      = [["{", "|-"],
             ["}", "-|"],
             ["<", "{"],
             [">", "}"]]
reps_back = [["|-", "{"],
             ["-|", "}"]];
def make_replacements(line, args) :
                    for r in reps : line = line.replace(*r)
                    line = line.format(*args)
                    for r in reps_back: line = line.replace(*r)
                    return line

# generate a name
fout_name = "model.auto.cfg"
while True:
                    fout_name = "model.auto.cfg."+str(random.randint(10000,99999))
                    if not os.path.isfile(fout_name) :
                                        break
fout = open(fout_name, "w")
fout.writelines("###############\n#\n# The file is created automatically using wavefitter.py script\n#\n###############\n\n")
# copy model content to the file
with open(model_name) as f:
    lines = f.readlines()
    lines = [make_replacements(l,results.model_args) for l in lines]
    fout.writelines(lines)

if results.fit:
                    fout.writelines("\n\n## FIT SECTION ##\n\n")
                    with open(fit_name) as f:
                        lines = f.readlines()
                        lines = [make_replacements(l, results.fit_args) for l in lines]
                        fout.writelines(lines)

if results.continuation:
                    fout.writelines("\n\n## CONTINUATION SECTION ##\n\n")
                    with open(continuation_name) as f:
                        lines = f.readlines()
                        lines = [make_replacements(l, results.continuation_args) for l in lines]
                        fout.writelines(lines)

if results.plot:
                    fout.writelines("\n\n## PLOT SECTION ##\n\n")
                    with open(plot_name) as f:
                        lines = f.readlines()
                        lines = [make_replacements(l,results.plot_args) for l in lines]
                        fout.writelines(lines)

print "--> config file " + fout_name + " has been created!"
print "   less " + fout_name
fout.close()

# run the program
print "Now execute wavefitter with the config file"
command = "  ./wavefitter " + fout_name
print command
os.system(command)
# print subprocess.call(["wavefitter " + fout_name])
