#!/usr/bin/env python

# a front end to run bin.ion.pl on a whole ensemble of .out files (must give range, e.g. 1 30)
# usage:
# binpy.pl [runname] [first member number] [last member number]

import os
import sys

#for arg in sys.argv:   
#       print arg

name = sys.argv[1]
start = int(sys.argv[2])
stop = int(sys.argv[3])

#print 'start, stop = ' + start + ', ' + stop


for i in range(start, stop+1):
    num = ".%.3d." % i
    file = name + num + 'out'
    outfile = name + num + 'bin.txt'
    cmd = 'bin.ion.pl ' + file + ' >& ' + outfile
    print cmd
    os.system(cmd)

