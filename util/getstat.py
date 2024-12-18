#!/usr/bin/env python

# a little front end to process multiple .stat or .out files with different run names.

import os
import sys

#for arg in sys.argv:   
#       print arg


basecmd = sys.argv[1]  # e.g., bin.ion.pl
infilesuf = sys.argv[2] # e.g., .000.000.out
outfilesuf = sys.argv[3] # e.g., .bin.txt
runnames = sys.argv[4:] # run names

for i in runnames:
#    print i
     infile = i + infilesuf
     outfile = i + outfilesuf
     cmd = basecmd + ' ' + infile + ' > ' + outfile
     print cmd
     os.system(cmd)

