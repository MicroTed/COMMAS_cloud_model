#!/usr/bin/env python

# a front end to run the sips converter on a set of pdf files to make png (appends .png)
# usage:
# sipspy.pl [file names]

import os
import sys

#for arg in sys.argv:   
#       print arg

#name = sys.argv[1]
#start = int(sys.argv[2])
#stop = int(sys.argv[3])

#print 'start, stop = ' + start + ', ' + stop


for file in sys.argv:
    #num = ".%.3d." % i
    #file = name + num + 'out'
    outfile = file+'.png'
    cmd = 'sips -s format png ' + file + ' --out ' + outfile
    print (cmd)
    os.system(cmd)

