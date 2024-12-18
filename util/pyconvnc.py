#!/usr/bin/env python
#
import sys

from optparse import OptionParser
import time as timer

import glob,os
import re

RunIt = True #  debug switch, False:  will print out ncks commands, but not run the compression

# ncks only converts netcdf3 to compressed (4-byte) netcdf4 -- nonlossy.
ncpdq_string = "ncks -7 -L 2"

# NCOMMAS fields
check_fields = ["W","TH", "QV", "WZ"]

# WRF OUTPUT
#check_fields = ["W"]


#---------------------------------------------------------------------------------------------------
# Use a regular expression to match the given prefix with files having 3 or 6 numeral digits and then ".nc"
#

def find_ncfiles(top_dir, prefix):

    files = []
    ncdf_files = os.listdir(top_dir)

    for item in ncdf_files:
        file = re.match(r'\b(?=\w)%s\b\.\d{3}.nc' % prefix, item, re.M|re.I)
        if file:
            files.append(file.group())
        file = re.match(r'\b(?=\w)%s\b\.\d{6}.nc' % prefix, item, re.M|re.I)
        if file:
            files.append(file.group())
   
    if len(files) > 0:
        return files
    else:
        return None
        
#---------------------------------------------------------------------------------------------------
# Compress a list of netcdf files....
#

def compress_ncfile(files, cmpdir, no_comp = False, stats=False):

    for file in files:
        infile  = os.path.join(root,file)
        outfile = os.path.join(cmpdir,file)
        
        if not no_comp:
            tt = timer.clock()
            cmd = "%s %s %s" % (ncpdq_string, infile, outfile)
            print "Executing....\n%s" % cmd
            
            if RunIt:  
                os.system(cmd)
                
            print round(timer.clock()-tt, 4),' secs to compress file....'
            
        if stats:
            infile_size = os.path.getsize(infile)
            outfile_size = os.path.getsize(outfile)
            print("Original file size:  %d mb / Compressed file size %d mb / Compression Ratio:  %d %s\n" \
                   % (infile_size/1000000, outfile_size/1000000, float(infile_size)/float(outfile_size)*100, "%"))
              
#---------------------------------------------------------------------------------------------------
# Main function defined to return correct sys.exit() calls
#

parser = OptionParser()
parser.add_option("-d", "--dir", dest="dir", type="string", default= None, \
                                  help="Name of netCDF file directory to compress")
parser.add_option("-p", "--prefix", dest="prefix", type="string", default= None, \
                                  help="Name of netCDF file prefixes to compress")
                                    
parser.add_option("-s","--stats", dest="stats", action="store_true", \
                              help="Compute error statistics from lossy compression")

parser.add_option("-n","--nocmp", dest="no_comp", action="store_true", \
                              help="Check differences without running compression (aka, after you created files)")

parser.add_option("-b","--byte", dest="byte", action="store_true", \
                              help="Compress down to 8 bit integers (default is 16 bit)")

(options, args) = parser.parse_args()

print """\n=====================================================================================================================

          pyconvnc:  A python script for compression of WRF/NCOMMAS netCDF ensembles using NetCDF Operators.

          This program will perform a nonlossy compression of netcdf files using the netCDF Operators (NCO) utility  
          'ncks' and uses python to walk down a directory structure to find files that match a prefix supplied 
          by the user. Only netCDF variables are compressed 
          Coordinate variables are not compressed.  The user can enable RMS statistics between the original fields 
          and the compressed fields using the '-s' command line flag.  
          
          Requirements:  One needs the netCDF Operators to be installed (which they are on many systems).  
                         A valid python installation with Jeff Whitaker's netCDF4 python interface is optimal,
                         but the program will run without netCDF4 installed without statistics capability. 
                         
          Example:  To compress files in a directory: 
             
             /Thor10TB/dwheatley/mybigrun/experiment_1 
                    
             using 16-bit integers for a series of files called wrfout and check the stats..

             Usage:  cmpdir.py -d /Thor10TB/dwheatley/mybigrun/experiment_1 -p wrfout -s 
          
          Example:  to compress the same directory using 8-bit integers and checking the stats:

             Usage:  cmpdir.py -d /Thor10TB/dwheatley/mybigrun/experiment_1 -p wrfout -s -b
          
          The script will also traverse down from a top level directory and compress all sub-directories
          that have files matching the user supplied prefix.  So large multi-experiment and ensemble directories can
          be automatically compressed.  
          
          If you are nervous about compression errors, the 16-bit compression is essentially lossless (except for variables
          that have a very large dynamic range - that may require some testing) and yet one still gains a factor
          of 2.5-3x in space.  8-bit compression will reduce the file sizes by a factor of 10 or more!
          
          Importantly, NCL and python netCDF, and most other tools can READ THE COMPRESSED FILES DIRECTLY - requiring 
          no decompression by the user prior to use.  Other tools may work this way as well.
          
          Lou Wicker, 19 Aug 2014. 
          
          =====================================================================================================================
"""

if options.dir == None:
    parser.print_help()
    print
    sys.exit(1)
else:
    if not os.path.exists(options.dir):
        print "\nError!  netCDF directory does not seem to exist?\n"
        print "Filename:  %s" % (options.dir)
        sys.exit(1)
    else:
        print "\n--->  Top level root directory to look for netCDF files:  %s" % options.dir
        
if options.prefix == None:
    print "\nError!  cmpdir requires a prefix name for compression....\n"
    parser.print_help()
    print
    sys.exit(1)    

if options.byte:
    ncpdq_string = ("%s %s" % (ncpdq_string, "-M flt_byt"))
    print "\n--->  Compressing to 8 BIT INTEGERS!!\n"

#---------------------------------------------------------------------------------------------------
# Use a recursive walk to get to subfolders      

for root, subFolders, files in os.walk(options.dir):
    
    if root.find("_cmp") == -1:   # Check for compression directories already created
        ncdf_files = find_ncfiles(root, options.prefix)   # Find netCDF files that match the prefix pattern
    else:
        continue                  # If its a compression directory, go back and walk to next directory

    if ncdf_files:                # Are there files that match the netCDF pattern?
    
        cmpdir = root+"_cmp"
        if not os.path.exists(cmpdir):   # If the compressed file directory does not exist, create it...
            os.mkdir(cmpdir)
             
        compress_ncfile(ncdf_files, cmpdir, no_comp = options.no_comp, stats=options.stats) # Compress files in list
    

        
        
        
        