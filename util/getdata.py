#!/opt/local/anaconda/bin/python

# the Scientific Python netCDF 3 interface
# http://dirac.cnrs-orleans.fr/ScientificPython/
# from Scientific.IO.NetCDF import NetCDFFile as Dataset
# the 'classic' version of the netCDF4 python interface
# http://code.google.com/p/netcdf4-python/
#from netCDF4_classic import Dataset
import netCDF4 as netcdf
#from numpy import arange # array module from http://numpy.scipy.org
#from numpy.testing import assert_array_equal, assert_array_almost_equal
import numpy as np
import numpy.ma as ma
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys

# Usage: read any number of RUNNAME from the command line
# The script will loop over the run names and read the netcdf files to calculate
# cold pool perturbations (T and theta-v)

filename = sys.argv[1]
print 'file name base = '+filename
# open data file as read only 'r'
narg = len(sys.argv)
print "Number of arguments: ", len(sys.argv)

for iname in range(1,narg):
    print 'name = ',sys.argv[iname]
    filename = sys.argv[iname]


    #choke
    
    #filename = 'may29zh1000cn10vt.000.nc'
    # filename = 'may29zvh1000cn13vtrr'
    # filename = 'wk10zhe100mpitrajr3a.000.000.traj1.nc'
    ncfile = netcdf.Dataset(filename+'.000.nc','r') 
    
    #print ncfile.file_format
    
    #print ncfile.dimensions.keys()
    
    #print ncfile.dimensions['TIME']
    
    times = ncfile.variables['TIME']
    
    nsteps = times.shape[0]
    
    print 'nsteps  = ',nsteps
    
    #print times[:nsteps]
    
    #print ncfile.dimensions['THETA-V']
    
    thv2d = ncfile.variables['THETA-V'][0,0,:,:]
    nx = thv2d.shape[1]
    ny = thv2d.shape[0]
    print 'thv2d shape ', thv2d.shape,ny,nx
    
    th2d = ncfile.variables['TH'][0,0,:,:]
    t2d = ncfile.variables['T'][0,0,:,:]
    th0 = th2d[0,0]
    t0 = t2d[0,0]
    thv0 = thv2d[0,0]
    
    aveval0 = np.average(thv2d)
    
    print aveval0, thv2d[0,0],th2d[0,0]
    
    thvthresh = thv0-0.5
    ththresh = th0-0.5
    tthresh = t0-0.5
    
    fileout = open(filename+'.coldpool.txt', 'w')
    fileout.write('Time (min),T-ave,T-area,T-sum,Theta-V-ave, Theta-V-area, Theta-V-sum\n')
    
    for t in range(nsteps):
    #   print 'time = ',times[t]
    #   if ( times[t] == 1440 ):
    #      continue
	 thv2d = ncfile.variables['THETA-V'][t,0,:,:]
	 th2d = ncfile.variables['TH'][t,0,:,:]
	 t2d = ncfile.variables['T'][t,0,:,:]
	 amask = np.where( thv2d[:,:] > thvthresh, 1,0 ) 
    #   thv2dm = ma.masked_array(thv2d,mask=amask)
	 thv2dm = ma.masked_greater(thv2d,thvthresh) - thv0
	 
	 num1 = np.ma.sum( 1 - amask ) # cold pool area (number of points)
	 aveval1 = 0
	 if ( num1 > 0 ):
	    sum1 = np.ma.sum( thv2dm )  # integrated cold pool
	    aveval1 = np.ma.average( thv2dm )
	 else:
	    aveval1 = 0
	    sum1 = 0
	    num1 = 0
    
    
	 amasktem = np.where( t2d[:,:] > tthresh, 1,0 ) 
	 t2dm = ma.masked_greater(t2d,tthresh) - t0
	 numtem = np.ma.sum( 1 - amasktem ) # cold pool area (number of points)
	 avevaltem = 0
	 if ( num1 > 0 ):
	    sumtem = np.ma.sum( t2dm )  # integrated cold pool
	    avevaltem = np.ma.average( t2dm ) # note: have to use "ma" version to honor the mask
	 else:
	    avevaltem = 0
	    sumtem = 0
	    numtem = 0
    
    
    #   aveval = np.average(thv2d - aveval0)
	 thv2dp = thv2d - aveval0
	 aveval = np.ma.average( ma.masked_array(thv2dp,mask=amask))
	 sum = 0
	 npoints = 0
	 sumth = 0
	 npointsth = 0
	 sumall = 0
	 nall = 0
    #    for i in range(nx):
    #       for j in range(ny):
    #          sumall = sumall + thv2d[j,i]
    #          nall = nall + 1
    #          if ( thv2d[j,i] < thvthresh ):
    #             sum = sum + thv2d[j,i] - thv0
    #             npoints = npoints + 1
    #          if ( th2d[j,i] < ththresh ):
    #             sumth = sumth + th2d[j,i] - th0
    #             npointsth = npointsth + 1
	 
    #   if ( nall == 0 ):
    #      nall = 1
    #   aveval2 = 0
    #   if ( npoints > 0 ):
    #      aveval2 = ( sum)/npoints 
    #   aveval2th = 0
    #   if ( npointsth > 0 ):
    #      aveval2th = ( sumth)/npointsth
    
    #   print times[t]/60, aveval,aveval1,(sumall/nall)-thv0,aveval2,aveval2th,npoints,npointsth,num1,sum1
    #   print times[t]/60, avevaltem,numtem,sumtem,aveval1,num1,sum1
	 fileout.write( str(times[t]/60)+','+str(avevaltem)+','+str(numtem)+','+str(sumtem)+','+str(aveval1)+','+str(num1)+','+str(sum1)+'\n')


