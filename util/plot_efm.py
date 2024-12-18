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
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_pdf import PdfPages

# open data file as read only 'r'

filename = 'wk10zhe100mpitraj002700.snapshot.traj1.nc'
# filename = 'wk10zhe100mpitrajr3a.000.000.traj1.nc'
ncfile = netcdf.Dataset(filename,'r') 

# constants
eps0 = 8.85e-12

# get the x,y, and z locations; factor of 0.001 converts meters to km
xloc = 0.001*ncfile.variables['X'][:]
yloc = 0.001*ncfile.variables['Y'][:]
zloc = 0.001*ncfile.variables['Z'][:]
times = ncfile.variables['TIME']
nx,ntraj = zloc.shape
nsteps = times.shape[0]
print 'nx,ntraj,nsteps  = ', nx,ntraj,nsteps
print 'x1,y1,z1 = ',xloc[0,0], yloc[0,0], zloc[0,0]
print 'time shape ', times.shape
print 'time1,2,3 = ',times[0],times[1],times[2]
#print 'x0 = ',xloc[0][:]

fig = plt.figure()
ax = fig.gca(projection='3d')

# make a plot of all sounding trajectories into the same plot:
for i in range (ntraj):
  xs = xloc[:,i]
  ys = yloc[:,i]
  zs = zloc[:,i]
  
  ax.plot(xs, ys, zs)

ax.set_xlabel("X (km)")
ax.set_ylabel("Y (km)")
ax.set_zlabel("Altitude (km)")
ax.set_title("Trajectories")
ax.set_zbound(0,None)
# plt.savefig('traj3test.png',dpi=150)
plt.savefig(filename+'.trajplot.pdf',format='pdf')
plt.close()

# get pre-lightning Ez in variable ez (factor 0.001 to convert V/m to kV/m):
emag = 0.001*ncfile.variables['EMAG-PreLightning'][:]
ex = 0.001*ncfile.variables['EX-PreLightning'][:]
ey = 0.001*ncfile.variables['EY-PreLightning'][:]

ez = 0.001*ncfile.variables['EZ-PreLightning'][:]
#ez = 0.001*ncfile.variables['EZ'][:]
print 'min,max ez = ',np.min(ez),np.max(ez)

# get pre-lightning potential (convert from V to MV)
# potential = 1.e-6*ncfile.variables['POTENTIAL-PreLightning'][:]
potential = 1.e-6*ncfile.variables['POTENTIAL'][:]
print 'min,max potential = ',np.min(potential),np.max(potential)

# get air temperature
airtemp = ncfile.variables['T'][:] - 273.15

# balloon rise rate
riserate = ncfile.variables['RISERATE'][:] 

# net charge density 
chgnet = ncfile.variables['CHGNET'][:] 
chgnetave = ncfile.variables['CHGNETAVE'][:] 
print 'min,max chgnet = ',np.min(chgnet),np.max(chgnet)

# number of flashes during the time step
numflash = ncfile.variables['NUMFLASH'][:]
print 'min,max numflash = ',np.min(numflash),np.max(numflash)

# test code to plot z and Ez for the ith trajectory:
i = 1

# loop for printing out values of altitude (z) and vertical e-field (ez) of the i-th sounding
#for j in xrange(0,nsteps,20):
#  print 'z,ez = ',  zloc[j][i], ', ',ez[j][i]

# This makes an x-y plot of Ez (x-axis) with altitude (zloc, y-axis)
plt.plot(ez[:,i], zloc[:,i])
plt.axis([-150,150,0,16])
# plt.show()
plt.close()

potentialcalc = 0.0*ncfile.variables['POTENTIAL-PreLightning'][:]
potentialcalc2 = 0.0*ncfile.variables['POTENTIAL-PreLightning'][:]
potentialcalc3d = 0.0*ncfile.variables['POTENTIAL-PreLightning'][:]
scnetcalc1d = 0.0*ncfile.variables['POTENTIAL-PreLightning'][:]
scnetcalc3d = 0.0*ncfile.variables['POTENTIAL-PreLightning'][:]

for i in range (ntraj):
#for i in range (70,73):
  for j in range  (1,nsteps):
    deltaz  = zloc[j,i] - zloc[j-1,i]
    deltax  = xloc[j,i] - xloc[j-1,i]
    deltay  = yloc[j,i] - yloc[j-1,i]
    exave = 0.5*(ex[j-1,i] + ex[j,i])
    eyave = 0.5*(ey[j-1,i] + ey[j,i])
    ezave = 0.5*(ez[j-1,i] + ez[j,i])
#    ezmagave = ezave/abs(ezave)*0.5*(emag[j-1,i] + emag[j,i])
    ezmagave = np.sign(ezave)*0.5*(emag[j-1,i] + emag[j,i])
    potentialcalc[j,i] = potentialcalc[j-1,i] - deltaz*ezave
    potentialcalc2[j,i] = potentialcalc2[j-1,i] - deltaz*ez[j-1,i]
    potentialcalc3d[j,i] = potentialcalc3d[j-1,i] - deltax*exave - deltay*eyave - deltaz*ezave


# calculate scnet from 1D Gauss Law
for i in range (ntraj):
#for i in range (70,73):
  for j in range  (3,nsteps-2):
    deltaz  = zloc[j,i] - zloc[j-1,i]
    deltax  = xloc[j,i] - xloc[j-1,i]
    deltay  = yloc[j,i] - yloc[j-1,i]

    deltaz2  = zloc[j+1,i] - zloc[j,i]
    deltax2  = xloc[j+1,i] - xloc[j,i]
    deltay2  = yloc[j+1,i] - yloc[j,i]

    dphi   = potentialcalc3d[j,i] - potentialcalc3d[j-1,i]
    dphi2  = potentialcalc3d[j+1,i] - potentialcalc3d[j,i]
    
    if (deltax == 0) :
       dphidx = 0
    else:
       dphidx = dphi/deltax
    
    if (deltay == 0) :
       dphidy = 0
    else:
       dphidy = dphi/deltay

    if (deltaz == 0) :
       dphidz = 0
    else:
       dphidz = dphi/deltaz



    if (deltax2 == 0) :
       dphidx2 = 0
    else:
       dphidx2 = dphi2/deltax2
    
    if (deltay2 == 0) :
       dphidy2 = 0
    else:
       dphidy2 = dphi2/deltay2

    if (deltaz2 == 0) :
       dphidz2 = 0
    else:
       dphidz2 = dphi2/deltaz2

#    dphidx2 = dphi2/deltax2
#    dphidy2 = dphi2/deltay2
#    dphidz2 = dphi2/deltaz2
    
    delta3d = math.sqrt((deltax + deltax2)**2 + (deltay + deltay2)**2 + (deltaz + deltaz2)**2 )
    
    if (deltax + deltax2 == 0) :
       phidx = 0
    else:
       phidx =  -(dphidx2 - dphidx)/(deltax + deltax2)
       
       
    if (deltay + deltay2 == 0 ):
       phidy = 0
    else:
       phidy =  -(dphidy2 - dphidy)/(deltay + deltay2)
    
    if (deltaz + deltaz2 == 0 ):
       phidz = 0
    else:
       phidz =  -(dphidz2 - dphidz)/(deltaz + deltaz2)
       phidl = -(dphidz2 - dphidz)/delta3d

# conversion to nC by 1.e9 factor

    if ( numflash[j-1,i] == 0 and numflash[j,i] == 0 and numflash[j+1,i] == 0 ):
       scnetcalc1d[j,i] = 1.e9*eps0*phidz
       scnetcalc3d[j,i] = 1.e9*eps0*phidl
    else:
       scnetcalc1d[j,i] = scnetcalc1d[j-1,i]
       scnetcalc3d[j,i] = scnetcalc3d[j-1,i]

    
    
#    scnetcalc3d[j,i] = 1.e9*eps0*(phidx + phidy + phidz)
    
# do it with electric field    
    exave = 0.5*(ex[j-1,i] + ex[j,i])
    eyave = 0.5*(ey[j-1,i] + ey[j,i])
    ezave = 0.5*(ez[j-1,i] + ez[j,i])

    exave2 = 0.5*(ex[j+1,i] + ex[j,i])
    eyave2 = 0.5*(ey[j+1,i] + ey[j,i])
    ezave2 = 0.5*(ez[j+1,i] + ez[j,i])


#    exave = 0.25*(ex[j-2,i] + 2*ex[j-1,i] + ex[j,i])
#    eyave = 0.25*(ey[j-2,i] + 2*ey[j-1,i] + ey[j,i])
#    ezave = 0.25*(ez[j-2,i] + 2*ez[j-1,i] + ez[j,i])

#    exave2 = 0.25*(ex[j+1,i] + ex[j-1,i] + 2*ex[j,i])
#    eyave2 = 0.25*(ey[j+1,i] + ey[j-1,i] + 2*ey[j,i])
#    ezave2 = 0.25*(ez[j+1,i] + ez[j-1,i] + 2*ez[j,i])
    
    if (deltax + deltax2 == 0) :
       phidx = 0
    else:
       phidx =  2*(exave2 - exave)/(deltax + deltax2)
       
       
    if (deltay + deltay2 == 0 ):
       phidy = 0
    else:
       phidy =  2*(eyave2 - eyave)/(deltay + deltay2)
    
    if (deltaz + deltaz2 == 0 ):
       phidz = 0
    else:
       phidz =   2*(ezave2 - ezave)/(deltaz + deltaz2)

# conversion to nC by 1.e9 factor
#    scnetcalc1d[j,i] = 1.e9*eps0*phidz
#    if (  numflash[j-1,i] == 0 and numflash[j,i] == 0 and numflash[j+1,i] == 0 ):
#       scnetcalc3d[j,i] = 1.e9*eps0*(phidx + phidy + phidz)
#    else:
#       scnetcalc3d[j,i] = scnetcalc3d[j-1,i]
    

print 'min,max chgnetcalc1 = ',np.min(scnetcalc1d),np.max(scnetcalc1d)
print 'min,max chgnetcalc3 = ',np.min(scnetcalc3d),np.max(scnetcalc3d)


# Here we make a multiple sounding plots as separate pages in a pdf file
# This just has the basics. 
# Some things to add might be:
#   The starting location and time for each plot (help keep them straight).
#   See if it is possible to have multiple x axes
#   Try making side-by-side plots to keep the single plot from getting too messy with so many variables
#   Try to generate an analysis value of chgnet using the 1-D Guass's law assumption
#      - ignore points where flashes occurred (

pp = PdfPages(filename+'.soundings_calcphi.pdf')

# use the first line to plot *all* of the soundings
for i in range (ntraj):
#for i in range (70,73):
  zs = zloc[:,i]
  plt.plot(ez[:,i], zloc[:,i], label='Ez')
  plt.plot(airtemp[:,i], zloc[:,i], label='Temperature')
  plt.plot(riserate[:,i], zloc[:,i], label='Rise rate (w+10m/s)')
  plt.plot(10*chgnet[:,i], zloc[:,i], label='Point net charge nCx10') # scale by 10 so it can be seen on the plot
  plt.plot(10*chgnetave[:,i], zloc[:,i], label='Ave. local charge nCx10') # scale by 10 so it can be seen on the plot
  plt.plot(potential[:,i], zloc[:,i], label='Potential (truth) MV')
  plt.plot(potentialcalc[:,i], zloc[:,i], label='Potential (calc 1) MV')
  plt.plot(potentialcalc2[:,i], zloc[:,i],':', label='Potential (calc 2) MV')
  plt.plot(potentialcalc3d[:,i], zloc[:,i],'--', label='Potential (calc 3D) MV')
  plt.plot(10*scnetcalc1d[:,i], zloc[:,i],dashes=[15, 2, 10, 2],  label='1D calc net charge nCx10') # scale by 10 so it can be seen on the plot
  plt.plot(10*scnetcalc3d[:,i], zloc[:,i],dashes=[10,5],  label='3D calc net charge nCx10') # scale by 10 so it can be seen on the plot
  plt.plot(numflash[:,i]*20 - 150, zloc[:,i], label='Flash indicator')
  plt.axis([-150,150,0,16]) # sets axes ranges
  plt.grid(axis='x')
  plt.legend(loc='upper right',fontsize='x-small')
  plt.ylabel("Altitude (km)")
  plt.suptitle('Sounding # '+str(i), fontsize=12)
  plt.savefig(pp, format='pdf')
  plt.close()
pp.close()


xmin = np.min(xloc)
xmax = np.max(xloc)
ymin = np.min(yloc)
ymax = np.max(yloc)
zmax = np.max(zloc)
print 'xmin,ymin,zmax = ',xmin,ymin,zmax
zmin = 0



ncfile.close()


