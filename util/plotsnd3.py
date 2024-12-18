#!/usr/bin/env /Users/Shared/miniconda3/bin/python3
# for EPD:   #!/usr/bin/env /Library/Frameworks/EPD64.framework/Versions/6.3/bin/python
#
import numpy as N
import sys, os
import Nio, string
import Ngl
from optparse import OptionParser

debug = 0
ncdebug = 1

#hodoPlotHeights  = [3000., 4000., 5000., 6000., 7000., 8000., 10000., 12000., 15000.]  # first height is to plot everything below that, then only plot levels above
hodoPlotHeights  = [ 15000.]  # first height is to plot everything below that, then only plot levels above
hodoLabelHeights = [0.0, 1000., 2000., 4000., 6000., 8000., 10000.,12000.,15000.]
hodoLabelOffset  = 1.25

Czero = 273.15
p1013 = 1013.25
Mw    = 18.0160 # molec mass of water
Md    = 28.9660 # molec mass of dry air

#============================================================================
# USER SPECIFIED VARIABLES FOR FILE READING
#============================================================================

#----------------------------------------------------------------------------
# BASE_VARIABLES:    These are the base state variables needed for the skewt
#                    U, V, PI (mean Exner), Theta, Qv (kg/kg)
#----------------------------------------------------------------------------

#base_variables   = ["UBAR","VBAR","PIBAR","THBAR","QVBAR"]
base_variables   = ["UINIT","VINIT","PIINIT","THINIT","QVINIT"]

#----------------------------------------------------------------------------
# THREED_VARIABLES:  These are the 3D variables that can be plotted 
#                    U, V, PI (pert Exner), Theta, Qv (kg/kg)
#----------------------------------------------------------------------------

threed_variables = ["U","V","PI","TH","QV"]

#----------------------------------------------------------------------------
# INPUT_UNITS:       Units of the input for the horizontal location
#                    IF "km" --> then input coordinates are in km
#                    IF ""   --> input coordinates in native grid units
#----------------------------------------------------------------------------

input_units      = "km"

#----------------------------------------------------------------------------
# COORD_VARIABLES:   Coordinate arrays of the model - must be 1D arrays
#----------------------------------------------------------------------------

#coord_variables  = ["xc","yc","zc","time"]
coord_variables  = ["XC","YC","ZC","TIME"]

#----------------------------------------------------------------------------
# CLASS SOUND:
#
# READS AN INPUT FILE FROM COMMAS...
#----------------------------------------------------------------------------
class sound:

    """ Class for reading soundings from COMMAS netcdf files """

    def __init__(self, filename='',  x=0,  y=0, time=0 ):

        if filename == '':
            print("No netcdf file or COMMAS input sounding file (e.g., may20.sound) supplied, please enter:")
            self.filename = sys.stdin.readline()
        else: 
            self.filename = filename

        if x == 0:
            self.x = 0
            self.y = 0
        else: self.x = x

        if y == 0:
            self.x = 0
            self.y = 0
        else: self.y = y

        self.time = time
        self.ti   = 0
        self.xi   = 0
        self.yi   = 0

        if x == 0:  print("Sounding:  No location supplied, base state sounding to be printed")

        self.base_variables   = base_variables
        self.threed_variables = threed_variables

        if input_units != "":
            self.input_units      = input_units
        else: 
            self.input_units      = "M"
            
        self.coord_variables  = coord_variables

        self.xc = None
        self.yc = None
        self.zc = None
        self.tc = None

        self.u  = None
        self.v  = None
        self.th = None
        self.tv = None
        self.qv = None
        self.pi = None
        self.t  = None
        self.td = None
        self.p  = None
        self.z  = None
        self.ri = None

    def info(self):
        print()
        print('Filename:  ',self.filename)
        print('X-location:',self.x, " ",self.input_units)
        print('Y-location:',self.y, " ",self.input_units)
        print('T-location:',self.time, "  seconds")
        print()

        if self.u == None:
            self.read()

        print('Done reading, self.nz = ',self.nz) 
        
        n = self.nz - 1

        print("  Z(m)     U(m/s)   V(m/s)  Pi ")
        print("=================================================================================================")
        while n > -1:
            print("%007.1f      %004.1f     %004.1f   %004.1f"  % \
                  (self.zc[n], self.u[n], self.v[n], self.pi[n] ))
            n -= 1
        print()

        n = self.nz - 1
        print("  Z(m)      P(mb)    T(K)   Td(K)    Theta(K)  Qv(g/kg)   U(m/s)   V(m/s)   SPD      WDIR     RI")
        print("=================================================================================================")
        while n > -1:
            print("%007.1f     %05.1f   %05.1f   %05.1f     %4.1f     %0004.1f       %004.1f     %004.1f    %004.1f     %004.1f    %004.1f" % \
                  (self.zc[n], self.p[n], self.t[n], self.td[n], self.th[n], self.qv[n]*1000., self.u[n], self.v[n], self.spd[n], self.dir[n], self.ri[n]))
            n -= 1
        print()

    def set(self,x=0, y=0, time=0):
        if x:
            self.x = x
        if y:
            self.y = y
        if time:
            self.time = time

        self.read()
        self.info()
        
    def read(self):

        print("Reading from file.....  ", self.filename)

#----------------------------------------------------------------------------
# Read from commas text file

        if self.filename.count(".commas") or self.filename.count(".sound"):      # Input is from a COMMAS model input file
 
            self.psfc, self.tsfc, self.qsfc = N.fromfile(self.filename, dtype='float', sep=" ",count=3)   # Load first line with surface values

            qscale = 1.0
            if self.qsfc >= 0.1:   qscale = 1.0 / 1000.
            self.qsfc = self.qsfc * qscale

            if self.psfc >= 1060.: self.psfc = self.psfc / 100.

            data = N.loadtxt(self.filename, skiprows=1)    # Load the rest of the file, skipping first line (nice function!)

            self.zc = data[:,0]
            self.th = data[:,1]
            self.qv = data[:,2] * qscale
            self.u  = data[:,3]
            self.v  = data[:,4]
            self.nz = N.size(self.zc)
            self.pi = N.zeros((self.nz),dtype='float')
            self.tv = N.zeros((self.nz),dtype='float')
            self.ri = N.zeros((self.nz),dtype='float')

            if debug:
                print()
                print('READING SOUND FILE:  ')
                print('INPUT ZC', self.zc)
                print('INPUT TH', self.th)
                print('INPUT QV', self.qv)
                print('INPUT U', self.u)
                print('INPUT V', self.v)
                print()

            self.tv  = self.th*(1.0 + 0.61*self.qv)

            if self.zc[0] < 1.0:                # if first point in sounding is essentially near ground, psfc --> pisfc
                self.pi[0] = (self.psfc/1000.)**0.285 
            else:
                self.tv[0] = self.tsfc*(1.0 + 0.61*self.qsfc)
                self.pi[0] = (self.psfc/1000.)**0.285 - 4.903/(self.tv[0]*1004.*self.zc[0])

            for n in N.arange(self.nz-1):      # Integrate the pi hydrostatic equation to get pressure
                dz    = self.zc[n+1] - self.zc[n]
                self.pi[n+1] = self.pi[n] - 9.806 / (0.5*(self.tv[n+1]+self.tv[n])*1004./dz)

            for n in N.arange(self.nz-1):      # Compute the Richardson number based on virtual temperature
                dz    = self.zc[n+1] - self.zc[n]
                dtdz  = 9.806 * (self.tv[n+1] - self.tv[n]) / (0.5*(self.tv[n+1] + self.tv[n])*dz)
                ush   = 0.5*(self.u[n+1] - self.u[n])/dz
                vsh   = 0.5*(self.v[n+1] - self.v[n])/dz
                bsh   = ush**2 + vsh**2
                self.ri[n+1] = dtdz / (1.0e-5 + ush**2 + vsh**2)
            self.ri[0] = -999.

            self.title = "File = " + self.filename + "  INPUT SOUNDING"

            self.p  = 1000.*(self.pi[:])**3.508
            self.t  = self.pi[:] * self.th[:]
    
            e       = N.clip( self.qv*self.p/(0.622+self.qv), 1.0e-5, 1000.)
            e       = N.log(e/6.112)
            self.td = Czero + 243.5 / ( 17.67/e - 1.0 )                 # Bolton's approximation        
            self.spd = N.sqrt( self.u[:]**2 + self.v[:]**2 )
            self.dir = 270.0 - 57.29578*N.arctan2(self.v[:], self.u[:])

            self.tc  = self.t - Czero
            self.tdc = self.td - Czero
            print(self.td, e)

#----------------------------------------------------------------------------
# Read from raw2 file from Univ of Wyoming (Text: List) -- add a '.raw2' extension
#  Note: must remove or comment (with #) the header lines
#  If the first data line has 2 values (pressure and height), then also have to 
#   fill out the rest of the values 
# 
#72249 FWD Ft Worth Observations at 12Z 14 Jun 2011
#
#-----------------------------------------------------------------------------
#   PRES   HGHT   TEMP   DWPT   RELH   MIXR   DRCT   SKNT   THTA   THTE   THTV
#    hPa     m      C      C      %    g/kg    deg   knot     K      K      K
#-----------------------------------------------------------------------------
# 1000.0     95   27.2   17.2     54  12.64    180     12  301.3  338.8  303.6
#  989.0    171   27.2   17.2     54  12.64    180     12  301.3  338.8  303.6
#  981.0    246   26.6   16.6     54  12.26    183     16  301.4  337.8  303.6

        elif self.filename.count(".raw2"):      # Input is from a Wyoming sound site input file
 
            print("file = " + self.filename)
            data = N.loadtxt(self.filename, dtype='float')    # Load the file

            self.p   = data[:,0]
            self.tc  = data[:,2] 
            self.tdc = data[:,3] 
            self.t   = data[:,2] + Czero
            self.td  = data[:,3] + Czero
            self.dir = data[:,6]
            self.spd = data[:,7]
            self.zc  = data[:,1]  
            self.qv  = data[:,5] 
#            self.u   = data[:,6] 
#            self.v   = data[:,6] 

            self.nz = N.size(self.zc)
            self.pi = N.zeros((self.nz),dtype='float')
            self.tv = N.zeros((self.nz),dtype='float')
            self.th = N.zeros((self.nz),dtype='float')
            self.ri = N.zeros((self.nz),dtype='float')

            # 0.514 to convert knots to m/s; Also have to rotate angle reference (0 is north, 90 east, 180 south, but for math, 0 is east, and angle increases counterclockwise)
            self.u  = 0.514*self.spd * N.cos(N.deg2rad(-self.dir + 270.))
            self.v  = 0.514*self.spd * N.sin(N.deg2rad(-self.dir + 270.))

            print("dir, spd, u,v = " + str(self.dir[1]) + ", " + str(self.spd[1]) + ", " + str(self.u[1]) + ", " + str(self.v[1]))
            
            self.psfc = self.p[0]
            if self.psfc >= 1060.: self.psfc = self.psfc / 100.

            self.tsfc = self.t[0]
            self.qsfc = self.qv[0]

            qscale = 1.0
            if self.qsfc >= 0.1:   qscale = 1.0 / 1000.
            self.qsfc = self.qsfc * qscale
            self.qv   = self.qv * qscale

            if self.zc[0] > 1.:
                self.zc = self.zc - self.zc[0]

            if self.zc[0] < 1.0:                # if first point in sounding is essentially near ground, psfc --> pisfc
                self.pi[0] = (self.psfc/1000.)**0.285 
            else:
                self.tv[0] = self.tsfc*(1.0 + 0.61*self.qsfc)
                self.pi[0] = (self.psfc/1000.)**0.285 - 4.903/(self.tv[0]*1004.*self.zc[0])


            self.tv[0] = self.t[0] * (1.0 + 0.61*self.qv[0]) / self.pi[0]
            self.th[0] = self.t[0] / self.pi[0]

            for n in N.arange(self.nz-1):                             # Integrate the pi hydrostatic equation to get pressure
                dz           = self.zc[n+1] - self.zc[n]
                self.pi[n+1] = self.pi[n] - 9.806 / (self.tv[n]*1004./dz)  # First guess so I can get Virtual-theta
                self.tv[n+1] = self.t[n+1] * (1.0 + 0.61*self.qv[n+1]) / self.pi[n+1]
                self.pi[n+1] = self.pi[n] - 9.806 / (0.5*(self.tv[n+1]+self.tv[n])*1004./dz)
                self.tv[n+1] = self.t[n+1] * (1.0 + 0.61*self.qv[n+1]) / self.pi[n+1]   # Recompute Virtual-theta
                self.th[n+1] = self.t[n+1] / self.pi[n+1]   # Recompute Virtual-theta

            for n in N.arange(self.nz-1):      # Compute the Richardson number based on virtual temperature
                dz    = self.zc[n+1] - self.zc[n]
                dtdz  = 9.806 * (self.tv[n+1] - self.tv[n]) / (0.5*(self.tv[n+1] + self.tv[n])*dz)
                ush   = 0.5*(self.u[n+1] - self.u[n])/dz
                vsh   = 0.5*(self.v[n+1] - self.v[n])/dz
                bsh   = ush**2 + vsh**2
                self.ri[n+1] = dtdz / (1.0e-5 + ush**2 + vsh**2)
            self.ri[0] = -999.

            self.nz = self.nz-1

            self.title = "File = " + self.filename + "  INPUT SOUNDING"

#----------------------------------------------------------------------------
# Read point forecast sounding from NSSL WRF

        elif self.filename.count(".pfs"):      # Input is from a point forecast sounding
 
            data = N.loadtxt(self.filename, delimiter = ",", dtype='float',skiprows=8)    # Load the file

            self.p   = data[:,0]
            self.tc  = data[:,2] 
            self.tdc = data[:,3] 
            self.t   = data[:,2] + Czero
            self.td  = data[:,3] + Czero
            self.dir = data[:,4]
            self.spd = data[:,5] * 0.514
            self.zc  = data[:,1]  
            
            self.qv  = data[:,6] 

            esat = p1013*10.**(10.79586*(1-Czero/(self.t))-5.02808*N.log10((self.t)/Czero)+ \
                1.50474e-4*(1.-10**(-8.29692*((self.t)/Czero-1.)))+ \
                0.42873e-3*(10.**(4.76955*(1.-Czero/(self.t)))-1.)-2.2195983)

            es   = p1013*10.**(10.79586*(1-Czero/(self.td))-5.02808*N.log10((self.td)/Czero)+ \
                1.50474e-4*(1.-10**(-8.29692*((self.td)/Czero-1.)))+ \
                0.42873e-3*(10.**(4.76955*(1.-Czero/(self.td)))-1.)-2.2195983)

            self.qv  = es * 1000. * Mw/Md / (self.p - esat) 

            self.nz = N.size(self.zc)
            self.pi = N.zeros((self.nz),dtype='float')
            self.tv = N.zeros((self.nz),dtype='float')
            self.th = N.zeros((self.nz),dtype='float')
            self.ri = N.zeros((self.nz),dtype='float')

            self.u  = self.spd * N.cos(N.deg2rad(-self.dir + 270.))
            self.v  = self.spd * N.sin(N.deg2rad(-self.dir + 270.))

            self.psfc = self.p[0]
            if self.psfc >= 1060.: self.psfc = self.psfc / 100.

            self.tsfc = self.t[0]
            self.qsfc = self.qv[0]

            qscale = 1.0
            if self.qsfc >= 0.1:   qscale = 1.0 / 1000.
            self.qsfc = self.qsfc * qscale
            self.qv   = self.qv * qscale

            if self.zc[0] > 1.:
                self.zc = self.zc - self.zc[0]

            if self.zc[0] < 1.0:                # if first point in sounding is essentially near ground, psfc --> pisfc
                self.pi[0] = (self.psfc/1000.)**0.285 
            else:
                self.tv[0] = self.tsfc*(1.0 + 0.61*self.qsfc)
                self.pi[0] = (self.psfc/1000.)**0.285 - 4.903/(self.tv[0]*1004.*self.zc[0])


            self.tv[0] = self.t[0] * (1.0 + 0.61*self.qv[0]) / self.pi[0]
            self.th[0] = self.t[0] / self.pi[0]

            for n in N.arange(self.nz-1):                             # Integrate the pi hydrostatic equation to get pressure
                dz           = self.zc[n+1] - self.zc[n]
                self.pi[n+1] = self.pi[n] - 9.806 / (self.tv[n]*1004./dz)  # First guess so I can get Virtual-theta
                self.tv[n+1] = self.t[n+1] * (1.0 + 0.61*self.qv[n+1]) / self.pi[n+1]
                self.pi[n+1] = self.pi[n] - 9.806 / (0.5*(self.tv[n+1]+self.tv[n])*1004./dz)
                self.tv[n+1] = self.t[n+1] * (1.0 + 0.61*self.qv[n+1]) / self.pi[n+1]   # Recompute Virtual-theta
                self.th[n+1] = self.t[n+1] / self.pi[n+1]   # Recompute Virtual-theta

            for n in N.arange(self.nz-1):      # Compute the Richardson number based on virtual temperature
                dz    = self.zc[n+1] - self.zc[n]
                dtdz  = 9.806 * (self.tv[n+1] - self.tv[n]) / (0.5*(self.tv[n+1] + self.tv[n])*dz)
                ush   = 0.5*(self.u[n+1] - self.u[n])/dz
                vsh   = 0.5*(self.v[n+1] - self.v[n])/dz
                bsh   = ush**2 + vsh**2
                self.ri[n+1] = dtdz / (1.0e-5 + ush**2 + vsh**2)
            self.ri[0] = -999.

            self.nz = self.nz-1

            self.title = "File = " + self.filename + "  INPUT SOUNDING"

#----------------------------------------------------------------------------
# Read from commas netcdf file

        elif self.filename.count(".nc"):
                
            print('reading netcdf file')
            
            cdf_file = Nio.open_file(self.filename, "r")
            
          #  print cdf_file

            self.xc = cdf_file.variables[self.coord_variables[0]][:]
            self.nx = N.size(self.xc)

            self.yc = cdf_file.variables[self.coord_variables[1]][:]
            self.ny = N.size(self.yc)

            self.zc = cdf_file.variables[self.coord_variables[2]][:]
            self.nz = N.size(self.zc)

            self.tc = cdf_file.variables[self.coord_variables[3]][:]
            self.nt = N.size(self.tc)

            self.ri = N.zeros((self.nz),dtype='float')

            print('self.zc = ',self.zc)
            
            
            if self.input_units == "km":
                self.xc = self.xc/1000.
                self.yc = self.yc/1000.

            if self.x == 0:                             # We plot the base state sounding....

                self.u  = cdf_file.variables[self.base_variables[0]][:]
                self.v  = cdf_file.variables[self.base_variables[1]][:]
                self.pi = cdf_file.variables[self.base_variables[2]][:]
                self.th = cdf_file.variables[self.base_variables[3]][:]
                self.qv = cdf_file.variables[self.base_variables[4]][:]
                self.tv  = self.th*(1.0 + 0.61*self.qv)

                self.title = "File = " + self.filename + "  BASE STATE SOUNDING"

                if ncdebug:
                   print() 
                   print('sound.read() BASE STATE')
                   print('self.x= ',self.x, 'self.y= ',self.y, 'self.t= ',self.t, 'self.time =', options.time)
                   print('options.xc = ',options.xc)
                   print(self.xi, self.yi, self.ti)
                   print('self.xindex= ',self.x, 'self.yindex= ',self.y, 'self.tindex= ',self.t)

                for n in N.arange(self.nz-1):      # Compute the Richardson number based on virtual temperature
                    dz    = self.zc[n+1] - self.zc[n]
                    dtdz  = 9.806 * (self.tv[n+1] - self.tv[n]) / (0.5*(self.tv[n+1] + self.tv[n])*dz)
                    ush   = 0.5*(self.u[n+1] - self.u[n])/dz
                    vsh   = 0.5*(self.v[n+1] - self.v[n])/dz
                    bsh   = ush**2 + vsh**2
#                    self.spd[n] = 
                    self.ri[n+1] = dtdz / (1.0e-5 + ush**2 + vsh**2)
                    self.title = "File = " + self.filename  \
                           + "  X=" + string.zfill(int(self.x),2) + " " + self.input_units \
                           + "  Y=" + string.zfill(int(self.y),2) + " " + self.input_units \
                           + "  T=" + string.zfill(int(self.time),4) + " s"

            else:                                       # Get the time and location and plot the sounding...

                x0 = N.nonzero( N.where(self.xc <= self.x,True,False) )
                y0 = N.nonzero( N.where(self.yc <= self.y,True,False) )
                t0 = N.nonzero( N.where(self.tc <= self.time,True,False) )

                self.xi = x0[-1][-1]
                self.x  = self.xc[x0[-1][-1]]
                self.yi = y0[-1][-1]
                self.y  = self.yc[y0[-1][-1]]
                self.ti = t0[-1][-1]           # this syntax first extracts the array from the tuple, then the last index from array
                self.t  = self.tc[t0[-1][-1]]  # this syntax first extracts the array from the tuple, then the last index from array

                if ncdebug:
                   print() 
                   print('sound.read() LOCATION')
                   print('self.x= ',self.x, 'self.y= ',self.y, 'self.t= ',self.t, 'self.time =', self.time)
                   print(self.xi, self.yi, self.ti)
                   print('self.xindex= ',self.xi, 'self.yindex= ',self.yi, 'self.tindex= ',self.ti)
                   print('size of base var (time): ',N.size(cdf_file.variables[self.base_variables[2]][:,0]))
                   print('size of base var (height): ',N.size(cdf_file.variables[self.base_variables[2]][self.ti,:]))

                
                
                self.u  = cdf_file.variables[self.threed_variables[0]][self.ti,:,self.yi,self.xi] 
                self.v  = cdf_file.variables[self.threed_variables[1]][self.ti,:,self.yi,self.xi]
                self.pi = cdf_file.variables[self.threed_variables[2]][self.ti,:,self.yi,self.xi]  + cdf_file.variables[self.base_variables[2]][self.ti,:]
                self.th = cdf_file.variables[self.threed_variables[3]][self.ti,:,self.yi,self.xi]
                self.qv = cdf_file.variables[self.threed_variables[4]][self.ti,:,self.yi,self.xi]
                self.tv  = self.th*(1.0 + 0.61*self.qv)
                
                print('self.u = ',self.u)
                print('self.v = ',self.v)
                print('self.pi = ',self.pi)

                for n in N.arange(self.nz-1):      # Compute the Richardson number based on virtual temperature
                    dz    = self.zc[n+1] - self.zc[n]
                    dtdz  = 9.806 * (self.tv[n+1] - self.tv[n]) / (0.5*(self.tv[n+1] + self.tv[n])*dz)
                    ush   = 0.5*(self.u[n+1] - self.u[n])/dz
                    vsh   = 0.5*(self.v[n+1] - self.v[n])/dz
                    bsh   = ush**2 + vsh**2
                    self.ri[n+1] = dtdz / (1.0e-5 + ush**2 + vsh**2)
                    self.title = "File = " + self.filename  \
                           + "  X=" + string.zfill(int(self.x),2) + " " + self.input_units \
                           + "  Y=" + string.zfill(int(self.y),2) + " " + self.input_units \
                           + "  T=" + string.zfill(int(self.time),4) + " s"

            self.p  = 1000.*(self.pi[:])**3.508
            self.t  = self.pi[:] * self.th[:]
    
            e       = N.clip( self.qv*self.p/(0.622+0.001*self.qv), 1.0e-6, 1000.)
            e       = N.log(e/6.112)
            self.td = N.minimum(self.t,Czero + 243.5 / ( 17.67/e - 1.0 )    )             # Bolton's approximation        
            print('self.td (Bolton) = ',self.td)

            e       = N.clip( self.qv, 1.0e-6, 1000.)*self.p*100/380.
            e       = N.log(e)
            # self.td = N.minimum(self.t,Czero + (35.86*e - 17.27*273.15)/(e - 17.27)    )             # 
            self.td = (35.86*e - 17.27*273.15)/(e - 17.27)           # 
            self.td = N.minimum(self.t,self.td)
            print('self.td = ',self.td)
            print('self.qv = ',self.qv)
            print('self.p = ',self.p)


            self.spd = N.sqrt( self.u[:]**2 + self.v[:]**2 )
            self.dir = 270.0 - 57.29578*N.arctan2(self.v[:], self.u[:])
            self.tc  = self.t - Czero
            self.tdc = self.td - Czero

            self.psfc = self.p[0]
            print('self.psfc = ', self.psfc)
            if self.psfc >= 1060.: self.psfc = self.psfc / 100.

            self.tsfc = self.t[0]
            self.qsfc = self.qv[0]

            print('done reading netcdf file')

#----------------------------------------------------------------------------
# Input sounding format is not supported...

        else:   # PROBLEM -> dont know what type of file you are reading!!!

            print("SOUND.READ:  FILE SUFFIX IS NOT RECOGNIZED")
            print("SOUND.READ:  NO INPUT")
            print("SOUND.READ:  ERROR")
            print("SOUND.READ:  Returning -1")
            sys.exit(-1)
            
# Debug code....

        if debug:
            print('Time Index:  ', self.ti)
            print('X-Index:     ', self.xi)
            print('Y-Index:     ', self.yi)
            print('U:  ', self.u[:])
            print('V:  ', self.v[:])
            print('PI:  ', self.pi[:])
            print('TH:  ', self.th[:])
            print('QV:  ', self.qv[:])

        if debug:
            print('Temp: ', self.t[:])
            print('Dewpt:', self.td[:])
            print('SPD:  ', self.spd[:])
            print('DIR:  ', self.dir[:])

    def write(self, filename=None):

        if self.filename == None:
            filename = "output.sound"
        else:
            filename = self.filename+".output.sound"

        if debug:  print("Writing file.....  ", filename)

        N.savetxt(filename, N.transpose((self.zc, self.th, self.qv*1000., self.u, self.v)), fmt='%8.2f')
        f = open(filename, 'r')
        text = f.read()
        f.close()
        f = open(filename, 'w')
        str = "%8.2f    %8.2f   %8.2f\n" % (self.psfc, self.tsfc, 1000.*self.qsfc)
        f.write(str)
        f.write(text)
        f.close()

    def skewt(self):

#        if self.u == None:
#            self.read()

        skewtOpts                              = Ngl.Resources()
        skewtOpts.sktTemperatureUnits          = "celsius"   # default is "fahrenheit"

        skewtOpts.sktHeightScaleOn             = True      # default is False
        skewtOpts.sktHeightScaleUnits          = "km"    # default is "feet"
        skewtOpts.sktColoredBandsOn            = False      # default is False
        skewtOpts.sktGeopotentialWindBarbColor = "black"
        skewtOpts.tiMainString                 = self.title
        skewtOpts.sktStdAtmosphereLineOn       = False # default is True

        dataOpts                               = Ngl.Resources()  # Options describing 
                                                                  # data and plotting.

        dataOpts.sktPressureWindBarbStride     = 2                # Plot every 2nd wind barb
        
# It seems the only way to force portrait orientation is before the open_wks command in the main section
#        skewtOpts.wkOrientation          = "portrait"
#        dataOpts.wkOrientation          = "portrait"
#        rlist                       = Ngl.Resources()
#        rlist.wkOrientation         = "portrait"
#        rlist.nglPaperOrientation         = "portrait"
#        skewtOpts.nglPaperOrientation   = "portrait"
#        dataOpts.nglPaperOrientation   = "portrait"
#        Ngl.set_values(a.wks,rlist)
        
        skewt_bkgd = Ngl.skewt_bkg(self.wks, skewtOpts)
        skewt_data = Ngl.skewt_plt(self.wks, skewt_bkgd, self.p[:], self.tc[:], self.tdc[:], self.zc[:],  self.spd[:], self.dir[:], dataOpts)
        Ngl.draw(skewt_bkgd)
        Ngl.draw(skewt_data)
        Ngl.frame(self.wks)

    def hodo(self):

#        if self.u == None:
#            self.read()

        hodoOpts                              = Ngl.Resources()

        maxwind = 1.1 * self.spd.max()
        circFr  = 7.5

# Make the hodograph background plot circles - from Shea wind rose
        MaxNumCircles = int(N.fix(maxwind/circFr + 1.))
        FrMaxNum      = circFr*MaxNumCircles   # max circle radius
        nCirc = N.arange(361)
        xCirc = N.zeros((MaxNumCircles,361), dtype='float')
        yCirc = N.zeros((MaxNumCircles,361), dtype='float')
        rad   = 4.*N.arctan(1.0)/180.            # degress to radians
        xcos  = N.cos(nCirc[:]*rad)*circFr
        xsin  = N.sin(nCirc[:]*rad)*circFr
        for n in N.arange(MaxNumCircles):        # plot coordinates for freq circles
            xCirc[n,:] = (n+1)*xcos[:]
            yCirc[n,:] = (n+1)*xsin[:]

# Specify data limits for X and Y axes.
        extraSpace          = max(3.,circFr/3.)         # Extra space beyond outer circle
        hodoOpts.trXMinF    = min(-10,self.u.min() - extraSpace)  # min X 
        hodoOpts.trXMaxF    = self.u.max() + extraSpace  # max X
        hodoOpts.trYMinF    = min(-10,self.v.min() - extraSpace)  # min X 
        hodoOpts.trYMaxF    =  self.v.max() + extraSpace  # max X
        hodoOpts.tmXTOn     = False                       # Tick marks? 
        hodoOpts.tmXBOn     = True 
        hodoOpts.tmYLOn     = True
        hodoOpts.tmYROn     = False
        hodoOpts.vpXF       = 0.1
        hodoOpts.vpYF       = 0.9
        hodoOpts.vpWidthF   = 0.7
        hodoOpts.vpHeightF  = 0.7*((hodoOpts.trYMaxF-hodoOpts.trYMinF) /  (hodoOpts.trXMaxF-hodoOpts.trXMinF))
        hodoOpts.tiMainFont = 21                             # Helvetica
        hodoOpts.tiMainFontHeightF    = 0.016
        hodoOpts.tiXAxisString        = "U winds m s~S~-1~N~"
        hodoOpts.tiXAxisFont          = 21
        hodoOpts.tiXAxisFontHeightF   = 0.018
        hodoOpts.tmXBLabelFontHeightF = 0.018
        hodoOpts.tmXBLabelFont        = 21
        hodoOpts.tiYAxisString        = "V winds m s~S~-1~N~"
        hodoOpts.tiYAxisFont          = 21
        hodoOpts.tiYAxisFontHeightF   = 0.018
        hodoOpts.tmYLLabelFontHeightF = 0.018
        hodoOpts.tmYLLabelFont        = 21
        hodoOpts.nglDraw              = False
        hodoOpts.nglFrame             = False
        hodoOpts.xyMonoDashPattern    = True                 # set all circles to solid
        hodoOpts.xyMonoLineThickness  = True
        hodoOpts.xyLineThicknessF     = 1.0
        hodoOpts.tiMainString         = self.title

# Draw circles

        hodo_bkgd = Ngl.xy(self.wks,xCirc,yCirc,hodoOpts)   # Draw circles

# Draw base lines

#       hodoOpts                   = Ngl.Resources()
#       hodoOpts.gsLineDashPattern = 0            # solid
#       hodoOpts.gsLineThicknessF  = 2.0
#       N.polyline (wks,hodo_bkgd,[-FrMaxNum-extraSpace, FrMaxNum+extraSpace], [0.0, 0.0], hodoOpts)
#       N.polyline (wks,hodo_bkgd,[0.0,0.0],[-FrMaxNum-extraSpace, FrMaxNum+extraSpace], hodoOpts)

# Draw hodograph line

        hodoOpts                      = Ngl.Resources()
        hodoOpts.gsLineThicknessF     = 2.0
        hodoOpts.gsLineColor          = "red" 

        n = self.zc.searchsorted(hodoPlotHeights[0])
        Ngl.add_polyline(self.wks,hodo_bkgd,self.u[:n], self.v[:n], hodoOpts)

        up, vp = [self.u[n-1]], [self.v[n-1]]
        for z in hodoPlotHeights[1:]:
            m = self.zc.searchsorted(z)
            up.append(self.u[m])
            vp.append(self.v[m])

        Ngl.add_polyline(self.wks,hodo_bkgd,up, vp, hodoOpts)

# Add in height markers

        hodoOpts                   = Ngl.Resources()
        hodoOpts.gsMarkerIndex     = 16       # dots
        hodoOpts.gsMarkerColor     = "Blue"
        hodoOpts.gsMarkerSizeF     = 0.014    # twice normal size

        hodoTxRes                  = Ngl.Resources()
        hodoTxRes.txFont           = 21                  # Change the default font.
        hodoTxRes.txFontHeightF    = 0.02                # Set the font height.
         
        for z in hodoLabelHeights:
            n = self.zc.searchsorted(z)
            Ngl.add_polymarker(self.wks, hodo_bkgd, self.u[n], self.v[n], hodoOpts)
            if z < 1000.:
                label = "%d m" % N.fix(self.zc[0]) 
            else:
                label = "%d km" % N.fix(z/1000.)

            Ngl.text(self.wks,hodo_bkgd,label, self.u[n]-hodoLabelOffset, self.v[n]+hodoLabelOffset, hodoTxRes) # Label plot

        Ngl.draw(hodo_bkgd)
        Ngl.frame(self.wks)

#-------------------------------------------------------------------------------
# Main program for testing...
# 
if __name__ == "__main__":

    if debug:
        a = sound(filename='./20110614.22Z/20110614.22Z.C34.pfs')
        a.info()
        a.wks = Ngl.open_wks('x11', "debug.snd" )
        a.skewt()
        a.hodo()
        a.write()

    else:

        usage = "usage: %prog [options] arg"
        parser = OptionParser(usage)
        parser.add_option("-f", "--file", dest="file", type="string", help="netCDF or text input file to be plotted, e.g., 'may20.000.nc'")
        parser.add_option("-o", "--out",  dest="out",  type="string", help="Name of output sounding file")
        parser.add_option("-w", "--wks",  dest="wks",  type="string", help="workstation type: x11, ps, ncgm, or pdf (default)")
        parser.add_option("-x", "--xc",   dest="xc",   type="int",    help="x-coordinate (km) of sounding to be plotted / no coordinate means base state plotted")
        parser.add_option("-y", "--yc",   dest="yc",   type="int",    help="y-coordinate (km) of sounding to be plotted / no coordinate means base state plotted")
        parser.add_option("-t", "--time", dest="time", type="float",  help="Time (seconds) of plot / default t=0 plotted")
        parser.add_option("-a", "--all",  dest="all",  action="store_true", help="Flag to plot all times in the file")
        parser.add_option("-c", "--commas",  dest="commas",  action="store_true", help="Boolean flag to write out commas file with prefix.commas")

        (options, args) = parser.parse_args()

        options.dir = os.getcwd()

        if options.file == None:
            print()
            parser.print_help()
            print()
            print("ERROR:  file not defined...EXITING PLOTTING")
            print()
            sys.exit(0)

        if options.wks == None:
            print()
            print("Setting workstation to PDF")
            print()
            options.wks = "pdf"
#            options.wks = "pdf"

        if options.xc == None:
            print()
            print("X-coordinate not specified, plotting base state sounding")
            print()
            options.xc = 0

        if options.yc == None:
            print()
            print("y-coordinate not specified, plotting base state sounding")
            print()
            options.yc = 0

        if options.time == None:
            options.time = 0.0

        a = sound(filename=options.file,x=options.xc,y=options.yc,time=options.time)
        print("RUN INFO")
        a.info()
        print("SET WKS")
        res = Ngl.Resources()
        res.wkOrientation = "portrait" # this forces portrait orientation
        a.wks = Ngl.open_wks(options.wks, options.file.split(".pfs")[0],res )
        print("RUN SKEWT")
        a.skewt()
        print("RUN HODO")
        a.hodo()
        Ngl.end()
        

        if options.commas:
          commasfile = options.file.split("pfs")[0] + "commas"
          a.write(filename=commasfile)
        else:
          a.write(filename=options.out)

# Read some version of the sounding data in to set
