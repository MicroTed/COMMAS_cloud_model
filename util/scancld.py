#!/usr/bin/env python
#
import string, sys
import numpy as N 
import Ngl

#--------------------------------------------------
def scancldout(filename, var):

      file = open(filename, 'r')
      time = []
      vmax = []
      vmin = []

      while 1:
            line = file.readline()
            if not line: break
            cols = string.split(line)
            if len(cols) >> 1:
                  if cols[0] == 'T':
                        t = eval(cols[2])/60.
                  if cols[0] == var:
                        time.append(t)
                        vmax.append(eval(cols[2]))
                        vmin.append(eval(cols[7]))

      file.close()
      return time, vmax, vmin, len(vmax)
#-------------------------------------------------------------------------------
# Main program for testing...
#
if __name__ == "__main__":

  args = sys.argv[1:]

  if len(args) < 4:

    print 'Insufficient arguments, Usage:  scancld -f filename1 filename2 ... -v W'
    print 'Example:  scancld -f L06M06.out L06M06.out.6 L06M06.out.12 -v Wz-sfc'

  else:
  
    wks = None  

    for i in range(len(args)):
      if args[i] == '-f':
        filepos = i
      if args[i] == '-v':
        varpos = i
      if args[i] == '-w':
        wks = str(args[i+1])
        
    if wks == None:
        wks = "x11"

    files = args[filepos+1:varpos]
    var   = str(args[varpos+1])

    max = []
    min = []

    wks    = Ngl.open_wks(wks,"scancld_plot")  # Open an X11 workstation.
    
    xyExplicitLegendLabels = []

    for i in range(len(files)):
      time, vmax, vmin, npts = scancldout(files[i], var)

      if i == 0:    # initialize arrays....
        dmax = N.zeros((len(files),npts),dtype="f")
        dmin = N.zeros((len(files),npts),dtype="f")

      dmax[i,:] = vmax
      dmin[i,:] = vmin
      
      xyExplicitLegendLabels.append(files[i].split(".")[0])

    print xyExplicitLegendLabels
    resources = Ngl.Resources()
    resources.tiMainString    = "Maximum " + var  # Title for the XY plot
    resources.tiXAxisString = "Time (min)"
    if var == "Wz-sfc":
      resources.tiYAxisString = var + " (x 1000 per sec)"
    elif var == "TH":
      resources.tiYAxisString = var + " (K)"
    elif var == "U" or var == "V" or var == "W":
      resources.tiYAxisString = var + " (m/s)"
    else:
      resources.tiYAxisString = var  
      
    resources.xyLineColors        = ["Red","Blue","Green","Purple"]  # Define line colors.
    resources.xyLineThicknesses   = [2.,2.,2.,2.]    # Define line thicknesses (1.0 is the default).
    
    resources.xyExplicitLegendLabels = xyExplicitLegendLabels
    resources.pmLegendDisplayMode    = "Always"     # Turn on the drawing
    resources.pmLegendZone           = 0            # Change the location
    resources.pmLegendOrthogonalPosF = 0.31         # of the legend
    resources.lgJustification        = "BottomLeft"
    resources.pmLegendWidthF         = 0.3          # Change width and
    resources.pmLegendHeightF        = 0.12         # height of legend.
    resources.pmLegendSide           = "Bottom"        # Change location of
    resources.lgPerimOn              = False        # legend and turn off
                                                  # the perimeter.
    plot = Ngl.xy(wks,time,dmax,resources)    # Draw the plot.

    resources = Ngl.Resources()
    resources.tiMainString    = "Minimum " + var  # Title for the XY plot
    resources.tiXAxisString = "Time (min)"
    if var == "Wz-sfc":
      resources.tiYAxisString = var + " (x 1000 per sec)"
    elif var == "TH":
      resources.tiYAxisString = var + " (K)"
    elif var == "U" or var == "V" or var == "W":
      resources.tiYAxisString = var + " (m/s)"
    else:
      resources.tiYAxisString = var 
      
    resources.xyLineColors        = ["Red","Blue","Green","Purple"]  # Define line colors.
    resources.xyLineThicknesses   = [2.,2.]    # Define line thicknesses (1.0 is the default).
    resources.xyLineThicknesses   = [2.,2.,2.,2.]    # Define line thicknesses (1.0 is the default).
    resources.xyExplicitLegendLabels = xyExplicitLegendLabels

    resources.xyExplicitLegendLabels = xyExplicitLegendLabels
    resources.pmLegendDisplayMode    = "Always"     # Turn on the drawing
    resources.pmLegendZone           = 0            # Change the location
    resources.pmLegendOrthogonalPosF = 0.31         # of the legend
    resources.lgJustification        = "BottomLeft"
    resources.pmLegendWidthF         = 0.3          # Change width and
    resources.pmLegendHeightF        = 0.12         # height of legend.
    resources.pmLegendSide           = "Bottom"        # Change location of
    resources.lgPerimOn              = False        # legend and turn off

    plot = Ngl.xy(wks,time,dmin,resources)    # Draw the plot.

# Clean up.
    del plot 
    del resources

    Ngl.end()

