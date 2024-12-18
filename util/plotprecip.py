#!/opt/local/anaconda/bin/python
# #!/opt/local/anaconda.updated/bin/python

# Purpose: Plot time series (bin.ion.pl) from multiple simulations

import csv
import sys
import pandas as pd

import numpy as np
import numpy.ma as ma
import math
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib import rc
#from matplotlib.ticker import ScalarFormatter

# class OOMFormatter(ticker.ScalarFormatter):
#     def __init__(self, order=0, fformat="%1.1f", offset=True, mathText=True):
#         self.oom = order
#         self.fformat = fformat
#         ticker.ScalarFormatter.__init__(self,useOffset=offset,useMathText=mathText)
#     def _set_orderOfMagnitude(self, nothing):
#         self.orderOfMagnitude = self.oom
#     def _set_format(self, vmin, vmax):
#         self.format = self.fformat
#         if self._useMathText:
#             self.format = '$%s$' % self.format # ticker._mathdefault(self.format)


rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

runhead = 'may29zvh1000cn'
tail = 'vtrr'

#runhead = 'may29tak1000cn'
#tail = 'vtr'

datatype = 3 # 1 = bin.ion; 2 = coldpool; 3 = precip; 4 = binproc
# set index
# datatype 1 (bin) : 21=graupel mass; 26= hail mass; 66 = graupel density
# datatype 2 (coldpool) : 1=T-ave, 2= T-area, 3 = T-sum, 4 = Theta-V-ave,  5= Theta-V-area, 6 = Theta-V-sum
# datatype 3 (precip): 2=rain rate ; 3 = hail rate, 4=rain accum. ; 5 = hail accum
# datatype 4 (binproc) : 19=rime density; 22=qcond; 23=rain evap; 37=total net riming; 38=gr+hl melting
index = 4

if ( datatype == 1 ):
   suffix = '.bin'
elif ( datatype == 2 ):
   suffix = '.coldpool'
elif ( datatype == 3 ):
   suffix = '.precip'
elif ( datatype == 4 ):
   suffix = '.binproc'

runtail =  tail + suffix + '.txt'

outputhead = runhead + 'XX' + tail 

vtlist = ['07','08','09','10','11','12','13']
label_list = ['0.7 Vt','0.8 Vt','0.9 Vt','1.0 Vt','1.1 Vt','1.2 Vt','1.3 Vt']
color_list = ['skyblue','dodgerblue','blue','black','red','lightcoral','violet']
#color_list = ['skyblue','dodgerblue','blue','black','crimson','red','violet']
files = []
numfiles = len(vtlist)

# overlay says whether to read in a second set to plot with dashed lines
overlay = 4
if ( overlay == 1 ):
   runhead2 = '../may29zvh1000cn6cdxrr/may29zvh1000cn' # icdx=6, variable density
   tail2 = 'vt6cdxrr'
   run2 = ''
elif ( overlay == 2 ):
   runhead2 = '../may29zvh1000cnvtgrr/may29zvh1000cn' # only graupelfallfac is changed
   tail2 = 'vtgrr'
   run2 = ''
elif ( overlay == 3 ):
   runhead2 = '../may29zh1000cnvtrr500h850hl/may29zh1000cn' # fixed density
   tail2 = 'vtrr500h850hl'
   run2 = 'zhXX500h850hl'
elif ( overlay == 4 ):
   runhead2 = '../may29tak1000cnvtr/may29tak1000cn' # Takahashi bin
   tail2 = 'vtr'
   run2 = 'takXX'
elif ( overlay == 5 ):
   runhead2 = '../vtxfso/may29zvh1000cn' # ifallsedonly = 1 (1fso) (sedimentation only)
   tail2 = 'vt1fso'
   run2 = 'zvh1fso'
elif ( overlay == 6 ):
   runhead2 = '../vtxfso/may29zvh1000cn' # ifallsedonly = 2 (2fso) (micro rates only)
   tail2 = 'vt2fso'
   run2 = 'zvh2fso'
else:
   tail2 = ''
   run2 = ''
   runhead2 = ''

runtail2 =  tail2 + suffix + '.txt'
files2 = []
if (overlay > 0):
    outputhead = outputhead + '.' + run2 + tail2



# set output file name depending on index

datscale = 1
if ( datatype == 1 ):
   if ( index == 21 ): # graupel mass
       pp = outputhead + '.bin.graupelmass.pdf'
   elif ( index == 26 ): # hail mass
       pp = outputhead + '.bin.hailmass.pdf'
   elif ( index == 29 ): # cloud ice mass
       pp = outputhead + '.bin.cloudicemass.pdf'
   elif ( index == 30 ): # droplet mass
       pp = outputhead + '.bin.dropletmass.pdf'
   elif ( index == 31 ): # rain mass
       pp = outputhead + '.bin.rainmass.pdf'
   elif ( index == 34 ): # snow mass
       pp = outputhead + '.bin.snowmass.pdf'
   elif ( index == 53 ): # updraft volume w > 10
       pp = outputhead + '.bin.udv10.pdf'
   elif ( index == 66 ): # graupel density
       pp = outputhead + '.bin.graupeldens.pdf'
   elif ( index == 67 ): # hail density
       pp = outputhead + '.bin.haildens.pdf'
   else:
       print 'index not defined yet'
       sys.exit
elif ( datatype == 2 ):
   if ( index == 3 ): # t-sum
       pp = outputhead + '.tsum.pdf'
       datscale = 1000
   elif ( index == 6 ): # thetav-sum
       pp = outputhead + '.thetavsum.pdf'
       datscale = 1000
   else:
       print 'index not defined yet'
       sys.exit
elif ( datatype == 3 ):
   if ( index == 4 ): # 
       pp = outputhead + '.rainaccum.pdf'
   elif ( index == 5 ): # 
       pp = outputhead + '.hailaccum.pdf'
   else:
       print 'index not defined yet'
       sys.exit

elif ( datatype == 4 ):
   # print 'set pp for datatype 4'
   if ( index == 1 ): # crfrzf
       pp = outputhead + '.binproc.crfrzf.pdf'
   elif ( index == 19 ): # rime density
       pp = outputhead + '.binproc.rimedens.pdf'
   elif ( index == 22 ): # pcond
       pp = outputhead + '.binproc.qcond.pdf'
   elif ( index == 23 ): # evaporation (pevap)
       pp = outputhead + '.binproc.pevap.pdf'
   elif ( index == 24 ): # pmlt
       pp = outputhead + '.binproc.pmlt.pdf'
   elif ( index == 25 ): # pfrz
       pp = outputhead + '.binproc.pfrz.pdf'
   elif ( index == 37 ): # net riming (qhacwrsh+qhlacwrsh)
       pp = outputhead + '.binproc.rimetot.pdf'
   elif ( index == 38 ): #  qh+hlmlr
       pp = outputhead + '.binproc.grhlmlt.pdf'
#   elif ( index == 66 ): # graupel density
#       pp = outputhead + '.binproc.graupeldens.pdf'
#   elif ( index == 67 ): # hail density
#       pp = outputhead + '.binproc.haildens.pdf'
   else:
       print 'index not defined yet'
       sys.exit

else:
    print 'datatype not defined yet'
    sys.exit

# print 'after setting pp'

#pp = PdfPages(runhead+'XXvtrr'+'.bin.hailmass.pdf')

plt.style.use('presentation')
# print plt.style.available

# create list of file names using vtlist
for num in vtlist:
   name = runhead+num+runtail
   files.append(name)
   if ( overlay > 0 ):
       name = runhead2+num+runtail2
       files2.append(name)
       

# create list 'alldat' to hold data from each file
alldat = []
alldat2 = []

# get data from all the files in the list
for filename in files:
   #f = open(filename ,  'r' )
   #reader1 = csv.reader(f,delimiter=',',skipinitialspace=True)
   
   #header1 = reader1.next()
   #header2 = reader1.next()
   #header3 = reader1.next()
   #print header1
   #print header2
   #print header3
   
   #titles = reader1.next()
   #numtitles = len(titles)
   #print numtitles
   #print titles
   #print titles[0]
   #for i in range(numtitles):
   #    print "title " + str(i) + ' = ' + titles[i]
   
   # now read the data:
   
   #f.close
   
   #dat = pd.read_table(f,delim_whitespace=True,engine='python')
   
   # header=3 tells it to skip the first 3 lines
   if ( datatype == 1 ):
      dat = pd.read_csv(filename,header=3,sep=',')
   elif ( datatype == 4 ):
      dat = pd.read_csv(filename,header=1,sep=',')
   else:
      dat = pd.read_csv(filename,header=0,sep=',')
   

   #print (dat.head(3))
   #print (dat.tail(3))
   
   alldat.append(dat)

if ( overlay > 0 ):
   for filename in files2:
      if ( datatype == 1 ):
         dat = pd.read_csv(filename,header=3,sep=',')
      elif ( datatype == 4 ):
         dat = pd.read_csv(filename,header=1,sep=',')
      else:
         dat = pd.read_csv(filename,header=0,sep=',')
      alldat2.append(dat)


titles = alldat[0].columns
numtitles = len(titles)
print titles[index]
#print numtitles

#dat1 =  alldat[0]
#time_min = dat1.iloc[:,0]
#graupelmass = dat1.iloc[:,index]
#print time_min
#print graupelmass

graupelmass2 = alldat[1].iloc[:,index]

if ( datatype == 1 or datatype == 2 or datatype == 4):
   time_min = alldat[0].iloc[:,0]
else:
   time_min = alldat[0].iloc[:,1]

ax = plt.axes()
for i in range(numfiles):
    graupelmass = alldat[i].iloc[:,index]
    plt.plot(time_min,graupelmass,label=label_list[i], linewidth=2, color=color_list[i])

if ( overlay > 0 ): 
    for i in range(numfiles):
        graupelmass = alldat2[i].iloc[:,index]
        plt.plot(time_min,graupelmass,'--',label=None, linewidth=2, color=color_list[i])

if ( datatype == 1 ):
   if ( index == 53 ):
       plt.ylabel('Updraft Volume ($>$ 10m/s) (km$^{3}$)')
   elif ( index == 31 ):
       plt.ylabel('Rain Mass (kg)')
   elif ( index == 34): # 34 = Snow Mass
       plt.ylabel('Snow Mass (kg)')
   elif ( index == 29): # 29 = Cloud ice Mass
       plt.ylabel('Cloud Ice Mass (kg)')
   elif ( index == 30): # 30 = Cloud droplet Mass
       plt.ylabel('Cloud Droplet Mass (kg)')
   elif ( index == 66): # graupel density
       plt.ylabel('Mean Graupel Density (kg\,m$^{-3}$)')
   elif ( index == 67): # hail density
       plt.ylabel('Mean Hail Density (kg\,m$^{-3}$)')
   else:
       plt.ylabel(titles[index])

elif ( datatype == 2 ): # coldpool
   if ( index == 3 ):
       plt.ylabel(r'Total ${\theta}$ Pert (K)') # use "r" at beginning to specify "raw" text so that \t is not changed to a tab. Otherwise can use \\theta to escape the \.
   elif ( index == 6 ):
       plt.ylabel(r'Total ${\theta}_{\rm v}$ Pert (K)')
   else:
       plt.ylabel(titles[index])

elif ( datatype == 3 ): # precip
   if ( index == 4 ):
       plt.ylabel(r'Rain Accumulation (kg)') # use "r" at beginning to specify "raw" text so that \t is not changed to a tab. Otherwise can use \\theta to escape the \.
   elif ( index == 5 ):
       plt.ylabel(r'Hail Accumulation (kg)') # use "r" at beginning to specify "raw" text so that \t is not changed to a tab. Otherwise can use \\theta to escape the \.
   else:
       plt.ylabel(titles[index])

elif ( datatype == 4 ): # binproc
   if ( index == 1 ):
       plt.ylabel(r'Drop freezing (#/min)') # use "r" at beginning to specify "raw" text so that \t is not changed to a tab. Otherwise can use \\theta to escape the \.
   elif ( index == 19 ):
       plt.ylabel(r'Avg. Rime Density (kg\,m$^{-3}$)') # use "r" at beginning to specify "raw" text so that \t is not changed to a tab. Otherwise can use \\theta to escape the \.
   elif ( index == 22 ):
       plt.ylabel(r'Condensation (kg/min)') # use "r" at beginning to specify "raw" text so that \t is not changed to a tab. Otherwise can use \\theta to escape the \.
   elif ( index == 23 ):
       plt.ylabel(r'Rain Evaporation (kg/min)') # use "r" at beginning to specify "raw" text so that \t is not changed to a tab. Otherwise can use \\theta to escape the \.
   elif ( index == 37 ):
       plt.ylabel(r'Graupel + Hail Net Riming (kg/min)') # use "r" at beginning to specify "raw" text so that \t is not changed to a tab. Otherwise can use \\theta to escape the \.
   elif ( index == 38 ):
       plt.ylabel(r'Graupel + Hail Melting (kg/min)') # use "r" at beginning to specify "raw" text so that \t is not changed to a tab. Otherwise can use \\theta to escape the \.
   else:
       plt.ylabel(titles[index])

else:
    plt.ylabel(titles[index])

if (  datscale != 1 ):
   #scale_y = datscale
   #ticks_y = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/scale_y))
   #ax.yaxis.set_major_formatter(ticks_y)
   #sf = ticker.ScalarFormatter()
   #sf.set_scientific(True)
   #sf.set_powerlimits((-1, 1))
   #print 'powerlim: ', sf.get_powerlimits()
   #ticks_y = ticker.EngFormatter(unit='K')
   ticks_y = ticker.ScalarFormatter(useOffset=False)
   ticks_y.set_scientific(True)
   ticks_y.set_powerlimits((-1, 3))
   #print 'useoffset: ', ticks_y.get_useOffset()
   ax.yaxis.set_major_formatter(ticks_y)
   #ax.yaxis.set_major_formatter(OOMFormatter(3))
   dum = ticks_y.get_offset()
   #plt.tight_layout()
   #dum = ax.yaxis.get_offset_text().get_text()
   #print 'offset = ',dum
   

plt.xlabel("Time (min)")
xmin, xmax, ymin, ymax = plt.axis()
plt.axis([20,180,ymin, ymax])
if ( index == 66 or index == 67 ):
    plt.axis([20,180,250, 950])
plt.minorticks_on()
#plt.legend(loc='upper left',bbox_to_anchor=(1,1),fontsize='x-small')
plt.legend(loc='best',fontsize='large')


#plt.plot(time_min,graupelmass)
### plt.axis([-150,150,0,16])

plt.savefig(pp, format='pdf')
plt.show()
plt.close()


# title 0 = Time (min)
# title 1 = IC Rate
# title 2 = IC/min
# title 3 = +CG Rate
# title 4 = -CG Rate
# title 5 = Tries
# title 6 = ICDIS
# title 7 = +CG charge
# title 8 = -CG charge
# title 9 = +CG cum. charge
# title 10 = -CG cum. charge
# title 11 = NONI. neg
# title 12 = NONI. pos
# title 13 = IND. neg
# title 14 = IND. pos.
# title 15 = W-max (m/s)
# title 16 = E-max (kV/m)
# title 17 = Net Charge
# title 18 = Net Pos. Chg.
# title 19 = Net Neg. Chg.
# title 20 = Graupel Volume (km**3)
# title 21 = Graupel Mass (kg)
# title 22 = Graupel Mass T < 0
# title 23 = Graupel Mass T > 0
# title 24 = Graupel Vol. T < 0
# title 25 = Graupel Vol. T > 0
# title 26 = Hail Mass (kg)
# title 27 = Hail Vol. T < 0
# title 28 = Hail Vol. T > 0
# title 29 = Cloud Ice Mass
# title 30 = Cloud Droplet Mass
# title 31 = Rain Mass Tot
# title 32 = Rain Mass (T < 0)
# title 33 = Rain Mass (T > 0)
# title 34 = Snow Mass
# title 35 = G-I Neg
# title 36 = G-I Pos
# title 37 = G-S Neg
# title 38 = G-S Pos
# title 39 = H-I Neg
# title 40 = H-I Pos
# title 41 = H-S Neg
# title 42 = H-S Pos
# title 43 = Updraft Mass Flux (T=0) (kg/s)
# title 44 = Updraft Mass Flux (T=-10) (kg/s)
# title 45 = Updraft Mass Flux (T=-20) (kg/s)
# title 46 = Updraft Mass Flux (T=-30) (kg/s)
# title 47 = Crystal Mass Flux (T=-10)
# title 48 = Crystal Mass Flux (T=-20)
# title 49 = Crystal Mass Flux (T=-30)
# title 50 = icnetchg
# title 51 = Downdraft Volume (< -5m/s)
# title 52 = Updraft Volume (> 5m/s)
# title 53 = Updraft Volume (> 10m/s)
# title 54 = Updraft Volume (> 20m/s)
# title 55 = Charge Fallout (C)
# title 56 = Advection Charge Loss (C)
# title 57 = Ion Drift Charge (C)
# title 58 = Other Ion (C)
# title 59 = Unaccounted Charge (C)
# title 60 = Ion Flux: Top
# title 61 = Ion Flux: Bottom
# title 62 = Ion Flux: Sides
# title 63 = Positive Sources
# title 64 = Negative Sources
# title 65 = Total Sources
# title 66 = Mean Graupel Density
# title 67 = Mean Hail Density
# title 68 = Total Liquid Fraction
# title 69 = Sed/Mic charge (C)
# title 70 = Sed Charge2 (C)


    
#reader1 = csv.reader(f,delimiter=' ',)
#for i in range(numtitles):
#    print "title " + str(i) + ' = ' + titles[i]
