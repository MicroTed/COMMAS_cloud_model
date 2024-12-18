#!/opt/local/anaconda/bin/python

# Purpose: Plot time series (bin.ion.pl) from multiple simulations

import csv
import sys
import pandas as pd

import numpy as np
import numpy.ma as ma
import math
import matplotlib.pyplot as plt
from matplotlib import rc

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

runhead = 'may29zvh1000cn'
tail = 'vtrr'
runtail =  tail + '.bin.txt'
outputhead = runhead + 'XX' + tail 

vtlist = ['07','08','09','10','11','12','13']
label_list = ['0.7 Vt','0.8 Vt','0.9 Vt','1.0 Vt','1.1 Vt','1.2 Vt','1.3 Vt']
color_list = ['skyblue','dodgerblue','blue','black','red','lightcoral','violet']
#color_list = ['skyblue','dodgerblue','blue','black','crimson','red','violet']
files = []
numfiles = len(vtlist)

# overlay says whether to read in a second set to plot with dashed lines
overlay = True
#runhead2 = '../may29zvh1000cn6cdxrr/may29zvh1000cn' # icdx=6, variable density
#tail2 = 'vt6cdxrr'
#run2 = ''
runhead2 = '../may29zvh1000cnvtgrr/may29zvh1000cn' # only graupelfallfac is changed
tail2 = 'vtgrr'
run2 = ''
runhead2 = '../may29zh1000cnvtrr/may29zh1000cn' # fixed density
tail2 = 'vtrr'
run2 = 'zhXX'

runtail2 =  tail2 + '.bin.txt'
files2 = []
if (overlay):
    outputhead = outputhead + '.' + run2 + tail2


# set index
index = 26

# set output file name depending on index

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

#pp = PdfPages(runhead+'XXvtrr'+'.bin.hailmass.pdf')

plt.style.use('presentation')
# print plt.style.available

# create list of file names using vtlist
for num in vtlist:
   name = runhead+num+runtail
   files.append(name)
   if ( overlay ):
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
   dat = pd.read_csv(filename,header=3,sep=',')
   

   #print (dat.head(3))
   #print (dat.tail(3))
   
   alldat.append(dat)

if ( overlay ):
   for filename in files2:
      dat = pd.read_csv(filename,header=3,sep=',')
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
time_min = alldat[0].iloc[:,0]

for i in range(numfiles):
    graupelmass = alldat[i].iloc[:,index]
    plt.plot(time_min,graupelmass,label=label_list[i], linewidth=2, color=color_list[i])

if ( overlay ): 
    for i in range(numfiles):
        graupelmass = alldat2[i].iloc[:,index]
        plt.plot(time_min,graupelmass,'--',label=None, linewidth=2, color=color_list[i])

    
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
