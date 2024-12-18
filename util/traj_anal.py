#This script creates time-height plots of various quantities along trajectories

import numpy as N
import matplotlib
import matplotlib.cm as cm
from matplotlib.ticker import *
import matplotlib.colors as colors
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import ImageGrid,make_axes_locatable
#matplotlib.use('Agg')
#from pylab import *
import sys
sys.path.append('/Users/ddawson/python_scripts/')
import arpsmodule as arps
import thermolib as thermo
import trajmodule
import trajcmodule as trajc
import raymond_lowpass as raymond
import matplotlib.pyplot as plt
import scipy.integrate as integ
import gc
import os

def mtokm(val,pos):
    """Convert m to km for formatting axes tick labels"""
    val=val/1000.0
    return '%i' % val

def plotsingle(x,y,xlim,ylim,field,clevels,cmap,ovrmap,savefig,ovrfieldopt,ovrfield,ovrfieldlvl):
    """Plot a field on its own plot"""
    fig = plt.figure(figsize=(8,6))
    axes = fig.add_subplot(111)
    plot = axes.contourf(x,y,field,levels=clevels,cmap=cmap)
     # Find a nice number of ticks for the colorbar
    cintv = clevels[1]-clevels[0]
    cintvs = N.arange(clevels[0],clevels[-1],cintv)
    while True:
        if(cintvs.size > 20):
            cintv = (cintvs[1]-cintvs[0])*2.
            cintvs = N.arange(cintvs[0],cintvs[-1],cintv)
        else:
            break
    clvllocator = MultipleLocator(base=cintv)
    plt.colorbar(plot,orientation='vertical',ticks=clvllocator)
    if(ovrfieldopt):
        plotovr = axes.contour(x,y,ovrfield,levels=ovrfieldlvl,colors='k')
    axes.set_xlim(xlim[0],xlim[1])
    axes.set_ylim(ylim[0],ylim[1])
    formatter = plt.FuncFormatter(mtokm)
    axes.xaxis.set_major_formatter(formatter)
    axes.yaxis.set_major_formatter(formatter)
    axes.xaxis.set_major_locator(plt.MultipleLocator(base=1000.0))
    axes.yaxis.set_major_locator(plt.MultipleLocator(base=1000.0))
    axes.set_aspect('equal')

    if(ovrmap): # Overlay map
        bgmap.readshapefile(track_shapefile_location,'track',drawbounds=True,linewidth=0.5,color='black',ax=axes)
        bgmap.readshapefile(county_shapefile_location,'counties',drawbounds=True, linewidth=0.5, color='gray',ax=axes)  #Draws US county boundaries.

    return fig,axes


def integrate_trajc(initial_value,field,step):
    """Simple forward integration along a trajectory given a field at each time, a time step value, and an initial value
       for the quantity being integrated.  Example: if initial_value=initial height, field=vertical velocity,
       then the function returns the integrated height along the trajectory."""
    
    fieldinteg = N.zeros_like(field)
    fieldinteg[0] = initial_value
    
    for t in xrange(1,field.size):
        fieldinteg[t] = fieldinteg[t-1] + field[t-1]*step
        
    return fieldinteg
    
def integrate_trajcs(initial_values,field,step):
    """Same as integrate_trajc except it does it for multiple trajectories at once"""
    
    fieldinteg = N.zeros_like(field)
    fieldinteg[0,:] = initial_values[:]
    
    for t in xrange(1,N.size(field,0)):
        fieldinteg[t,:] = fieldinteg[t-1,:] + field[t-1,:]*step
        
    return fieldinteg

def integrate_trajc_trap(initial_value,field,step):
    """Integrates along a trajectory using the composite trapezoidal rule."""
    
    fieldinteg = integ.trapz(field,dx=step)+initial_value
    
    return fieldinteg

def integrate_trajcs_trap(initial_values,field,step):
    """Integrates along trajectories using the composite trapezoidal rule."""
    
    fieldinteg = integ.trapz(field,dx=step,axis=0)+initial_values
    
    return fieldinteg
        
Rd=287.0                #Gas constant for dry air (J/kg/K)
cp=1005.6               #Specific heat at constant pressure of dry air (J/kg/K)
Lv=2.501e6              #Latent heat of vaporization (J/kg)
Lf=3.34e5               #Latent heat of freezing (J/kg)
Ls=Lv+Lf                #Latent heat of sublimation (J/kg)

# Intensifying stage
# tstart_list = N.array([1830.])
# tstop_list = N.array([3030.])
# tref_list = N.array([2730.])

# Peak intensity
tstart_list = N.arange(3480.0,3510.0,30.0)
tstop_list = N.arange(4680.0,4710.0,30.0)
tref_list = N.arange(4380.0,4410.0,30.0)

# Weakening phase
# tstart_list = N.array([4200.])
# tstop_list = N.array([5400.])
# tref_list = N.array([5100.])

# Full time period
# tstart_list = N.arange(900.0,5700.0,30.0)      #Start time of trajectories
# tstop_list = N.arange(2100.0,6900.0,30.0)    #End time of trajectories
# tref_list = N.arange(1800.0,6600.0,30.0)      #Reference time of trajectories

# tstart_list = N.arange(3510.0,3930.0,30.0)
# tstop_list = N.arange(4710.0,5130.0,30.0)
# tref_list = N.arange(4410.0,4830.0,30.0)

# tref_list = N.arange(4260.0,4500.0,30.0)
# tstart_list = N.arange(4260.0-900.0,4500.0-900.0,30.0)
# tstop_list = N.arange(4260.0+300.0,4500.0+300.0,30.0)

tintv=30.0                                  #Interval of model data (possibly time interpolated)
trajintv=2.0                                #Interval of trajectories

runname='250m0503992305_1km3DVARCA225010min90_3km03001hr9MY3rv'
trailer='.520x560x053'
mphyopt = 11
dir='250m_Drobo/'+runname+'_subdomain/'
dir2='250m_scr/'+runname+'_subdomain/'
dir_extra='extra/'
dir_wmomentum='wmomentum_fields/'
dir_DSD='DSD_data/'
dir_images='trajanal_images_paper/'
dir_3dimages = 'traj3d_images/'
dir_npzs_read='trajsfccalc_npzs/'
dir_npzs='trajsfcanal_npzs/'
trajdir='/Users/ddawson/torcases/may0399/bugfix_realdata/arpstrajc_sfcanal_work_2ssub/'
#trajdir='/Volumes/scr/torcases/may0399/bugfix_realdata/arpstrajc_grid50magl_work_2ssub/'
#trajdir='/Users/ddawson/torcases/may0399/bugfix_realdata/arpstrajc_grid1kmagl_work_2ssub/'
#trajdir='/Users/ddawson/torcases/may0399/bugfix_realdata/arpstrajc_sfcanal_work_2ssub_test/'
trajxoffset = -54750.0
trajyoffset = -64750.0
varnames=['w','pte','vor','pprt_r','ptprtr','qvprt']
#varnames=['vor','pprt_r','ptprtr','qvprt']
#varnames=['w']
wforce_varnames=['bforce','ppgf_d','waccel','buoy__','buoy_t','buoy_p','buoy_q','w_ppgf','wbandp']
DSD_varnames=['qc','qr','qh','evapqc','evapqr','meltqh','Dmr___','Dmh___','alphar','alphah']
#DSD_varnames=[]
npzreaddir = dir2+dir_npzs_read
npzsavedir = dir2+dir_npzs

plot_var                = True  # Analyze and create plots for "normal" variables?
plot_wforce             = True  # w-forcing variables?
plot_DSD                = False  # DSD variables?

plot_minmaxavg          = True  # Create XY plots of trajectory quantities at reference time?
plot_TS                 = False # Create timeseries plots of trajectory quantities?
plot_hist               = False # Create histograms of trajectory quantities?
plot_scatter            = True # Create scatter plots of trajectory quantities vs. height in TLV?
plot3D                  = True # Make some 3D scatterplots of trajectories?
plot_torpoints          = False # Plot all points that end up in tornado on 3D plot?

# Arrays for height bins for scatterplots
height_bins = N.arange(0.0,5100.0,100.0)    # Every 100 m
height_middle = N.arange(50.0,5000.0,100.0) # Middle of each bin

w_lim3d = [-50.0,50.0]
pte_lim3d = [315.0,355.0]
vor_lim3d = [0.0,0.5]
pprt_r_lim3d = [-4000.0,4000.0]
ptprtr_lim3d = [-10.0,10.0]
qvprt_lim3d = [-0.005,0.005]
wforce_lim3d = [-2.5,2.5]
varlims3D = [w_lim3d,pte_lim3d,vor_lim3d,pprt_r_lim3d,ptprtr_lim3d,qvprt_lim3d]
qr_lim3d = [0.0,5.0e-3]
qc_lim3d = [0.0,5.0e-3]
qh_lim3d = [0.0,5.0e-3]
evapqr_lim3d = [0.0,0.01e-3]
evapqc_lim3d = [0.0,0.01e-3]
meltqh_lim3d = [0.0,0.01e-3]
Dmr_lim3d = [0.0,5.0]
Dmh_lim3d = [0.0,5.0]
alphar_lim3d = [0.0,15.0]
alphah_lim3d = [0.0,15.0]
DSDlims3D = [qc_lim3d,qr_lim3d,qh_lim3d,evapqc_lim3d,evapqr_lim3d,meltqh_lim3d,
                Dmr_lim3d,Dmh_lim3d,alphar_lim3d,alphah_lim3d]


plot_integrated_wforcing = False    # It appears this is pretty inaccurate EDIT: for 15 min period, not bad for 5 min
plot_integrated_theta    = False    # Ditto

saveTS                   = False    # Save timeseries?
loadTS                   = True    # Load timeseries?

masked_tornado = True         # Set to true to mask trajectories at each reference time
                                    # by those that enter the tornado in the subsequent 5 min
                                    # according to the vorticity criteria given below
maxvor_thresh = 0.1                 # Threshold of max vorticity to determine tornado trajectory
maxw_thresh = 20.0                  # Threshold of max vertical velocity to determine tornado trajectory (not used currently)
tormask = None

usetorreftime = True            # If True, also compute "before" and "after" periods relative to 
                                # when trajectory enters tornado

# Construct some needed directories to save plots
imagesavedir = dir2+dir_images
if (not os.path.exists(imagesavedir)):
    os.mkdir(imagesavedir)

image3ddir = dir2+dir_3dimages
if (not os.path.exists(imagesavedir)):
    os.mkdir(image3ddir)

# Start trajectory file loop

for t,tstart,tstop,tref in zip(xrange(tref_list.size),tstart_list,tstop_list,tref_list):
    
    if(loadTS and not (plot_integrated_wforcing or plot_minmaxavg or plot_integrated_theta)):
        break
    
    timestring = "%06d" % tstart
    timestringref = "%06d" % tref
    hisfilename=dir+runname+'.hdf'+timestring+trailer
    hisfilenameref=dir+runname+'.hdf'+timestringref+trailer
    
    if(tstart == tstart_list[0]):
        #Read the grid arrays from the first arps history file
        print "Reading grid and base state data from file "+hisfilename
        nx,ny,nz,dx,dy,x,y,zp,xs,ys,zs,zpagl,zsagl = arps.readarpsgrid(hisfilename)
    
    ntimes=int((tstop-tstart)/tintv)+1
    timerange=N.arange(tstart,tstop+tintv,tintv)
    trajtimestring = '.trajc_'+'%06d' % int(tstart)+'-'+'%06d' % int(tstop)+'_'+'%06d' % int(tref)
    
    # Construct directories to save images

    if(tstart == tstart_list[0]):
        imagetrajdir = imagesavedir+runname+trajtimestring
        if (not os.path.exists(imagetrajdir)):
            os.mkdir(imagetrajdir)

    # Read in the trajectory file
    
    trajfilename = trajdir+runname+trajtimestring
    numtrajcs,npoints,starttime,endtime,ttrajc,trajdata = trajc.readtrajc2(trajfilename,tintv,trajintv)
    print "Number of trajectories in trajectory file is "+str(numtrajcs)
    
    ivals = N.empty((npoints,numtrajcs))
    jvals = N.empty((npoints,numtrajcs))
    #kvals = empty((npoints,numtrajcs),dtype=int)
    
    xvals=N.empty((npoints,numtrajcs))
    yvals=N.empty((npoints,numtrajcs))
    zvals=N.empty((npoints,numtrajcs))
    zvals_agl=N.empty((npoints,numtrajcs))
    size=N.empty((npoints,numtrajcs))
    ztrajc_max = 0
    
    for traj in range(numtrajcs):
        curtrajdata=trajdata[:,:,traj]
        xtrajc = curtrajdata[:,0]+trajxoffset
        ytrajc = curtrajdata[:,1]+trajyoffset
        ztrajc = curtrajdata[:,2]
        
        # Save index of trajectory with greatest initial height (but interpolation will be done for
        # all trajectories)
        if (ztrajc[0] > ztrajc_max):
            traj_maxz = traj
            ztrajc_max = ztrajc[0]
        
        for i in range(npoints):        # Note npoints should be same as ntimes for right now
            #itrajc = max(0,min(nx-1,around((xtrajc[i]-xs[0])/dx)))
            #jtrajc = max(0,min(ny-1,around((ytrajc[i]-ys[0])/dx)))
                
            itrajc = (xtrajc[i]-xs[0])/dx
            jtrajc = (ytrajc[i]-ys[0])/dy
            
            ivals[i,traj] = itrajc
            jvals[i,traj] = jtrajc
            xvals[i,traj] = xtrajc[i]
            yvals[i,traj] = ytrajc[i]
            zvals[i,traj] = ztrajc[i]
            zvals_agl[i,traj] = ztrajc[i]-zp[min(itrajc,nx-1),min(jtrajc,ny-1),1]
            
            #Find the fractional k-level of the trajectory height
            #flag=0
            #for k in range(1,kend):
                #print 'ztrajc,zp',ztrajc[i],zp[itrajc,jtrajc,k]
            #    if(zp[itrajc,jtrajc,k] >= ztrajc[i]):
            #        flag=1
            #       if(abs(ztrajc[i]-zp[itrajc,jtrajc,k]) < 
            #        abs(ztrajc[i]-zp[itrajc,jtrajc,k-1])):
            #            ktrajc = k
            #        else:
            #            ktrajc = k-1
            #        break
            #    if(flag == 0):
            #        ktrajc = k
            #kvals[i,traj] = ktrajc
    
    print 'The trajectory number with the maximum initial height is '+str(traj_maxz)
    print 'That height is '+str(zvals_agl[0,traj_maxz])+' m AGL'
    
    #Load interpolated fields from npz file
        
    npzfilename = npzreaddir+'/'+runname+'_intrpfields_'+trajtimestring+'.npz'
    print "Loading interpolated fields from file "+npzfilename
    trajcalc_file=N.load(npzfilename)
    
    trefindex = [i for i,time in enumerate(timerange) if time == tref][0]   #Index of reference time   
    treftimestring = "%06d" % tref   
    
    # Initialize arrays to hold interpolated quantities
    varintrp = N.zeros((len(varnames),ntimes,numtrajcs))
    for i,var in enumerate(varnames):
        varintrp[i,:] = trajcalc_file[var]
    
    varbtref = varintrp[:,:trefindex+1,:]
    varatref = varintrp[:,trefindex:,:]
    
    zmaxbtref = zvals[:trefindex+1,:].max(axis=0)
    zmaxatref = zvals[trefindex:,:].max(axis=0)
    
    varminbtref = varbtref.min(axis=1)
    varminatref = varatref.min(axis=1)

    varmaxbtref = varbtref.max(axis=1)
    varmaxatref = varatref.max(axis=1)

    varavgbtref = varbtref.mean(axis=1)
    varavgatref = varatref.mean(axis=1)    
    
    if(plot_wforce):
        wforceintrp = N.zeros((len(wforce_varnames),ntimes,numtrajcs))
        for i,var in enumerate(wforce_varnames):
            wforceintrp[i,:] = trajcalc_file[var]
            
        wforcebtref = wforceintrp[:,:trefindex+1,:]
        wforceatref = wforceintrp[:,trefindex:,:]

        wforceminbtref = wforcebtref.min(axis=1)
        wforceminatref = wforceatref.min(axis=1)

        wforcemaxbtref = wforcebtref.max(axis=1)
        wforcemaxatref = wforceatref.max(axis=1)

        wforceavgbtref = wforcebtref.mean(axis=1)
        wforceavgatref = wforceatref.mean(axis=1)
        
    if(plot_DSD or plot_var):
        DSDintrp = N.zeros((len(DSD_varnames),ntimes,numtrajcs))
        for i,var in enumerate(DSD_varnames):
            DSDintrp[i,:] = trajcalc_file[var]
        
        # Set up some masks
        qrintrp = DSDintrp[DSD_varnames.index('qr'),:]*1000.0 # Multiply by 1000 to convert to g/kg
        qrmask = N.where(qrintrp > 0.0,False,True)
        qcintrp = DSDintrp[DSD_varnames.index('qc'),:]*1000.0
        qcmask = N.where(qcintrp > 0.0,False,True)
        qhintrp = DSDintrp[DSD_varnames.index('qh'),:]*1000.0
        qhmask = N.where(qhintrp > 0.0,False,True)
    
        # Mask arrays by nonzero hydrometeor content
        print "Masking arrays by nonzero hydrometeor content"
    
#       evapqcintrp = N.ma.masked_array(DSDintrp[DSD_varnames.index('evapqc'),:],mask=qcmask)*1000.0   # Convert to g/kg/s
#       evapqrintrp = N.ma.masked_array(DSDintrp[DSD_varnames.index('evapqr'),:],mask=qrmask)*1000.0
#       meltqhintrp = N.ma.masked_array(DSDintrp[DSD_varnames.index('meltqh'),:],mask=qhmask)*1000.0
    
        evapqcintrp = DSDintrp[DSD_varnames.index('evapqc'),:]*1000.0   # Convert to g/kg/s
        evapqrintrp = DSDintrp[DSD_varnames.index('evapqr'),:]*1000.0
        meltqhintrp = DSDintrp[DSD_varnames.index('meltqh'),:]*1000.0
    
        Dmrintrp = N.ma.masked_array(DSDintrp[DSD_varnames.index('Dmr___'),:],mask=qrmask)
        #print "Dmrintrp: "+str(Dmrintrp[:,traj_maxz])
        Dmhintrp = N.ma.masked_array(DSDintrp[DSD_varnames.index('Dmh___'),:],mask=qhmask)
        alpharintrp = N.ma.masked_array(DSDintrp[DSD_varnames.index('alphar'),:],mask=qrmask)
        alphahintrp = N.ma.masked_array(DSDintrp[DSD_varnames.index('alphah'),:],mask=qhmask)
        
        qrbtref = qrintrp[:trefindex+1,:]
        qcbtref = qcintrp[:trefindex+1,:]
        qhbtref = qhintrp[:trefindex+1,:]
        evapqcbtref = evapqcintrp[:trefindex+1,:]
        evapqrbtref = evapqrintrp[:trefindex+1,:]
        meltqhbtref = meltqhintrp[:trefindex+1,:]
        Dmrbtref = Dmrintrp[:trefindex+1,:]
        Dmhbtref = Dmhintrp[:trefindex+1,:]
        alpharbtref = alpharintrp[:trefindex+1,:]
        alphahbtref = alphahintrp[:trefindex+1,:]
        
        qrmaxbtref = qrbtref.max(axis=0)
        qcmaxbtref = qcbtref.max(axis=0)
        qhmaxbtref = qhbtref.max(axis=0)
        evapqcmaxbtref = evapqcbtref.max(axis=0)
        evapqrmaxbtref = evapqrbtref.max(axis=0)
        meltqhmaxbtref = meltqhbtref.max(axis=0)
        Dmrmaxbtref = Dmrbtref.max(axis=0)
        Dmhmaxbtref = Dmhbtref.max(axis=0)
        alpharmaxbtref = alpharbtref.max(axis=0)
        alphahmaxbtref = alphahbtref.max(axis=0)
    
        qravgbtref = qrbtref.mean(axis=0)
        qcavgbtref = qcbtref.mean(axis=0)
        qhavgbtref = qhbtref.mean(axis=0)
        evapqcavgbtref = evapqcbtref.mean(axis=0)
        evapqravgbtref = evapqrbtref.mean(axis=0)
        meltqhavgbtref = meltqhbtref.mean(axis=0)
        Dmravgbtref = Dmrbtref.mean(axis=0)
        Dmhavgbtref = Dmhbtref.mean(axis=0)
        alpharavgbtref = alpharbtref.mean(axis=0)
        alphahavgbtref = alphahbtref.mean(axis=0)
    
        Dmrminbtref = Dmrbtref.min(axis=0)
        Dmhminbtref = Dmhbtref.min(axis=0)
        alpharminbtref = alpharbtref.min(axis=0)
        alphahminbtref = alphahbtref.min(axis=0)
        
    if(masked_tornado):    # Create masks for trajectories that enter tornado and for points along 
                           # trajectories that are inside the tornado.
        tormask = N.where(varmaxatref[varnames.index('vor'),:] > maxvor_thresh,False,True)
        xvalstor = N.ma.masked_array(xvals[trefindex,:],mask=tormask)
        yvalstor = N.ma.masked_array(yvals[trefindex,:],mask=tormask)
        numtortrajcs = N.ma.count(xvalstor)
        print "There are "+str(numtortrajcs)+" trajectories that enter the tornado in the next 5 min."
        
        if(usetorreftime): # Masks points outside TLV (TLV trajectory mask applied later)
            torrefmask = N.where(varintrp[varnames.index('vor'),trefindex:,:] > 0.025,False,True)
            nontorrefmask = ~torrefmask
            torrefmaskfull = N.where(varintrp[varnames.index('vor'),...] > 0.025,False,True)
            nontorrefmaskfull = ~torrefmaskfull
                        
            indexattor = N.argmin(torrefmask,axis=0)    # argmin on boolean array returns index of first "False" 
                                                        # along given axis! Neat, huh?
                                                        # But if all True, then returns first element,
                                                        # so have to be careful
        
            xvalstorref = xvals[indexattor+trefindex,N.arange(numtrajcs)]
            xvalstorref = xvalstorref[~tormask]
        
            yvalstorref = yvals[indexattor+trefindex,N.arange(numtrajcs)]
            yvalstorref = yvalstorref[~tormask]
            
            # Use index at TLV to remake mask for points within and without TLV
            #torrefmaskfull = N.zeros_like(torrefmaskfull,dtype=N.bool)
            torrefmaskfull = N.array([[i < (indexattor[j]+trefindex-1) for j in xrange(numtrajcs)] for i in xrange(npoints)])
            nontorrefmaskfull = ~torrefmaskfull
            
            # Find trajectories that enter TLV for entire 5-min period (those closest to location of TLV at the reference time).
            tor5minmask = N.where((indexattor == 0) & ~tormask,False,True)
            xvalstor5min = xvals[trefindex,~tor5minmask]
            yvalstor5min = yvals[trefindex,~tor5minmask]
            
            
            if(plot_var):
                dims = torrefmask.shape
                torrefmask_var = torrefmask.reshape(1,dims[0],dims[1]).repeat(len(varnames),0)  # repeats mask for every variable
                tormask_var = tormask.reshape(1,numtrajcs).repeat(len(varnames),0)
                
                dims = torrefmaskfull.shape
                torrefmaskfull_var = torrefmaskfull.reshape(1,dims[0],dims[1]).repeat(len(varnames),0)  # repeats mask for every variable
                nontorrefmaskfull_var = nontorrefmaskfull.reshape(1,dims[0],dims[1]).repeat(len(varnames),0)
                
                varbtortref = N.ma.masked_array(varintrp,mask=nontorrefmaskfull_var)
                varatortref = N.ma.masked_array(varintrp,mask=torrefmaskfull_var)
            
                #Min/max/avg before and after entering tornado
                varminbtortref = varbtortref.min(axis=1)
                varminatortref = varatortref.min(axis=1)
                varmaxbtortref = varbtortref.max(axis=1)
                varmaxatortref = varatortref.max(axis=1)
                varavgbtortref = varbtortref.mean(axis=1)
                varavgatortref = varatortref.mean(axis=1)
                
                valattor = varatref[:,indexattor,N.arange(numtrajcs)]       # Grabs value of each variable at the tornado
                                                        # reference time for each trajectory
                
#                 
#                 # Finally mask out trajectories that don't enter tornado (have to do this separately
#                 # because of some strange quirk that I haven't figured out yet).
#                 for i,var in enumerate(varnames):
#                     valattor[i,:].mask=tormask
#                     varminbtortref[i,:].mask=tormask
#                     varmaxbtortref[i,:].mask=tormask
#                     varavgbtortref[i,:].mask=tormask
#                     varminatortref[i,:].mask=tor5minmask   # Change to those within TLV for entire 5 min
#                     varmaxatortref[i,:].mask=tor5minmask
#                     varavgatortref[i,:].mask=tor5minmask
            
            
            if(plot_wforce):
                dims = torrefmask.shape
                torrefmask_wforce = torrefmask.reshape(1,dims[0],dims[1]).repeat(len(wforce_varnames),0)  # repeats mask for every variable
                tormask_wforce = tormask.reshape(1,numtrajcs).repeat(len(wforce_varnames),0)
            
                dims = torrefmaskfull.shape
                torrefmaskfull_wforce = torrefmaskfull.reshape(1,dims[0],dims[1]).repeat(len(wforce_varnames),0)  # repeats mask for every variable
                nontorrefmaskfull_wforce = ~torrefmaskfull_wforce
                
                wforcebtortref = N.ma.masked_array(wforceintrp,mask=~torrefmaskfull_wforce)
                wforceatortref = N.ma.masked_array(wforceintrp,mask=torrefmaskfull_wforce)
            
                #Min/max/avg before and after entering tornado
                wforceminbtortref = wforcebtortref.min(axis=1)
                wforceminatortref = wforceatortref.min(axis=1)
                wforcemaxbtortref = wforcebtortref.max(axis=1)
                wforcemaxatortref = wforceatortref.max(axis=1)
                wforceavgbtortref = wforcebtortref.mean(axis=1)
                wforceavgatortref = wforceatortref.mean(axis=1)

                wforceattor = wforceatref[:,indexattor,N.arange(numtrajcs)]       # Grabs value of each variable at the tornado
                                                            # reference time for each trajectory
#                 
#                 # Finally mask out trajectories that don't enter tornado (have to do this separately
#                 # because of some strange quirk that I haven't figured out yet).
#                 for i,var in enumerate(wforce_varnames):
#                     wforceattor[i,:] = N.ma.masked_array(wforceattor[i,:],mask=tormask)
#                     wforceminbtortref[i,:] = N.ma.masked_array(wforceminbtortref[i,:],mask=tormask)
#                     wforcemaxbtortref[i,:] = N.ma.masked_array(wforcemaxbtortref[i,:],mask=tormask)
#                     wforceavgbtortref[i,:] = N.ma.masked_array(wforceavgbtortref[i,:],mask=tormask)
#                     wforceminatortref[i,:] = N.ma.masked_array(wforceminatortref[i,:],mask=tor5minmask)
#                     wforcemaxatortref[i,:] = N.ma.masked_array(wforcemaxatortref[i,:],mask=tor5minmask)
#                     wforceavgatortref[i,:] = N.ma.masked_array(wforceavgatortref[i,:],mask=tor5minmask)
            
            if(plot_DSD):
                qrbtortref = N.ma.masked_array(qrintrp,mask=~torrefmaskfull)
                qcbtortref = N.ma.masked_array(qcintrp,mask=~torrefmaskfull)
                qhbtortref = N.ma.masked_array(qhintrp,mask=~torrefmaskfull)
                evapqcbtortref = N.ma.masked_array(evapqcintrp,mask=~torrefmaskfull)
                evapqrbtortref = N.ma.masked_array(evapqrintrp,mask=~torrefmaskfull)
                meltqhbtortref = N.ma.masked_array(meltqhintrp,mask=~torrefmaskfull)
                Dmrbtortref = N.ma.masked_array(Dmrintrp,mask=~torrefmaskfull)
                Dmhbtortref = N.ma.masked_array(Dmhintrp,mask=~torrefmaskfull)
                alpharbtortref = N.ma.masked_array(alpharintrp,mask=~torrefmaskfull)
                alphahbtortref = N.ma.masked_array(alphahintrp,mask=~torrefmaskfull)
        
                #Min/max/avg before and after entering tornado
                qrmaxbtortref = qrbtortref.max(axis=0)
                qcmaxbtortref = qcbtortref.max(axis=0)
                qhmaxbtortref = qhbtortref.max(axis=0)
                evapqrmaxbtortref = evapqrbtortref.max(axis=0)
                evapqcmaxbtortref = evapqcbtortref.max(axis=0)
                meltqhmaxbtortref = meltqhbtortref.max(axis=0)
                Dmrmaxbtortref = Dmrbtortref.max(axis=0)
                Dmhmaxbtortref = Dmhbtortref.max(axis=0)
                alpharmaxbtortref = alpharbtortref.max(axis=0)
                alphahmaxbtortref = alphahbtortref.max(axis=0)
            
                qravgbtortref = qrbtortref.mean(axis=0)
                qcavgbtortref = qcbtortref.mean(axis=0)
                qhavgbtortref = qhbtortref.mean(axis=0)
                evapqravgbtortref = evapqrbtortref.mean(axis=0)
                evapqcavgbtortref = evapqcbtortref.mean(axis=0)
                meltqhavgbtortref = meltqhbtortref.mean(axis=0)
                Dmravgbtortref = Dmrbtortref.mean(axis=0)
                Dmhavgbtortref = Dmhbtortref.mean(axis=0)
                alpharavgbtortref = alpharbtortref.mean(axis=0)
                alphahavgbtortref = alphahbtortref.mean(axis=0)
            
                Dmrminbtortref = Dmrbtortref.min(axis=0)
                Dmhminbtortref = Dmhbtortref.min(axis=0)
                alpharminbtortref = alpharbtortref.min(axis=0)
                alphahminbtortref = alphahbtortref.min(axis=0)
            
    # Create empty arrays for scatterplot timeseries and initialize height bins
    if(plot_scatter and tstart == tstart_list[0]):
        if(plot_var):
            npoints_TS = N.empty((tstart_list.size,height_middle.size))
            npointsb_TS = N.empty_like(npoints_TS)
            var_avg_TS = N.empty((len(varnames),tstart_list.size,height_middle.size))
            varb_avg_TS = N.empty_like(var_avg_TS)
        if(plot_wforce):
            wforce_avg_TS = N.empty((len(wforce_varnames),tstart_list.size,height_middle.size))
            wforceb_avg_TS = N.empty_like(wforce_avg_TS)
        if(plot_DSD):
            DSD_avg_TS = N.empty((len(DSD_varnames),tstart_list.size,height_middle.size))
            DSDb_avg_TS = N.empty_like(DSD_avg_TS)
    
    
    # Set up plot bounds, etc.
    # Determine plot boundaries
    xmin = xvals[trefindex,:].min()     #These values should lie on scalar points at the reference time
    xmax = xvals[trefindex,:].max()
    ymin = yvals[trefindex,:].min()
    ymax = yvals[trefindex,:].max()
    
    # Determine plot boundary indices
    try:
        imin = N.where(xs == xmin)[0][0]
    except:
        imin = 1
    try:
        imax = N.where(xs == xmax)[0][0]
    except:
        imax = nx-2
    try:
        jmin = N.where(ys == ymin)[0][0]
    except:
        jmin = 1
    try:
        jmax = N.where(ys == ymax)[0][0]
    except:
        jmax = ny-2
    
    #print "Test: Plot east bound index: "+str(imax-imin+1)
    #print "Test: Plot north bound index: "+str(jmax-jmin+1)
 
    xsplt = xs[imin:imax+1]
    ysplt = ys[jmin:jmax+1]
    xlim = (xsplt[0],xsplt[-1])
    ylim = (ysplt[0],ysplt[-1])    
            
    if(plot_minmaxavg):
        # Reshape arrays for plotting

        if(plot_var):
            zmaxbtref2D = zmaxbtref.reshape((int(N.sqrt(numtrajcs)),-1))[:jmax-jmin+1,:imax-imin+1]
            zmaxatref2D = zmaxatref.reshape((int(N.sqrt(numtrajcs)),-1))[:jmax-jmin+1,:imax-imin+1]

            varminbtref2D = varminbtref.reshape((len(varnames),int(N.sqrt(numtrajcs)),-1))[:,:jmax-jmin+1,:imax-imin+1]
            varminatref2D = varminatref.reshape((len(varnames),int(N.sqrt(numtrajcs)),-1))[:,:jmax-jmin+1,:imax-imin+1]

            varmaxbtref2D = varmaxbtref.reshape((len(varnames),int(N.sqrt(numtrajcs)),-1))[:,:jmax-jmin+1,:imax-imin+1]
            varmaxatref2D = varmaxatref.reshape((len(varnames),int(N.sqrt(numtrajcs)),-1))[:,:jmax-jmin+1,:imax-imin+1]

            varavgbtref2D = varavgbtref.reshape((len(varnames),int(N.sqrt(numtrajcs)),-1))[:,:jmax-jmin+1,:imax-imin+1]
            varavgatref2D = varavgatref.reshape((len(varnames),int(N.sqrt(numtrajcs)),-1))[:,:jmax-jmin+1,:imax-imin+1]

            # Now for the plotting!
        
            clevels = N.arange(0.0,4100.0,100.0)
            cmap = cm.Reds
        
            # Reshape interpolated trajectory array at reference time into 2D (assumes trajectory reference domain is square)
        
            qrintrp2D = DSDintrp[DSD_varnames.index('qr'),trefindex,:].reshape((int(N.sqrt(numtrajcs)),-1))
        
            fig,ax = plotsingle(xsplt,ysplt,xlim,ylim,zmaxbtref2D,clevels,cmap,False,False,True,qrintrp2D*1000.,[0.01,0.1,0.5,1.0,5.0])
            #Overlay wind vectors
            u = arps.hdfreadvar3d(hisfilenameref,'u')
            v = arps.hdfreadvar3d(hisfilenameref,'v')
            us = N.empty_like(u)
            vs = N.empty_like(v)
            us[:-1,:,:] = 0.5*(u[:-1,:,:]+u[1:,:,:])
            vs[:,:-1,:] = 0.5*(v[:,:-1,:]+v[:,1:,:])
            usplt = us.swapaxes(0,1)[jmin:jmax+1,imin:imax+1,1]
            vsplt = vs.swapaxes(0,1)[jmin:jmax+1,imin:imax+1,1]
            Q = plt.quiver(xsplt,ysplt,usplt,vsplt,pivot='middle', minshaft=1, scale=150, units='inches',headwidth=4, width=.01,color='k')
        
            if(masked_tornado and numtortrajcs > 0):
                #ax.scatter(xvalstor,yvalstor,facecolors='none',edgecolors='g',marker='o')
                if(usetorreftime):
                    ax.scatter(xvalstorref,yvalstorref,facecolors='none',edgecolors='y',marker='s')
        
            figfile = runname+'_'+treftimestring+'_height_max_before.eps'
            plt.savefig(imagesavedir+'/'+figfile,dpi=200)
            plt.close()
        
            for i,var in enumerate(varnames):
                print "Creating plots for variable "+var+" and reference time "+str(tref)
                if(var == 'w'):
                    clevels = N.arange(-40.0,40.5,0.5)
                    cmap = cm.RdBu_r
                elif(var == 'pt'):
                    clevels = N.arange(250.0,350.0,1.0)
                    cmap = cm.RdBu_r
                elif(var == 'qv'):
                    clevels = N.arange(0.0,16.2,0.2)
                    cmap = cm.Blues
                elif(var == 'pte'):
                    clevels = N.arange(300.0,400.0,1.0)
                    cmap = cm.RdBu_r
                elif(var == 'vor'):
                    clevels = N.arange(-0.5,0.51,0.01)
                    cmap = cm.RdBu_r
                elif(var == 'pprt_r'):
                    clevels = N.arange(-5000.0,5100.0,100.0)
                    cmap = cm.RdBu_r
                elif(var == 'ptprtr'):
                    clevels = N.arange(-20.0,20.5,0.5)
                    cmap = cm.RdBu_r
                elif(var == 'qvprt'):
                    clevels = N.arange(-0.01,0.011,0.001)
                    cmap = cm.RdBu_r

                fig,ax = plotsingle(xsplt,ysplt,xlim,ylim,varminbtref2D[i,:],clevels,cmap,False,False,True,qrintrp2D*1000.,[0.01,0.1,0.5,1.0,5.0])
                Q = plt.quiver(xsplt,ysplt,usplt,vsplt,pivot='middle', minshaft=1, scale=150, units='inches',headwidth=4, width=.01,color='k')
                figfile = runname+'_'+treftimestring+'_'+var+'_min_before.eps'
                if(masked_tornado and numtortrajcs > 0):
                    ax.scatter(xvalstor,yvalstor,facecolors='none',edgecolors='g',marker='o')
                plt.savefig(imagesavedir+'/'+figfile,dpi=200)
                plt.close()

                fig,ax = plotsingle(xsplt,ysplt,xlim,ylim,varmaxbtref2D[i,:],clevels,cmap,False,False,True,qrintrp2D*1000.,[0.01,0.1,0.5,1.0,5.0])
                Q = plt.quiver(xsplt,ysplt,usplt,vsplt,pivot='middle', minshaft=1, scale=150, units='inches',headwidth=4, width=.01,color='k')
                if(masked_tornado and numtortrajcs > 0):
                    ax.scatter(xvalstor,yvalstor,facecolors='none',edgecolors='g',marker='o')
                figfile = runname+'_'+treftimestring+'_'+var+'_max_before.eps'
                if(var == 'vor' and masked_tornado and numtortrajcs > 0 and usetorreftime):
                    ax.scatter(xvalstor5min,yvalstor5min,facecolors='none',edgecolors='k',marker='o')
                plt.savefig(imagesavedir+'/'+figfile,dpi=200)
                plt.close()

                fig,ax = plotsingle(xsplt,ysplt,xlim,ylim,varavgbtref2D[i,:],clevels,cmap,False,False,True,qrintrp2D*1000.,[0.01,0.1,0.5,1.0,5.0])
                Q = plt.quiver(xsplt,ysplt,usplt,vsplt,pivot='middle', minshaft=1, scale=150, units='inches',headwidth=4, width=.01,color='k')
                if(masked_tornado and numtortrajcs > 0):
                    ax.scatter(xvalstor,yvalstor,facecolors='none',edgecolors='g',marker='o')
            
                figfile = runname+'_'+treftimestring+'_'+var+'_avg_before.eps'
                plt.savefig(imagesavedir+'/'+figfile,dpi=200)
                plt.close()
        
        if(plot_DSD):
            print "Creating plots for DSD parameters for reference time "+str(tref)
        
            qrmaxbtref2D = qrmaxbtref.reshape((int(N.sqrt(numtrajcs)),-1))[:jmax-jmin+1,:imax-imin+1]
            qcmaxbtref2D = qcmaxbtref.reshape((int(N.sqrt(numtrajcs)),-1))[:jmax-jmin+1,:imax-imin+1]
            qhmaxbtref2D = qhmaxbtref.reshape((int(N.sqrt(numtrajcs)),-1))[:jmax-jmin+1,:imax-imin+1]
            evapqcmaxbtref2D = evapqcmaxbtref.reshape((int(N.sqrt(numtrajcs)),-1))[:jmax-jmin+1,:imax-imin+1]
            evapqrmaxbtref2D = evapqrmaxbtref.reshape((int(N.sqrt(numtrajcs)),-1))[:jmax-jmin+1,:imax-imin+1]
            meltqhmaxbtref2D = meltqhmaxbtref.reshape((int(N.sqrt(numtrajcs)),-1))[:jmax-jmin+1,:imax-imin+1]
            Dmrmaxbtref2D = Dmrmaxbtref.reshape((int(N.sqrt(numtrajcs)),-1))[:jmax-jmin+1,:imax-imin+1]
            Dmhmaxbtref2D = Dmhmaxbtref.reshape((int(N.sqrt(numtrajcs)),-1))[:jmax-jmin+1,:imax-imin+1]
            alpharmaxbtref2D = alpharmaxbtref.reshape((int(N.sqrt(numtrajcs)),-1))[:jmax-jmin+1,:imax-imin+1]
            alphahmaxbtref2D = alphahmaxbtref.reshape((int(N.sqrt(numtrajcs)),-1))[:jmax-jmin+1,:imax-imin+1]
        
            qravgbtref2D = qravgbtref.reshape((int(N.sqrt(numtrajcs)),-1))[:jmax-jmin+1,:imax-imin+1]
            qcavgbtref2D = qcavgbtref.reshape((int(N.sqrt(numtrajcs)),-1))[:jmax-jmin+1,:imax-imin+1]
            qhavgbtref2D = qhavgbtref.reshape((int(N.sqrt(numtrajcs)),-1))[:jmax-jmin+1,:imax-imin+1]
            evapqcavgbtref2D = evapqcavgbtref.reshape((int(N.sqrt(numtrajcs)),-1))[:jmax-jmin+1,:imax-imin+1]
            evapqravgbtref2D = evapqravgbtref.reshape((int(N.sqrt(numtrajcs)),-1))[:jmax-jmin+1,:imax-imin+1]
            meltqhavgbtref2D = meltqhavgbtref.reshape((int(N.sqrt(numtrajcs)),-1))[:jmax-jmin+1,:imax-imin+1]
            Dmravgbtref2D = Dmravgbtref.reshape((int(N.sqrt(numtrajcs)),-1))[:jmax-jmin+1,:imax-imin+1]
            Dmhavgbtref2D = Dmhavgbtref.reshape((int(N.sqrt(numtrajcs)),-1))[:jmax-jmin+1,:imax-imin+1]
            alpharavgbtref2D = alpharavgbtref.reshape((int(N.sqrt(numtrajcs)),-1))[:jmax-jmin+1,:imax-imin+1]
            alphahavgbtref2D = alphahavgbtref.reshape((int(N.sqrt(numtrajcs)),-1))[:jmax-jmin+1,:imax-imin+1]
        
            Dmrminbtref2D = Dmrminbtref.reshape((int(N.sqrt(numtrajcs)),-1))[:jmax-jmin+1,:imax-imin+1]
            Dmhminbtref2D = Dmhminbtref.reshape((int(N.sqrt(numtrajcs)),-1))[:jmax-jmin+1,:imax-imin+1]
            alpharminbtref2D = alpharminbtref.reshape((int(N.sqrt(numtrajcs)),-1))[:jmax-jmin+1,:imax-imin+1]
            alphahminbtref2D = alphahminbtref.reshape((int(N.sqrt(numtrajcs)),-1))[:jmax-jmin+1,:imax-imin+1]
        
            # Now for the plotting!
        
            clevels = N.arange(0.0,15.1,0.1)
            cmap = cm.Blues
        
            # Maximum qc,qr,qh prior to reference time (g/kg)
            fig,ax = plotsingle(xsplt,ysplt,xlim,ylim,qrmaxbtref2D,clevels,cmap,False,False,False,None,None)
            figfile=runname+'_'+treftimestring+'_qr_max_before.eps'
            plt.savefig(imagesavedir+'/'+figfile,dpi=200)
            plt.close()
        
            fig,ax = plotsingle(xsplt,ysplt,xlim,ylim,qcmaxbtref2D,clevels,cmap,False,False,False,None,None)
            figfile=runname+'_'+treftimestring+'_qc_max_before.eps'
            plt.savefig(imagesavedir+'/'+figfile,dpi=200)
            plt.close()
        
            fig,ax = plotsingle(xsplt,ysplt,xlim,ylim,qhmaxbtref2D,clevels,cmap,False,False,False,None,None)
            figfile=runname+'_'+treftimestring+'_qh_max_before.eps'
            plt.savefig(imagesavedir+'/'+figfile,dpi=200)
            plt.close()
        
            clevels = N.arange(0.0,0.011,0.001)
            cmap = cm.Blues
        
            # Maximum evapqc,evapqr,meltqh prior to reference time (g/kg/s)
            fig,ax = plotsingle(xsplt,ysplt,xlim,ylim,evapqrmaxbtref2D,clevels,cmap,False,False,False,None,None)
            figfile=runname+'_'+treftimestring+'_evapqr_max_before.eps'
            plt.savefig(imagesavedir+'/'+figfile,dpi=200)
            plt.close()
        
            fig,ax = plotsingle(xsplt,ysplt,xlim,ylim,evapqcmaxbtref2D,clevels,cmap,False,False,False,None,None)
            figfile=runname+'_'+treftimestring+'_evapqc_max_before.eps'
            plt.savefig(imagesavedir+'/'+figfile,dpi=200)
            plt.close()
        
            fig,ax = plotsingle(xsplt,ysplt,xlim,ylim,meltqhmaxbtref2D,clevels,cmap,False,False,False,None,None)
            figfile=runname+'_'+treftimestring+'_meltqh_max_before.eps'
            plt.savefig(imagesavedir+'/'+figfile,dpi=200)
            plt.close()
        
            clevels = N.arange(0.0,10.0,0.1)
            cmap = cm.Blues
        
            # Maximum Dmr, Dmh, alphar, alphah prior to reference time (mm, unitless)
            fig,ax = plotsingle(xsplt,ysplt,xlim,ylim,Dmrmaxbtref2D,clevels,cmap,False,False,False,None,None)
            figfile=runname+'_'+treftimestring+'_Dmr_max_before.eps'
            plt.savefig(imagesavedir+'/'+figfile,dpi=200)
            plt.close()
        
            fig,ax = plotsingle(xsplt,ysplt,xlim,ylim,Dmhmaxbtref2D,clevels,cmap,False,False,False,None,None)
            figfile=runname+'_'+treftimestring+'_Dmh_max_before.eps'
            plt.savefig(imagesavedir+'/'+figfile,dpi=200)
            plt.close()
        
            fig,ax = plotsingle(xsplt,ysplt,xlim,ylim,alpharmaxbtref2D,clevels,cmap,False,False,False,None,None)
            figfile=runname+'_'+treftimestring+'_alphar_max_before.eps'
            plt.savefig(imagesavedir+'/'+figfile,dpi=200)
            plt.close()
        
            fig,ax = plotsingle(xsplt,ysplt,xlim,ylim,alphahmaxbtref2D,clevels,cmap,False,False,False,None,None)
            figfile=runname+'_'+treftimestring+'_alphah_max_before.eps'
            plt.savefig(imagesavedir+'/'+figfile,dpi=200)
            plt.close()
        
            clevels = N.arange(0.0,15.1,0.1)
            cmap = cm.Blues
        
            # Average qc,qr,qh prior to reference time (g/kg)
            fig,ax = plotsingle(xsplt,ysplt,xlim,ylim,qravgbtref2D,clevels,cmap,False,False,False,None,None)
            figfile=runname+'_'+treftimestring+'_qr_avg_before.eps'
            plt.savefig(imagesavedir+'/'+figfile,dpi=200)
            plt.close()
        
            fig,ax = plotsingle(xsplt,ysplt,xlim,ylim,qcavgbtref2D,clevels,cmap,False,False,False,None,None)
            figfile=runname+'_'+treftimestring+'_qc_avg_before.eps'
            plt.savefig(imagesavedir+'/'+figfile,dpi=200)
            plt.close()
        
            fig,ax = plotsingle(xsplt,ysplt,xlim,ylim,qhavgbtref2D,clevels,cmap,False,False,False,None,None)
            figfile=runname+'_'+treftimestring+'_qh_avg_before.eps'
            plt.savefig(imagesavedir+'/'+figfile,dpi=200)
            plt.close()
        
            clevels = N.arange(0.0,0.011,0.001)
            cmap = cm.Blues
        
            # Average evapqc,evapqr,meltqh prior to reference time (g/kg/s)
            fig,ax = plotsingle(xsplt,ysplt,xlim,ylim,evapqravgbtref2D,clevels,cmap,False,False,False,None,None)
            figfile=runname+'_'+treftimestring+'_evapqr_avg_before.eps'
            plt.savefig(imagesavedir+'/'+figfile,dpi=200)
            plt.close()
        
            fig,ax = plotsingle(xsplt,ysplt,xlim,ylim,evapqcavgbtref2D,clevels,cmap,False,False,False,None,None)
            figfile=runname+'_'+treftimestring+'_evapqc_avg_before.eps'
            plt.savefig(imagesavedir+'/'+figfile,dpi=200)
            plt.close()
        
            fig,ax = plotsingle(xsplt,ysplt,xlim,ylim,meltqhavgbtref2D,clevels,cmap,False,False,False,None,None)
            figfile=runname+'_'+treftimestring+'_meltqh_avg_before.eps'
            plt.savefig(imagesavedir+'/'+figfile,dpi=200)
            plt.close()
        
            clevels = N.arange(0.0,10.0,0.1)
            cmap = cm.Blues
        
            # Average Dmr, Dmh, alphar, alphah prior to reference time (mm, unitless)
            fig,ax = plotsingle(xsplt,ysplt,xlim,ylim,Dmravgbtref2D,clevels,cmap,False,False,False,None,None)
            figfile=runname+'_'+treftimestring+'_Dmr_avg_before.eps'
            plt.savefig(imagesavedir+'/'+figfile,dpi=200)
            plt.close()
        
            fig,ax = plotsingle(xsplt,ysplt,xlim,ylim,Dmhavgbtref2D,clevels,cmap,False,False,False,None,None)
            figfile=runname+'_'+treftimestring+'_Dmh_avg_before.eps'
            plt.savefig(imagesavedir+'/'+figfile,dpi=200)
            plt.close()
        
            fig,ax = plotsingle(xsplt,ysplt,xlim,ylim,alpharavgbtref2D,clevels,cmap,False,False,False,None,None)
            figfile=runname+'_'+treftimestring+'_alphar_avg_before.eps'
            plt.savefig(imagesavedir+'/'+figfile,dpi=200)
            plt.close()
        
            fig,ax = plotsingle(xsplt,ysplt,xlim,ylim,alphahavgbtref2D,clevels,cmap,False,False,False,None,None)
            figfile=runname+'_'+treftimestring+'_alphah_avg_before.eps'
            plt.savefig(imagesavedir+'/'+figfile,dpi=200)
            plt.close()
        
            clevels = N.arange(0.0,10.0,0.1)
            cmap = cm.Blues
        
            # Min Dmr, Dmh, alphar, alphah prior to reference time (mm, unitless)
            fig,ax = plotsingle(xsplt,ysplt,xlim,ylim,Dmrminbtref2D,clevels,cmap,False,False,False,None,None)
            figfile=runname+'_'+treftimestring+'_Dmr_min_before.eps'
            plt.savefig(imagesavedir+'/'+figfile,dpi=200)
            plt.close()
        
            fig,ax = plotsingle(xsplt,ysplt,xlim,ylim,Dmhminbtref2D,clevels,cmap,False,False,False,None,None)
            figfile=runname+'_'+treftimestring+'_Dmh_min_before.eps'
            plt.savefig(imagesavedir+'/'+figfile,dpi=200)
            plt.close()
        
            fig,ax = plotsingle(xsplt,ysplt,xlim,ylim,alpharminbtref2D,clevels,cmap,False,False,False,None,None)
            figfile=runname+'_'+treftimestring+'_alphar_min_before.eps'
            plt.savefig(imagesavedir+'/'+figfile,dpi=200)
            plt.close()
        
            fig,ax = plotsingle(xsplt,ysplt,xlim,ylim,alphahminbtref2D,clevels,cmap,False,False,False,None,None)
            figfile=runname+'_'+treftimestring+'_alphah_min_before.eps'
            plt.savefig(imagesavedir+'/'+figfile,dpi=200)
            plt.close()
    
        if(plot_wforce):
            print "Creating plots for vertical velocity forcing for time "+str(tref)
            # Reshape arrays for plotting
        
            wforceminbtref2D = wforceminbtref.reshape((len(wforce_varnames),int(N.sqrt(numtrajcs)),-1))[:,:jmax-jmin+1,:imax-imin+1]
            wforceminatref2D = wforceminatref.reshape((len(wforce_varnames),int(N.sqrt(numtrajcs)),-1))[:,:jmax-jmin+1,:imax-imin+1]
        
            wforcemaxbtref2D = wforcemaxbtref.reshape((len(wforce_varnames),int(N.sqrt(numtrajcs)),-1))[:,:jmax-jmin+1,:imax-imin+1]
            wforcemaxatref2D = wforcemaxatref.reshape((len(wforce_varnames),int(N.sqrt(numtrajcs)),-1))[:,:jmax-jmin+1,:imax-imin+1]
        
            wforceavgbtref2D = wforceavgbtref.reshape((len(wforce_varnames),int(N.sqrt(numtrajcs)),-1))[:,:jmax-jmin+1,:imax-imin+1]
            wforceavgatref2D = wforceavgatref.reshape((len(wforce_varnames),int(N.sqrt(numtrajcs)),-1))[:,:jmax-jmin+1,:imax-imin+1]
        
            # Now for the plotting!
        
            clevels = N.arange(-0.5,0.51,0.05)
            cmap = cm.RdBu_r

            # Minimum buoyancy acceleration prior to reference time (m/s^2)
            fig,ax = plotsingle(xsplt,ysplt,xlim,ylim,wforceminbtref2D[0,:],clevels,cmap,False,False,True,qrintrp2D*1000.,[0.01,0.1,0.5,1.0,5.0])
            Q = plt.quiver(xsplt,ysplt,usplt,vsplt,pivot='middle', minshaft=1, scale=150, units='inches',headwidth=4, width=.01,color='k')
            if(masked_tornado and numtortrajcs > 0):
                ax.scatter(xvalstor,yvalstor,facecolors='none',edgecolors='g',marker='o')
            figfile=runname+'_'+treftimestring+'_bforce_min_before.eps'
            plt.savefig(imagesavedir+'/'+figfile,dpi=200)
            plt.close()
        
            # Minimum dynamic pressure gradient acceleration prior to reference time (m/s^2)
            fig,ax = plotsingle(xsplt,ysplt,xlim,ylim,wforceminbtref2D[1,:],clevels,cmap,False,False,True,qrintrp2D*1000.,[0.01,0.1,0.5,1.0,5.0])
            Q = plt.quiver(xsplt,ysplt,usplt,vsplt,pivot='middle', minshaft=1, scale=150, units='inches',headwidth=4, width=.01,color='k')
            if(masked_tornado and numtortrajcs > 0):
                ax.scatter(xvalstor,yvalstor,facecolors='none',edgecolors='g',marker='o')
            figfile=runname+treftimestring+'_ppgfd_min_before.eps'
            plt.savefig(imagesavedir+'/'+figfile,dpi=200)
            plt.close()

            # Minimum total vertical acceleration prior to reference time (m/s^2)
            fig,ax = plotsingle(xsplt,ysplt,xlim,ylim,wforceminbtref2D[2,:],clevels,cmap,False,False,True,qrintrp2D*1000.,[0.01,0.1,0.5,1.0,5.0])
            Q = plt.quiver(xsplt,ysplt,usplt,vsplt,pivot='middle', minshaft=1, scale=150, units='inches',headwidth=4, width=.01,color='k')
            if(masked_tornado and numtortrajcs > 0):
                ax.scatter(xvalstor,yvalstor,facecolors='none',edgecolors='g',marker='o')
            figfile=runname+'_'+treftimestring+'_waccel_min_before.eps'
            plt.savefig(imagesavedir+'/'+figfile,dpi=200)
            plt.close()
        
            clevels = N.arange(-0.5,0.51,0.05)
            cmap = cm.RdBu_r

            # Maximum buoyancy acceleration prior to reference time (m/s^2)
            fig,ax = plotsingle(xsplt,ysplt,xlim,ylim,wforcemaxbtref2D[0,:],clevels,cmap,False,False,True,qrintrp2D*1000.,[0.01,0.1,0.5,1.0,5.0])
            Q = plt.quiver(xsplt,ysplt,usplt,vsplt,pivot='middle', minshaft=1, scale=150, units='inches',headwidth=4, width=.01,color='k')
            if(masked_tornado and numtortrajcs > 0):
                ax.scatter(xvalstor,yvalstor,facecolors='none',edgecolors='g',marker='o')
            figfile=runname+'_'+treftimestring+'_bforce_max_before.eps'
            plt.savefig(imagesavedir+'/'+figfile,dpi=200)
            plt.close()
        
            # Maximum dynamic pressure gradient acceleration prior to reference time (m/s^2)
            fig,ax = plotsingle(xsplt,ysplt,xlim,ylim,wforcemaxbtref2D[1,:],clevels,cmap,False,False,True,qrintrp2D*1000.,[0.01,0.1,0.5,1.0,5.0])
            Q = plt.quiver(xsplt,ysplt,usplt,vsplt,pivot='middle', minshaft=1, scale=150, units='inches',headwidth=4, width=.01,color='k')
            if(masked_tornado and numtortrajcs > 0):
                ax.scatter(xvalstor,yvalstor,facecolors='none',edgecolors='g',marker='o')
            figfile=runname+'_'+treftimestring+'_ppgfd_max_before.eps'
            plt.savefig(imagesavedir+'/'+figfile,dpi=200)
            plt.close()
 
            # Maximum total vertical acceleration prior to reference time (m/s^2)
            fig,ax = plotsingle(xsplt,ysplt,xlim,ylim,wforcemaxbtref2D[2,:],clevels,cmap,False,False,True,qrintrp2D*1000.,[0.01,0.1,0.5,1.0,5.0])
            Q = plt.quiver(xsplt,ysplt,usplt,vsplt,pivot='middle', minshaft=1, scale=150, units='inches',headwidth=4, width=.01,color='k')
            if(masked_tornado and numtortrajcs > 0):
                ax.scatter(xvalstor,yvalstor,facecolors='none',edgecolors='g',marker='o')
            figfile=runname+'_'+treftimestring+'_waccel_max_before.eps'
            plt.savefig(imagesavedir+'/'+figfile,dpi=200)
            plt.close()
        
            clevels = N.arange(-0.1,0.105,0.005)
            cmap = cm.RdBu_r
            norm = matplotlib.colors.normalize(-0.1,0.1)

            # Average buoyancy acceleration prior to reference time (m/s^2)
            fig,ax = plotsingle(xsplt,ysplt,xlim,ylim,wforceavgbtref2D[0,:],clevels,cmap,False,False,True,qrintrp2D*1000.,[0.01,0.1,0.5,1.0,5.0])
            Q = plt.quiver(xsplt,ysplt,usplt,vsplt,pivot='middle', minshaft=1, scale=150, units='inches',headwidth=4, width=.01,color='k')
            if(masked_tornado and numtortrajcs > 0):
                ax.scatter(xvalstor,yvalstor,facecolors='none',edgecolors='g',marker='o')
            figfile=runname+'_'+treftimestring+'_bforce_avg_before.eps'
            plt.savefig(imagesavedir+'/'+figfile,dpi=200)
            plt.close()

            # Average dynamic pressure gradient acceleration prior to reference time (m/s^2)
            fig,ax = plotsingle(xsplt,ysplt,xlim,ylim,wforceavgbtref2D[1,:],clevels,cmap,False,False,True,qrintrp2D*1000.,[0.01,0.1,0.5,1.0,5.0])
            Q = plt.quiver(xsplt,ysplt,usplt,vsplt,pivot='middle', minshaft=1, scale=150, units='inches',headwidth=4, width=.01,color='k')
            if(masked_tornado and numtortrajcs > 0):
                ax.scatter(xvalstor,yvalstor,facecolors='none',edgecolors='g',marker='o')
            figfile=runname+'_'+treftimestring+'_ppgfd_avg_before.eps'
            plt.savefig(imagesavedir+'/'+figfile,dpi=200)
            plt.close()
 
            # Average total vertical acceleration prior to reference time (m/s^2)
            fig,ax = plotsingle(xsplt,ysplt,xlim,ylim,wforceavgbtref2D[2,:],clevels,cmap,False,False,True,qrintrp2D*1000.,[0.01,0.1,0.5,1.0,5.0])
            Q = plt.quiver(xsplt,ysplt,usplt,vsplt,pivot='middle', minshaft=1, scale=150, units='inches',headwidth=4, width=.01,color='k')
            if(masked_tornado and numtortrajcs > 0):
                ax.scatter(xvalstor,yvalstor,facecolors='none',edgecolors='g',marker='o')
            figfile=runname+'_'+treftimestring+'_waccel_avg_before.eps'
            plt.savefig(imagesavedir+'/'+figfile,dpi=200)
            plt.close()
    
    if(plot_integrated_theta):
        # Computed integrated theta
        
        evapqccoolintrp = -DSDintrp[DSD_varnames.index('evapqc'),:trefindex+1,:]*Lv/cp
        evapqrcoolintrp = -DSDintrp[DSD_varnames.index('evapqr'),:trefindex+1,:]*Lv/cp
        meltqhcoolintrp = -DSDintrp[DSD_varnames.index('meltqh'),:trefindex+1,:]*Lf/cp
        
        evapqccoolinteg = integrate_trajcs(varintrp[varnames.index('pt'),0,:],
                        evapqccoolintrp,tintv)[-1]
        
        evapqrcoolinteg  = integrate_trajcs(varintrp[varnames.index('pt'),0,:],
                        evapqrcoolintrp,tintv)[-1]
        
        meltqhcoolinteg =integrate_trajcs(varintrp[varnames.index('pt'),0,:],
                        meltqhcoolintrp,tintv)[-1]
                        
        cooltotinteg = integrate_trajcs(varintrp[varnames.index('pt'),0,:],
                        (evapqccoolintrp+evapqrcoolintrp+meltqhcoolintrp),tintv)[-1]
        
        # Find trajectory with min integrated theta
        
        traj_minthetaacc = N.argmin(cooltotinteg)
        mintheta = N.min(cooltotinteg)
        
        print "Number of trajectory with minimum integrated theta is ",traj_minthetaacc
        print "The min theta value is: ",mintheta
                        
        # Reshape interpolated trajectory array at reference time into 2D (assumes trajectory reference domain is square)
        
        varintrp2D = varintrp[:,trefindex,:].reshape((len(varnames),int(N.sqrt(numtrajcs)),-1))
        
        evapqccoolinteg2D = evapqccoolinteg[:].reshape(int(N.sqrt(numtrajcs)),-1)
        evapqrcoolinteg2D  = evapqrcoolinteg[:].reshape(int(N.sqrt(numtrajcs)),-1)
        meltqhcoolinteg2D = meltqhcoolinteg[:].reshape(int(N.sqrt(numtrajcs)),-1)
        cooltotinteg2D = cooltotinteg[:].reshape(int(N.sqrt(numtrajcs)),-1)
    
        # Filter with raymond 2D lowpass?
        #bforceinteg2D = raymond.raymond2d_lowpass(bforceinteg2D,1000.0)
        #ppgfdinteg2D = raymond.raymond2d_lowpass(ppgfdinteg2D,1000.0)
        #wacctotinteg2D = raymond.raymond2d_lowpass(wacctotinteg2D,1000.0)

        # Now make contour plots of these integrated forcings
        
        clevels = N.arange(290.0,331.0,1.0)
        cmap = cm.jet
        norm = matplotlib.colors.normalize(290.0,331.0)
        
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(511)
        ax2 = fig1.add_subplot(512)
        ax3 = fig1.add_subplot(513)
        ax4 = fig1.add_subplot(514)
        ax5 = fig1.add_subplot(515)
        
        # Interpolated theta (K)
        CS1 = ax1.contourf(xsplt,ysplt,varintrp2D[varnames.index('pt')],levels=clevels,cmap=cmap,norm=norm)
        print varintrp2D[varnames.index('pt')]
        
        # Integrated cloud evaporative cooling (K)
        CS2 = ax2.contourf(xsplt,ysplt,evapqccoolinteg2D,levels=clevels,cmap=cmap,norm=norm)
        #plt.colorbar(CS2,ax=ax2)
        
        # Integrated rain evaporative cooling (K)
        CS3 = ax3.contourf(xsplt,ysplt,evapqrcoolinteg2D,levels=clevels,cmap=cmap,norm=norm)
        
        # Integrated hail melting (K)
        CS4 = ax4.contourf(xsplt,ysplt,meltqhcoolinteg2D,levels=clevels,cmap=cmap,norm=norm)
        #plt.colorbar(CS4,ax=ax4)
        
        # Integrated total cooling (K)
        CS5 = ax5.contourf(xsplt,ysplt,cooltotinteg2D,levels=clevels,cmap=cmap,norm=norm)
        plt.colorbar(CS5,ax=ax5)
        
        plt.scatter([xvals[trefindex,traj_minthetaacc]],[yvals[trefindex,traj_minthetaacc]])
        
        formatter=plt.FuncFormatter(mtokm)
        ax1.xaxis.set_major_formatter(formatter)
        ax1.xaxis.set_major_locator(plt.MultipleLocator(base=1000.0))
        ax1.yaxis.set_major_formatter(formatter)
        ax1.yaxis.set_major_locator(plt.MultipleLocator(base=1000.0))
        ax1.set_aspect('equal')
        ax2.xaxis.set_major_formatter(formatter)
        ax2.xaxis.set_major_locator(plt.MultipleLocator(base=1000.0))
        ax2.yaxis.set_major_formatter(formatter)
        ax2.yaxis.set_major_locator(plt.MultipleLocator(base=1000.0))
        ax2.set_aspect('equal')
        ax3.xaxis.set_major_formatter(formatter)
        ax3.xaxis.set_major_locator(plt.MultipleLocator(base=1000.0))
        ax3.yaxis.set_major_formatter(formatter)
        ax3.yaxis.set_major_locator(plt.MultipleLocator(base=1000.0))
        ax3.set_aspect('equal')
        ax4.xaxis.set_major_formatter(formatter)
        ax4.xaxis.set_major_locator(plt.MultipleLocator(base=1000.0))
        ax4.yaxis.set_major_formatter(formatter)
        ax4.yaxis.set_major_locator(plt.MultipleLocator(base=1000.0))
        ax4.set_aspect('equal')
        ax5.xaxis.set_major_formatter(formatter)
        ax5.xaxis.set_major_locator(plt.MultipleLocator(base=1000.0))
        ax5.yaxis.set_major_formatter(formatter)
        ax5.yaxis.set_major_locator(plt.MultipleLocator(base=1000.0))
        ax5.set_aspect('equal')
        
        figfile=runname+'_'+treftimestring+'_thetainteg_before.eps'
        plt.savefig(imagesavedir+'/'+figfile,dpi=200)
             
    if(plot_integrated_wforcing):
        # Compute integrated vertical velocity forcing along trajectories up to reference time
        if(False):
            # Filter with raymond low-pass filter first?
            print "Filtering in 1D with raymond_lowpass filter"
            bforcefilt=[]
            ppgfdfilt=[]
            waccfilt=[]
            for traj in xrange(numtrajcs):
                bforcefilt.append(raymond.raymond1d_lowpass(wforceintrp[wforce_varnames.index('bforce'),:trefindex,traj],1000.0))
                ppgfdfilt.append(raymond.raymond1d_lowpass(wforceintrp[wforce_varnames.index('ppgf_d'),:trefindex,traj],1000.0))
                waccfilt.append(raymond.raymond1d_lowpass(wforceintrp[wforce_varnames.index('waccel'),:trefindex,traj],1000.0)) 

            bforcefilt = N.array(bforcefilt).swapaxes(0,1)
            ppgfdfilt = N.array(ppgfdfilt).swapaxes(0,1)
            waccfilt = N.array(waccfilt).swapaxes(0,1)
        else:
            bforcefilt = wforceintrp[wforce_varnames.index('bforce'),:trefindex,:]
            ppgfdfilt =  wforceintrp[wforce_varnames.index('ppgf_d'),:trefindex,:]
            waccfilt =   wforceintrp[wforce_varnames.index('waccel'),:trefindex,:]

        bforceinteg = integrate_trajcs(varintrp[varnames.index('w'),0,:],
                        bforcefilt,tintv)[-1]
        
        ppgfdinteg  = integrate_trajcs(varintrp[varnames.index('w'),0,:],
                        ppgfdfilt,tintv)[-1]
        
        wacctotinteg =integrate_trajcs(varintrp[varnames.index('w'),0,:],
                        waccfilt,tintv)[-1]
                        
        # Find trajectory with max integrated w
        
        traj_minwacc = N.argmin(wacctotinteg)
        
        print "Number of trajectory with minimum integrated vertical velocity is ",traj_minwacc
                        
        # Reshape interpolated trajectory array at reference time into 2D (assumes trajectory reference domain is square)
        
        varintrp2D = varintrp[:,trefindex,:].reshape((len(varnames),int(N.sqrt(numtrajcs)),-1))
        
        bforceinteg2D = bforceinteg[:].reshape(int(N.sqrt(numtrajcs)),-1)
        ppgfdinteg2D  = ppgfdinteg[:].reshape(int(N.sqrt(numtrajcs)),-1)
        wacctotinteg2D = wacctotinteg[:].reshape(int(N.sqrt(numtrajcs)),-1)
    
        # Filter with raymond 2D lowpass?
        #bforceinteg2D = raymond.raymond2d_lowpass(bforceinteg2D,1000.0)
        #ppgfdinteg2D = raymond.raymond2d_lowpass(ppgfdinteg2D,1000.0)
        #wacctotinteg2D = raymond.raymond2d_lowpass(wacctotinteg2D,1000.0)

        # Now make contour plots of these integrated forcings
        
        clevels = N.arange(-10.0,10.1,0.1)
        cmap = cm.RdBu_r
        norm = matplotlib.colors.normalize(-10.0,10.0)
        
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(221)
        ax2 = fig1.add_subplot(222)
        ax3 = fig1.add_subplot(223)
        ax4 = fig1.add_subplot(224)
        
        # Interpolated vertical velocity (m/s)
        CS1 = ax1.contourf(xsplt,ysplt,varintrp2D[0],levels=clevels,cmap=cmap,norm=norm)
        
        # Integrated buoyancy forcing (m/s)
        CS2 = ax2.contourf(xsplt,ysplt,bforceinteg2D,levels=clevels,cmap=cmap,norm=norm)
        plt.colorbar(CS2,ax=ax2)
        
        # Integrated dynamic perturbation pressure gradient forcing (m/s)
        CS3 = ax3.contourf(xsplt,ysplt,ppgfdinteg2D,levels=clevels,cmap=cmap,norm=norm)
        
        # Integrated total vertical velocity forcing (m/s)
        CS4 = ax4.contourf(xsplt,ysplt,wacctotinteg2D,levels=clevels,cmap=cmap,norm=norm)
        plt.colorbar(CS4,ax=ax4)
        plt.scatter([xvals[trefindex,traj_minwacc]],[yvals[trefindex,traj_minwacc]])
        
        formatter=plt.FuncFormatter(mtokm)
        ax1.xaxis.set_major_formatter(formatter)
        ax1.xaxis.set_major_locator(plt.MultipleLocator(base=1000.0))
        ax1.yaxis.set_major_formatter(formatter)
        ax1.yaxis.set_major_locator(plt.MultipleLocator(base=1000.0))
        ax1.set_aspect('equal')
        ax2.xaxis.set_major_formatter(formatter)
        ax2.xaxis.set_major_locator(plt.MultipleLocator(base=1000.0))
        ax2.yaxis.set_major_formatter(formatter)
        ax2.yaxis.set_major_locator(plt.MultipleLocator(base=1000.0))
        ax2.set_aspect('equal')
        ax3.xaxis.set_major_formatter(formatter)
        ax3.xaxis.set_major_locator(plt.MultipleLocator(base=1000.0))
        ax3.yaxis.set_major_formatter(formatter)
        ax3.yaxis.set_major_locator(plt.MultipleLocator(base=1000.0))
        ax3.set_aspect('equal')
        ax4.xaxis.set_major_formatter(formatter)
        ax4.xaxis.set_major_locator(plt.MultipleLocator(base=1000.0))
        ax4.yaxis.set_major_formatter(formatter)
        ax4.yaxis.set_major_locator(plt.MultipleLocator(base=1000.0))
        ax4.set_aspect('equal')
        
        figfile=runname+'_'+treftimestring+'_winteg_before.eps'
        plt.savefig(imagesavedir+'/'+figfile,dpi=200)
    
    if(plot_TS and not loadTS):
        if(plot_var):
            print "Computing time series of trajectory averages for variables"
            if(t == 0):
                zmaxbTS = N.zeros((tref_list.size))
                zmaxaTS = N.zeros((tref_list.size))
            
                varminbTS = N.zeros((tref_list.size,len(varnames)))
                varminaTS = N.zeros_like(varminbTS)
                varmaxbTS = N.zeros_like(varminbTS)
                varmaxaTS = N.zeros_like(varminbTS)
                varavgbTS = N.zeros_like(varminbTS)
                varavgaTS = N.zeros_like(varminbTS)
            
                varminbTS = N.zeros((tref_list.size,len(varnames)))
                varminaTS = N.zeros_like(varminbTS)
                varmaxbTS = N.zeros_like(varminbTS)
                varmaxaTS = N.zeros_like(varminbTS)
                varavgbTS = N.zeros_like(varminbTS)
                varavgaTS = N.zeros_like(varminbTS)
            
                numtortrajcsTS = N.zeros_like(zmaxbTS)
        
                if(usetorreftime):
                    varminattorTS = N.zeros_like(varminbTS)
                    varavgattorTS = N.zeros_like(varminbTS)
                    varmaxattorTS = N.zeros_like(varminbTS)
                    varminbtorTS = N.zeros_like(varminbTS)
                    varavgbtorTS = N.zeros_like(varminbTS)
                    varmaxbtorTS = N.zeros_like(varminbTS)
                    varminatorTS = N.zeros_like(varminbTS)
                    varavgatorTS = N.zeros_like(varminbTS)
                    varmaxatorTS = N.zeros_like(varminbTS)
        
            if(masked_tornado): # Define masked arrays that mask out trajectories that don't enter tornado and perform averaging
                numtortrajcsTS[t] = numtortrajcs
            
                zmaxbTS[t] = N.ma.masked_array(zmaxbtref,mask=tormask).mean()
                zmaxaTS[t] = N.ma.masked_array(zmaxatref,mask=tormask).mean()
            
                for i,var in enumerate(varnames):
                    varminbTS[t,i] = N.ma.masked_array(varminbtref[i,:],mask=tormask).mean()
                    varminaTS[t,i] = N.ma.masked_array(varminatref[i,:],mask=tormask).mean()
                    varmaxbTS[t,i] = N.ma.masked_array(varmaxbtref[i,:],mask=tormask).mean()
                    varmaxaTS[t,i] = N.ma.masked_array(varmaxatref[i,:],mask=tormask).mean()
                    varavgbTS[t,i] = N.ma.masked_array(varavgbtref[i,:],mask=tormask).mean()
                    varavgaTS[t,i] = N.ma.masked_array(varavgatref[i,:],mask=tormask).mean()
            
                # Compute min/mean/max relative to trajectory entering tornado
                if(usetorreftime):
                    for i,var in enumerate(varnames):
                        varminattorTS[t,i] = N.ma.masked_array(valattor[i,:],mask=tormask).min()
                        varavgattorTS[t,i] = N.ma.masked_array(valattor[i,:],mask=tormask).mean()
                        varmaxattorTS[t,i] = N.ma.masked_array(valattor[i,:],mask=tormask).max()
                        varminbtorTS[t,i] = N.ma.masked_array(varminbtortref[i,:],mask=tormask).mean()
                        varavgbtorTS[t,i] = N.ma.masked_array(varavgbtortref[i,:],mask=tormask).mean()
                        varmaxbtorTS[t,i] = N.ma.masked_array(varmaxbtortref[i,:],mask=tormask).mean()
                        varminatorTS[t,i] = N.ma.masked_array(varminatortref[i,:],mask=tor5minmask).mean()
                        varavgatorTS[t,i] = N.ma.masked_array(varavgatortref[i,:],mask=tor5minmask).mean()
                        varmaxatorTS[t,i] = N.ma.masked_array(varmaxatortref[i,:],mask=tor5minmask).mean()
                                                                                                                
        if(plot_wforce):
            print "Computing time series of trajectory averages for vertical velocity forcing"
            if(t == 0):
                wforceminbTS = N.zeros((tref_list.size,len(wforce_varnames)))
                wforceminaTS = N.zeros_like(wforceminbTS)
                wforcemaxbTS = N.zeros_like(wforceminbTS)
                wforcemaxaTS = N.zeros_like(wforceminbTS)
                wforceavgbTS = N.zeros_like(wforceminbTS)
                wforceavgaTS = N.zeros_like(wforceminbTS)
            
                if(usetorreftime):
                    wforceminattorTS = N.zeros_like(wforceminbTS)
                    wforceavgattorTS = N.zeros_like(wforceminbTS)
                    wforcemaxattorTS = N.zeros_like(wforceminbTS)
                    wforceminbtorTS = N.zeros_like(wforceminbTS)
                    wforceavgbtorTS = N.zeros_like(wforceminbTS)
                    wforcemaxbtorTS = N.zeros_like(wforceminbTS)
                    wforceminatorTS = N.zeros_like(wforceminbTS)
                    wforceavgatorTS = N.zeros_like(wforceminbTS)
                    wforcemaxatorTS = N.zeros_like(wforceminbTS)
            
            if(masked_tornado): # Define masked arrays that mask out trajectories that don't enter tornado and perform averaging
                for i,var in enumerate(wforce_varnames):
                    wforceminbTS[t,i] = N.ma.masked_array(wforceminbtref[i,:],mask=tormask).mean()
                    wforceminaTS[t,i] = N.ma.masked_array(wforceminatref[i,:],mask=tormask).mean()
                    wforcemaxbTS[t,i] = N.ma.masked_array(wforcemaxbtref[i,:],mask=tormask).mean()
                    wforcemaxaTS[t,i] = N.ma.masked_array(wforcemaxatref[i,:],mask=tormask).mean()
                    wforceavgbTS[t,i] = N.ma.masked_array(wforceavgbtref[i,:],mask=tormask).mean()
                    wforceavgaTS[t,i] = N.ma.masked_array(wforceavgatref[i,:],mask=tormask).mean()
    
                if(usetorreftime):
                    for i,var in enumerate(wforce_varnames):    
                        wforceminattorTS[t,i] = N.ma.masked_array(wforceattor[i,:],mask=tormask).min()
                        wforceavgattorTS[t,i] = N.ma.masked_array(wforceattor[i,:],mask=tormask).mean()
                        wforcemaxattorTS[t,i] = N.ma.masked_array(wforceattor[i,:],mask=tormask).max()
                        wforceminbtorTS[t,i] = N.ma.masked_array(wforceminbtortref[i,:],mask=tormask).mean()
                        wforceavgbtorTS[t,i] = N.ma.masked_array(wforceavgbtortref[i,:],mask=tormask).mean()
                        wforcemaxbtorTS[t,i] = N.ma.masked_array(wforcemaxbtortref[i,:],mask=tormask).mean()
                        wforceminatorTS[t,i] = N.ma.masked_array(wforceminatortref[i,:],mask=tor5minmask).mean()
                        wforceavgatorTS[t,i] = N.ma.masked_array(wforceavgatortref[i,:],mask=tor5minmask).mean()
                        wforcemaxatorTS[t,i] = N.ma.masked_array(wforcemaxatortref[i,:],mask=tor5minmask).mean()
                    
        if(plot_DSD):
            print "Computing time series of trajectory averages for DSD variables"
            if(t == 0):
                qrmaxbTS = N.zeros((tref_list.size))
                qcmaxbTS = N.zeros_like(qrmaxbTS)
                qhmaxbTS = N.zeros_like(qrmaxbTS)
                evapqcmaxbTS = N.zeros_like(qrmaxbTS)
                evapqrmaxbTS = N.zeros_like(qrmaxbTS)
                meltqhmaxbTS = N.zeros_like(qrmaxbTS)
                DmrmaxbTS = N.zeros_like(qrmaxbTS)
                DmhmaxbTS = N.zeros_like(qrmaxbTS)
                alpharmaxbTS = N.zeros_like(qrmaxbTS)
                alphahmaxbTS = N.zeros_like(qrmaxbTS)
                qravgbTS = N.zeros_like(qrmaxbTS)
                qcavgbTS = N.zeros_like(qrmaxbTS)
                qhavgbTS = N.zeros_like(qrmaxbTS)
                evapqcavgbTS = N.zeros_like(qrmaxbTS)
                evapqravgbTS = N.zeros_like(qrmaxbTS)
                meltqhavgbTS = N.zeros_like(qrmaxbTS)
                DmravgbTS = N.zeros_like(qrmaxbTS)
                DmhavgbTS = N.zeros_like(qrmaxbTS)
                alpharavgbTS = N.zeros_like(qrmaxbTS)
                alphahavgbTS = N.zeros_like(qrmaxbTS)
                DmrminbTS = N.zeros_like(qrmaxbTS)
                DmhminbTS = N.zeros_like(qrmaxbTS)
                alpharminbTS = N.zeros_like(qrmaxbTS)
                alphahminbTS = N.zeros_like(qrmaxbTS)
            
                if(usetorreftime):
                    qrminbtorTS = N.zeros_like(qrmaxbTS) 
                    qrmaxbtorTS = N.zeros_like(qrmaxbTS)
                    qravgbtorTS = N.zeros_like(qrmaxbTS)
                    qcminbtorTS = N.zeros_like(qrmaxbTS)
                    qcmaxbtorTS = N.zeros_like(qrmaxbTS)
                    qcavgbtorTS = N.zeros_like(qrmaxbTS)
                    qhminbtorTS = N.zeros_like(qrmaxbTS)
                    qhmaxbtorTS = N.zeros_like(qrmaxbTS)
                    qhavgbtorTS = N.zeros_like(qrmaxbTS)        
                    evapqcmaxbtorTS = N.zeros_like(qrmaxbTS)
                    evapqcavgbtorTS = N.zeros_like(qrmaxbTS)
                    evapqrmaxbtorTS = N.zeros_like(qrmaxbTS)
                    evapqravgbtorTS = N.zeros_like(qrmaxbTS)
                    meltqhmaxbtorTS = N.zeros_like(qrmaxbTS)
                    meltqhavgbtorTS = N.zeros_like(qrmaxbTS)
                    DmrmaxbtorTS = N.zeros_like(qrmaxbTS)
                    DmravgbtorTS = N.zeros_like(qrmaxbTS)
                    DmrminbtorTS = N.zeros_like(qrmaxbTS)
                    DmhmaxbtorTS = N.zeros_like(qrmaxbTS)
                    DmhavgbtorTS = N.zeros_like(qrmaxbTS)
                    DmhminbtorTS = N.zeros_like(qrmaxbTS)
                    alpharmaxbtorTS = N.zeros_like(qrmaxbTS)
                    alpharavgbtorTS = N.zeros_like(qrmaxbTS)
                    alpharminbtorTS = N.zeros_like(qrmaxbTS)
                    alphahmaxbtorTS = N.zeros_like(qrmaxbTS)
                    alphahavgbtorTS = N.zeros_like(qrmaxbTS)
                    alphahminbtorTS = N.zeros_like(qrmaxbTS)
                                
            qrmaxbTS[t] = N.ma.masked_array(qrmaxbtref,mask=tormask).mean()
            qcmaxbTS[t] = N.ma.masked_array(qcmaxbtref,mask=tormask).mean()
            qhmaxbTS[t] = N.ma.masked_array(qhmaxbtref,mask=tormask).mean()
            evapqcmaxbTS[t] = N.ma.masked_array(evapqcmaxbtref,mask=tormask).mean()
            evapqrmaxbTS[t] = N.ma.masked_array(evapqrmaxbtref,mask=tormask).mean()
            meltqhmaxbTS[t] = N.ma.masked_array(meltqhmaxbtref,mask=tormask).mean()
            DmrmaxbTS[t] = N.ma.masked_array(Dmrmaxbtref,mask=tormask).mean()
            DmhmaxbTS[t] = N.ma.masked_array(Dmhmaxbtref,mask=tormask).mean()
            alpharmaxbTS[t] = N.ma.masked_array(alpharmaxbtref,mask=tormask).mean()
            alphahmaxbTS[t] = N.ma.masked_array(alphahmaxbtref,mask=tormask).mean() 
            qravgbTS[t] = N.ma.masked_array(qravgbtref,mask=tormask).mean()
            qcavgbTS[t] = N.ma.masked_array(qcavgbtref,mask=tormask).mean() 
            qhavgbTS[t] = N.ma.masked_array(qhavgbtref,mask=tormask).mean() 
            evapqcavgbTS[t] = N.ma.masked_array(evapqcavgbtref,mask=tormask).mean() 
            evapqravgbTS[t] = N.ma.masked_array(evapqravgbtref,mask=tormask).mean()
            meltqhavgbTS[t] = N.ma.masked_array(meltqhavgbtref,mask=tormask).mean() 
            DmravgbTS[t] = N.ma.masked_array(Dmravgbtref,mask=tormask).mean()
            DmhavgbTS[t] = N.ma.masked_array(Dmhavgbtref,mask=tormask).mean() 
            alpharavgbTS[t] = N.ma.masked_array(alpharavgbtref,mask=tormask).mean() 
            alphahavgbTS[t] = N.ma.masked_array(alphahavgbtref,mask=tormask).mean() 
            DmrminbTS[t] = N.ma.masked_array(Dmrminbtref,mask=tormask).mean() 
            DmhminbTS[t] = N.ma.masked_array(Dmhminbtref,mask=tormask).mean() 
            alpharminbTS[t] = N.ma.masked_array(alpharminbtref,mask=tormask).mean() 
            alphahminbTS[t] = N.ma.masked_array(alphahminbtref,mask=tormask).mean()
        
            if(usetorreftime):
                qrmaxbtorTS[t] = N.ma.masked_array(qrmaxbtortref,mask=tormask).mean()
                qravgbtorTS[t] = N.ma.masked_array(qravgbtortref,mask=tormask).mean()
                qcmaxbtorTS[t] = N.ma.masked_array(qcmaxbtortref,mask=tormask).mean()
                qcavgbtorTS[t] = N.ma.masked_array(qcavgbtortref,mask=tormask).mean()
                qhmaxbtorTS[t] = N.ma.masked_array(qhmaxbtortref,mask=tormask).mean()
                qhavgbtorTS[t] = N.ma.masked_array(qhavgbtortref,mask=tormask).mean()  
                evapqcmaxbtorTS[t] = N.ma.masked_array(evapqcmaxbtortref,mask=tormask).mean()
                evapqcavgbtorTS[t] = N.ma.masked_array(evapqcavgbtortref,mask=tormask).mean()
                evapqrmaxbtorTS[t] = N.ma.masked_array(evapqrmaxbtortref,mask=tormask).mean()
                evapqravgbtorTS[t] = N.ma.masked_array(evapqravgbtortref,mask=tormask).mean()
                meltqhmaxbtorTS[t] = N.ma.masked_array(meltqhmaxbtortref,mask=tormask).mean()
                meltqhavgbtorTS[t] = N.ma.masked_array(meltqhavgbtortref,mask=tormask).mean()
                DmrmaxbtorTS[t] = N.ma.masked_array(Dmrmaxbtortref,mask=tormask).mean()
                DmravgbtorTS[t] = N.ma.masked_array(Dmravgbtortref,mask=tormask).mean()
                DmrminbtorTS[t] = N.ma.masked_array(Dmrminbtortref,mask=tormask).mean()
                DmhmaxbtorTS[t] = N.ma.masked_array(Dmhmaxbtortref,mask=tormask).mean()
                DmhavgbtorTS[t] = N.ma.masked_array(Dmhavgbtortref,mask=tormask).mean()
                DmhminbtorTS[t] = N.ma.masked_array(Dmhminbtortref,mask=tormask).mean()
                alpharmaxbtorTS[t] = N.ma.masked_array(alpharmaxbtortref,mask=tormask).mean()
                alpharavgbtorTS[t] = N.ma.masked_array(alpharavgbtortref,mask=tormask).mean()
                alpharminbtorTS[t] = N.ma.masked_array(alpharminbtortref,mask=tormask).mean()
                alphahmaxbtorTS[t] = N.ma.masked_array(alphahmaxbtortref,mask=tormask).mean()
                alphahavgbtorTS[t] = N.ma.masked_array(alphahavgbtortref,mask=tormask).mean()
                alphahminbtorTS[t] = N.ma.masked_array(alphahminbtortref,mask=tormask).mean()
    
    if(plot_hist):   # Histograms
        print "Creating histograms of a bunch of stuff"
        
        if(plot_var):
            zvals_hist = zvals_agl[0,:].flatten()
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.hist(zvals_hist,bins=N.linspace(0.,3000.,30),normed=False,log=True,color='b',alpha=1.0,orientation='horizontal')
            if(masked_tornado):
                zvals_tor_hist = N.ma.masked_array(zvals_agl[0,:],mask=tormask).flatten().compressed()
                if(numtortrajcs > 0):
                    ax.hist(zvals_tor_hist,bins=N.linspace(0.,3000.,30),normed=False,log=True,color='r',alpha=0.75,orientation='horizontal')
            ax.set_xlim(1.0,600.0)
            ax.set_xscale('log')
            ax.set_ylim(0.0,3000.0)
            ax.set_ylabel(r'$height (m)$')
            ax.set_xlabel(r'Trajectory count')

            fig.savefig(imagesavedir+'/'+runname+'_'+treftimestring+'_Z_hist.png',dpi=200,bbox_inches='tight')
        
            thetae_hist = varminbtref[varnames.index('pte'),:,:].flatten()
            ptprtr_hist = varminbtref[varnames.index('ptprtr'),:,:].flatten()
            qvprt_hist = varminbtref[varnames.index('qvprt'),:,:].flatten()
        
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.hist(thetae_hist,bins=N.linspace(310.,350.,40.),normed=False,log=False,color='b',alpha=1.0,orientation='vertical')
            if(masked_tornado):
                thetae_tor_hist = varminbtref[varnames.index('pte'),:,~tormask].flatten()
                if(numtortrajcs > 0):
                    ax.hist(thetae_tor_hist,bins=N.linspace(310.,350.,40.),normed=False,log=False,color='r',alpha=0.75,orientation='vertical')
            #ax.set_ylim(1.0,600.0)
            #ax.set_yscale('log')
            ax.set_xlim(310.,350.)
            ax.set_xlabel(r'$\theta_{e} (K)$')
            ax.set_ylabel(r'Trajectory points')

            fig.savefig(imagesavedir+'/'+runname+'_'+treftimestring+'_pte_hist.png',dpi=200,bbox_inches='tight')
        
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.hist(ptprtr_hist,bins=N.linspace(-10.0,10.0,40.),normed=False,log=False,color='b',alpha=1.0,orientation='vertical')
            if(masked_tornado):
                ptprtr_tor_hist = varminbtref[varnames.index('ptprtr'),:,~tormask].flatten()
                if(numtortrajcs > 0):
                    ax.hist(ptprtr_tor_hist,bins=N.linspace(-10.0,10.0,40.),normed=False,log=False,color='r',alpha=0.75,orientation='vertical')
            #ax.set_ylim(1.0,600.0)
            #ax.set_yscale('log')
            ax.set_xlim(-10.,10.)
            ax.set_xlabel(r"$\theta' (K)$")
            ax.set_ylabel(r'Trajectory points')

            fig.savefig(imagesavedir+'/'+runname+'_'+treftimestring+'_ptprtr_hist.png',dpi=200,bbox_inches='tight')
        
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.hist(qvprt_hist,bins=N.linspace(-0.01,0.01,40.),normed=False,log=False,color='b',alpha=1.0,orientation='vertical')
            if(masked_tornado):
                qvprt_tor_hist = varminbtref[varnames.index('qvprt'),:,~tormask].flatten()
                if(numtortrajcs > 0):
                    ax.hist(qvprt_tor_hist,bins=N.linspace(-0.01,0.01,40.),normed=False,log=False,color='r',alpha=0.75,orientation='vertical')
            #ax.set_ylim(1.0,600.0)
            #ax.set_yscale('log')
            ax.set_xlim(-0.01,0.01)
            ax.set_xlabel(r"$q_{v}' (kg kg^{-1})$'")
            ax.set_ylabel(r'Trajectory points')

            fig.savefig(imagesavedir+'/'+runname+'_'+treftimestring+'_pte_hist.png',dpi=200,bbox_inches='tight')
    
        if(plot_wforce):   # Histograms of w-forcing
            print "Creating histograms of w-forcing"
            for i,var in enumerate(wforce_varnames):
        
                varb_hist = wforcebtref[wforce_varnames.index(var),:,:].flatten()
                vara_hist = wforceatref[wforce_varnames.index(var),:,:].flatten()
        
                fig = plt.figure()
                ax = fig.add_subplot(111)
                #ax.hist(varb_hist,bins=N.linspace(-0.5,0.5,100.),normed=False,log=False,color='b',alpha=1.0,orientation='vertical')
                if(masked_tornado):
                    varb_tor_hist = wforcebtref[wforce_varnames.index(var),:,~tormask].flatten()
                    if(numtortrajcs > 0):
                        ax.hist(varb_tor_hist,bins=N.linspace(-0.5,0.5,100.),normed=False,log=False,color='r',alpha=0.75,orientation='vertical')
                #ax.set_ylim(1.0,600.0)
                #ax.set_yscale('log')
                ax.set_xlim(-0.5,0.5)
                ax.set_xlabel(var+r'$acceleration (m s^{-2})$ (before reference time)')
                ax.set_ylabel(r'Trajectory points')

                fig.savefig(imagesavedir+'/'+runname+'_'+treftimestring+'_'+var+'b_hist.png',dpi=200,bbox_inches='tight')
            
                fig = plt.figure()
                ax = fig.add_subplot(111)
                #ax.hist(vara_hist,bins=N.linspace(-0.5,0.5,100.),normed=False,log=False,color='b',alpha=1.0,orientation='vertical')
                if(masked_tornado):
                    vara_tor_hist = wforceatref[wforce_varnames.index(var),:,~tormask].flatten()
                    if(numtortrajcs > 0):
                        ax.hist(vara_tor_hist,bins=N.linspace(-0.5,0.5,100.),normed=False,log=False,color='r',alpha=0.75,orientation='vertical')
                #ax.set_ylim(1.0,600.0)
                #ax.set_yscale('log')
                ax.set_xlim(-0.5,0.5)
                ax.set_xlabel(var+r'$acceleration (m s^{-2})$ (after reference time)')
                ax.set_ylabel(r'Trajectory points')

                fig.savefig(imagesavedir+'/'+runname+'_'+treftimestring+'_'+var+'a_hist.png',dpi=200,bbox_inches='tight')
            
        if(plot_DSD):   # Histograms of DSD quantities
            print "Creating histograms of DSD quantities"
                            
            qr_hist = qrbtref.flatten()
            qc_hist = qcbtref.flatten()
            qh_hist = qhbtref.flatten()
            evapqc_hist = evapqcbtref.flatten()
            evapqr_hist = evapqrbtref.flatten()
            meltqh_hist = meltqhbtref.flatten()
            Dmr_hist = Dmrbtref.flatten()
            Dmh_hist = Dmhbtref.flatten()
        
            fig = plt.figure()
            ax = fig.add_subplot(111)
            #ax.hist(varb_hist,bins=N.linspace(-0.5,0.5,100.),normed=False,log=False,color='b',alpha=1.0,orientation='vertical')
            if(masked_tornado):
                qr_tor_hist = qrbtref[:,~tormask].flatten()
                if(numtortrajcs > 0):
                    ax.hist(qr_tor_hist,bins=N.linspace(-5.,5.,100.),normed=False,log=False,color='r',alpha=0.75,orientation='vertical')
            #ax.set_ylim(1.0,600.0)
            #ax.set_yscale('log')
            ax.set_xlim(-5.,5.)
            ax.set_xlabel(r'$q_{r} (g kg^{-1})$ (before reference time)')
            ax.set_ylabel(r'Trajectory points')

            fig.savefig(imagesavedir+'/'+runname+'_'+treftimestring+'_qr_hist.png',dpi=200,bbox_inches='tight')
        
            fig = plt.figure()
            ax = fig.add_subplot(111)
            #ax.hist(varb_hist,bins=N.linspace(-0.5,0.5,100.),normed=False,log=False,color='b',alpha=1.0,orientation='vertical')
            if(masked_tornado):
                qc_tor_hist = qcbtref[:,~tormask].flatten()
                if(numtortrajcs > 0):
                    ax.hist(qc_tor_hist,bins=N.linspace(-5.,5.,100.),normed=False,log=False,color='r',alpha=0.75,orientation='vertical')
            #ax.set_ylim(1.0,600.0)
            #ax.set_yscale('log')
            ax.set_xlim(-5.,5.)
            ax.set_xlabel(r'$q_{c} (g kg^{-1})$ (before reference time)')
            ax.set_ylabel(r'Trajectory points')

            fig.savefig(imagesavedir+'/'+runname+'_'+treftimestring+'_qc_hist.png',dpi=200,bbox_inches='tight')
        
            fig = plt.figure()
            ax = fig.add_subplot(111)
            #ax.hist(varb_hist,bins=N.linspace(-0.5,0.5,100.),normed=False,log=False,color='b',alpha=1.0,orientation='vertical')
            if(masked_tornado):
                qh_tor_hist = qhbtref[:,~tormask].flatten()
                if(numtortrajcs > 0):
                    ax.hist(qh_tor_hist,bins=N.linspace(-5.,5.,100.),normed=False,log=False,color='r',alpha=0.75,orientation='vertical')
            #ax.set_ylim(1.0,600.0)
            #ax.set_yscale('log')
            ax.set_xlim(-5.,5.)
            ax.set_xlabel(r'$q_{h} (g kg^{-1})$ (before reference time)')
            ax.set_ylabel(r'Trajectory points')

            fig.savefig(imagesavedir+'/'+runname+'_'+treftimestring+'_qh_hist.png',dpi=200,bbox_inches='tight')
        
            fig = plt.figure()
            ax = fig.add_subplot(111)
            #ax.hist(varb_hist,bins=N.linspace(-0.5,0.5,100.),normed=False,log=False,color='b',alpha=1.0,orientation='vertical')
            if(masked_tornado):
                evapqr_tor_hist = evapqrbtref[:,~tormask].flatten()
                if(numtortrajcs > 0):
                    ax.hist(evapqr_tor_hist,bins=N.linspace(-0.01,0.01,100.),normed=False,log=False,color='r',alpha=0.75,orientation='vertical')
            #ax.set_ylim(1.0,600.0)
            #ax.set_yscale('log')
            ax.set_xlim(-0.01,0.01)
            ax.set_xlabel(r'$evapq_{r} (g kg^{-1} s^{-1})$ (before reference time)')
            ax.set_ylabel(r'Trajectory points')

            fig.savefig(imagesavedir+'/'+runname+'_'+treftimestring+'_evapqr_hist.png',dpi=200,bbox_inches='tight')
        
            fig = plt.figure()
            ax = fig.add_subplot(111)
            #ax.hist(varb_hist,bins=N.linspace(-0.5,0.5,100.),normed=False,log=False,color='b',alpha=1.0,orientation='vertical')
            if(masked_tornado):
                evapqc_tor_hist = evapqcbtref[:,~tormask].flatten()
                if(numtortrajcs > 0):
                    ax.hist(evapqc_tor_hist,bins=N.linspace(-0.01,0.01,100.),normed=False,log=False,color='r',alpha=0.75,orientation='vertical')
            #ax.set_ylim(1.0,600.0)
            #ax.set_yscale('log')
            ax.set_xlim(-0.01,0.01)
            ax.set_xlabel(r'$evapq_{c} (g kg^{-1} s^{-1})$ (before reference time)')
            ax.set_ylabel(r'Trajectory points')

            fig.savefig(imagesavedir+'/'+runname+'_'+treftimestring+'_evapqc_hist.png',dpi=200,bbox_inches='tight')
        
            fig = plt.figure()
            ax = fig.add_subplot(111)
            #ax.hist(varb_hist,bins=N.linspace(-0.5,0.5,100.),normed=False,log=False,color='b',alpha=1.0,orientation='vertical')
            if(masked_tornado):
                meltqh_tor_hist = meltqhbtref[:,~tormask].flatten()
                if(numtortrajcs > 0):
                    ax.hist(meltqh_tor_hist,bins=N.linspace(-0.01,0.01,100.),normed=False,log=False,color='r',alpha=0.75,orientation='vertical')
            #ax.set_ylim(1.0,600.0)
            #ax.set_yscale('log')
            ax.set_xlim(-0.01,0.01)
            ax.set_xlabel(r'$meltq_{h} (g kg^{-1} s^{-1})$ (before reference time)')
            ax.set_ylabel(r'Trajectory points')

            fig.savefig(imagesavedir+'/'+runname+'_'+treftimestring+'_meltqh_hist.png',dpi=200,bbox_inches='tight')
        
            fig = plt.figure()
            ax = fig.add_subplot(111)
            #ax.hist(varb_hist,bins=N.linspace(-0.5,0.5,100.),normed=False,log=False,color='b',alpha=1.0,orientation='vertical')
            if(masked_tornado):
                Dmr_tor_hist = Dmrbtref[:,~tormask].flatten().compressed()
                if(numtortrajcs > 0):
                    ax.hist(Dmr_tor_hist,bins=N.linspace(0.0,8.0,100.),normed=False,log=False,color='r',alpha=0.75,orientation='vertical')
            #ax.set_ylim(1.0,600.0)
            #ax.set_yscale('log')
            ax.set_xlim(0.0,8.0)
            ax.set_xlabel(r'$D_{mr} (mm)$ (before reference time)')
            ax.set_ylabel(r'Trajectory points')

            fig.savefig(imagesavedir+'/'+runname+'_'+treftimestring+'_Dmr_hist.png',dpi=200,bbox_inches='tight')
    
            fig = plt.figure()
            ax = fig.add_subplot(111)
            #ax.hist(varb_hist,bins=N.linspace(-0.5,0.5,100.),normed=False,log=False,color='b',alpha=1.0,orientation='vertical')
            if(masked_tornado):
                Dmh_tor_hist = Dmhbtref[:,~tormask].flatten().compressed()
                if(numtortrajcs > 0):
                    ax.hist(Dmh_tor_hist,bins=N.linspace(0.0,20.0,100.),normed=False,log=False,color='r',alpha=0.75,orientation='vertical')
            #ax.set_ylim(1.0,600.0)
            #ax.set_yscale('log')
            ax.set_xlim(0.0,20.0)
            ax.set_xlabel(r'$D_{mh} (mm)$ (before reference time)')
            ax.set_ylabel(r'Trajectory points')

            fig.savefig(imagesavedir+'/'+runname+'_'+treftimestring+'_Dmr_hist.png',dpi=200,bbox_inches='tight')
    
    if(plot3D):
        if(masked_tornado):
            vorintrp = varintrp[varnames.index('vor'),...]
        
        xmin3d = xvals[trefindex,:].min()
        xmax3d = xvals[trefindex,:].max()
        ymin3d = yvals[trefindex,:].min()
        ymax3d = yvals[trefindex,:].max()

        image3ddir_sub = image3ddir+runname+'_'+treftimestring+'/'
        
        if (not os.path.exists(image3ddir_sub)):
            os.mkdir(image3ddir_sub)
        
        # Create a scatter plot of trajectory points vs. time, colored by a given field
        
        if(plot_var):
            for i,var,varlim in zip(xrange(len(varnames)),varnames,varlims3D):
                fig = plt.figure()
                ax = fig.gca(projection='3d')
                ax.view_init(elev=30,azim=300) # From the southeast
                
                print "Min/max of "+var+" for trajectories for reference time "+treftimestring
                print varintrp[i,...].min(),varintrp[i,...].max()
                                
                #norm = colors.Normalize(varintrp[i,...].min(),varintrp[i,...].max())
                norm = colors.Normalize(varlim[0],varlim[1])
                
                #scalarMap = cm.ScalarMappable(norm=norm,cmap=cm.RdBu_r)
                
                if(masked_tornado and usetorreftime and plot_torpoints):
                    xvalstor3d = N.ma.masked_array(xvals[trefindex:,:],mask=torrefmask)
                    yvalstor3d = N.ma.masked_array(yvals[trefindex:,:],mask=torrefmask)
                    zvalstor3d = N.ma.masked_array(zvals[trefindex:,:],mask=torrefmask)
                    vartor3d = N.ma.masked_array(varintrp[i,trefindex:,:],mask=torrefmask)
                    sc = ax.scatter(xvalstor3d,yvalstor3d,zvalstor3d,c=vartor3d,cmap=arps.blue_red1,norm=norm,edgecolors='none',s=2)
                    plt.colorbar(sc)
                    ax.set_xlim3d(xmin3d,xmax3d)
                    ax.set_ylim3d(ymin3d,ymax3d)
                    ax.set_zlim3d(0.0,5000.0)
                    plt.savefig(image3ddir_sub+runname+'_'+var+'_3dtraj_torpoints.png',dpi=300)
                else:
                    for ttraj in xrange(ntimes):
                        ax.cla()
                        if(masked_tornado):
                            xvalstor3d = xvals[ttraj,~tormask]
                            yvalstor3d = yvals[ttraj,~tormask]
                            zvalstor3d = zvals[ttraj,~tormask]
                            vartor3d = varintrp[i,ttraj,~tormask]
                            sc = ax.scatter(xvalstor3d,yvalstor3d,zvalstor3d,c=vartor3d,cmap=cm.RdBu_r,norm=norm,edgecolors='none',s=2,marker='s')
                        
                            xvalsnontor3d = xvals[ttraj,tormask]
                            yvalsnontor3d = yvals[ttraj,tormask]
                            zvalsnontor3d = zvals[ttraj,tormask]
                            varnontor3d = varintrp[i,ttraj,tormask]
                            sc = ax.scatter(xvalsnontor3d,yvalsnontor3d,zvalsnontor3d,c=varnontor3d,cmap=cm.RdBu_r,norm=norm,edgecolors='none',s=2)
                       
                        else:
                            sc = ax.scatter(xvals[ttraj,:],yvals[ttraj,:],zvals[ttraj,:],c=varintrp[i,ttraj,:],cmap=cm.RdBu_r,norm=norm,edgecolors='none',s=2)
                        if(ttraj == 0):    
                            plt.colorbar(sc)
                        ax.set_xlim3d(xmin3d,xmax3d)
                        ax.set_ylim3d(ymin3d,ymax3d)
                        ax.set_zlim3d(0.0,5000.0)
                        plt.savefig(image3ddir_sub+runname+'_'+var+'_3dtraj_'+'%02d' % ttraj+'.png',dpi=300)
        
        if(plot_wforce):
            for i,var in enumerate(wforce_varnames):
                fig = plt.figure()
                ax = fig.gca(projection='3d')
                ax.view_init(elev=30,azim=300) # From the southeast
                
                print "Min/max of "+var+" for trajectories for reference time "+treftimestring
                print wforceintrp[i,...].min(),wforceintrp[i,...].max()
                
                norm = colors.Normalize(wforce_lim3d[0],wforce_lim3d[1])
                #scalarMap = cm.ScalarMappable(norm=norm,cmap=cm.RdBu_r)
                
                if(masked_tornado and usetorreftime and plot_torpoints):
                    xvalstor3d = N.ma.masked_array(xvals[trefindex:,:],mask=torrefmask)
                    yvalstor3d = N.ma.masked_array(yvals[trefindex:,:],mask=torrefmask)
                    zvalstor3d = N.ma.masked_array(zvals[trefindex:,:],mask=torrefmask)
                    wforcetor3d = N.ma.masked_array(wforceintrp[i,trefindex:,:],mask=torrefmask)
                    sc = ax.scatter(xvalstor3d,yvalstor3d,zvalstor3d,c=wforcetor3d,cmap=arps.blue_red1,norm=norm,edgecolors='none',s=2)
                    plt.colorbar(sc)
                    ax.set_xlim3d(xmin3d,xmax3d)
                    ax.set_ylim3d(ymin3d,ymax3d)
                    ax.set_zlim3d(0.0,5000.0)
                    plt.savefig(image3ddir_sub+runname+'_'+var+'_3dtraj_torpoints.png',dpi=300)
                else:
                    for ttraj in xrange(ntimes):
                        ax.cla()
                        if(masked_tornado):
                            xvalstor3d = xvals[ttraj,~tormask]
                            yvalstor3d = yvals[ttraj,~tormask]
                            zvalstor3d = zvals[ttraj,~tormask]
                            vartor3d = wforceintrp[i,ttraj,~tormask]
                            sc = ax.scatter(xvalstor3d,yvalstor3d,zvalstor3d,c=vartor3d,cmap=cm.RdBu_r,norm=norm,edgecolors='none',s=2,marker='s')
                        
                            xvalsnontor3d = xvals[ttraj,tormask]
                            yvalsnontor3d = yvals[ttraj,tormask]
                            zvalsnontor3d = zvals[ttraj,tormask]
                            varnontor3d = wforceintrp[i,ttraj,tormask]
                            sc = ax.scatter(xvalsnontor3d,yvalsnontor3d,zvalsnontor3d,c=varnontor3d,cmap=cm.RdBu_r,norm=norm,edgecolors='none',s=2)
                       
                        else:
                            sc = ax.scatter(xvals[ttraj,:],yvals[ttraj,:],zvals[ttraj,:],c=wforceintrp[i,ttraj,:],cmap=cm.RdBu_r,norm=norm,edgecolors='none',s=2)
                        if(ttraj == 0):    
                            plt.colorbar(sc)
                        ax.set_xlim3d(xmin3d,xmax3d)
                        ax.set_ylim3d(ymin3d,ymax3d)
                        ax.set_zlim3d(0.0,5000.0)
                        plt.savefig(image3ddir_sub+runname+'_'+var+'_3dtraj_'+'%02d' % ttraj+'.png',dpi=300)

        if(plot_DSD):
            for i,var,varlim in zip(xrange(len(DSD_varnames)),DSD_varnames,DSDlims3d):
                fig = plt.figure()
                ax = fig.gca(projection='3d')
                ax.view_init(elev=30,azim=300) # From the southeast
                
                print "Min/max of "+var+" for trajectories for reference time "+treftimestring
                print DSDintrp[i,...].min(),DSDintrp[i,...].max()
                
                norm = colors.Normalize(varlim[0],varlim[1])
                #scalarMap = cm.ScalarMappable(norm=norm,cmap=cm.RdBu_r)
                if(masked_tornado and usetorreftime and plot_torpoints):
                    xvalstor3d = N.ma.masked_array(xvals[trefindex:,:],mask=torrefmask)
                    yvalstor3d = N.ma.masked_array(yvals[trefindex:,:],mask=torrefmask)
                    zvalstor3d = N.ma.masked_array(zvals[trefindex:,:],mask=torrefmask)
                    DSDtor3d = N.ma.masked_array(DSDintrp[i,trefindex:,:],mask=torrefmask)
                    sc = ax.scatter(xvalstor3d,yvalstor3d,zvalstor3d,c=DSDtor3d,cmap=arps.blue_red1,norm=norm,edgecolors='none',s=2)
                    plt.colorbar(sc)
                    ax.set_xlim3d(xmin3d,xmax3d)
                    ax.set_ylim3d(ymin3d,ymax3d)
                    ax.set_zlim3d(0.0,5000.0)
                    plt.savefig(image3ddir_sub+runname+'_'+var+'_3dtraj_torpoints.png',dpi=300)
                else:
                    for ttraj in xrange(ntimes):
                        ax.cla()
                        if(masked_tornado):
                            xvalstor3d = xvals[ttraj,~tormask]
                            yvalstor3d = yvals[ttraj,~tormask]
                            zvalstor3d = zvals[ttraj,~tormask]
                            vartor3d = DSDintrp[i,ttraj,~tormask]
                            sc = ax.scatter(xvalstor3d,yvalstor3d,zvalstor3d,c=vartor3d,cmap=cm.RdBu_r,norm=norm,edgecolors='none',s=2,marker='s')
                        
                            xvalsnontor3d = xvals[ttraj,tormask]
                            yvalsnontor3d = yvals[ttraj,tormask]
                            zvalsnontor3d = zvals[ttraj,tormask]
                            varnontor3d = DSDintrp[i,ttraj,tormask]
                            sc = ax.scatter(xvalsnontor3d,yvalsnontor3d,zvalsnontor3d,c=varnontor3d,cmap=cm.RdBu_r,norm=norm,edgecolors='none',s=2)
                       
                        else:
                            sc = ax.scatter(xvals[ttraj,:],yvals[ttraj,:],zvals[ttraj,:],c=wforceintrp[i,ttraj,:],cmap=cm.RdBu_r,norm=norm,edgecolors='none',s=2)
                        if(ttraj == 0):    
                            plt.colorbar(sc)
                        ax.set_xlim3d(xmin3d,xmax3d)
                        ax.set_ylim3d(ymin3d,ymax3d)
                        ax.set_zlim3d(0.0,5000.0)
                        plt.savefig(image3ddir_sub+runname+'_'+var+'_3dtraj_'+'%02d' % ttraj+'.png',dpi=300)

    if(plot_scatter and not loadTS):   # Similar to the 3D plot, but plotting quantities vs. height outside and inside TLV
        
        imagescatterdir = imagesavedir+runname+'_'+treftimestring+'/'
        if (not os.path.exists(imagescatterdir)):
            os.mkdir(imagescatterdir)
        
        if(masked_tornado and usetorreftime):
            zvalstor = N.ma.masked_array(zvals_agl,mask=torrefmaskfull)
            zvalstor = zvalstor[:,~tormask]
            zvalsbtor = N.ma.masked_array(zvals_agl,mask=~torrefmaskfull)
            zvalsbtor = zvalsbtor[:,~tormask]
            
            # Flatten the height array and then digitize into
            # 100-m bins.  Using the index information we can compute the average of the variable
            # array for points in each of the height bins.
            
            zvalsinitbtor = zvalsbtor[0,:].flatten().compressed()
            zvalsfintor = zvalstor[-1,:].flatten().compressed()
            
            zvalstor = zvalstor.flatten()
            zvalsbtor = zvalsbtor.flatten()
            if(numtortrajcs > 0):
                height_indices = N.digitize(zvalstor,height_bins)
                heightb_indices = N.digitize(zvalsbtor,height_bins)
                
            
            
        if(plot_var):
            if(numtortrajcs > 0):
                # Count number of initial/final trajectory points in each height bin
                npoints,dummy = N.histogram(zvalsfintor,bins=height_bins)
                npointsb,dummy = N.histogram(zvalsinitbtor,bins=height_bins)
                npoints_TS[t,:] = npoints
                npointsb_TS[t,:] = npointsb
            else:
                npoints_TS[t,:] = N.nan
                npointsb_TS[t,:] = N.nan
                
            for i,var,varlim in zip(xrange(len(varnames)),varnames,varlims3D):
                figtor = plt.figure()
                axtor = figtor.add_subplot(111)
                plt.title(var+' vs. height AGL within TLV for reference time '+treftimestring+' s')
                figbtor = plt.figure()
                axbtor = figbtor.add_subplot(111)
                plt.title(var+' vs. height AGL before TLV for reference time '+treftimestring+' s')
                print "Min/max of "+var+" for trajectories for reference time "+treftimestring
                print varintrp[i,...].min(),varintrp[i,...].max()
                                
                #norm = colors.Normalize(varintrp[i,...].min(),varintrp[i,...].max())
                norm = colors.Normalize(varlim[0],varlim[1])
                
                #scalarMap = cm.ScalarMappable(norm=norm,cmap=cm.RdBu_r)
                
                if(masked_tornado and usetorreftime):
                    vartor = N.ma.masked_array(varintrp[i,...],mask=torrefmaskfull)                    
                    vartor = vartor[:,~tormask]                 
                    varbtor = N.ma.masked_array(varintrp[i,...],mask=~torrefmaskfull)               
                    varbtor = varbtor[:,~tormask]
                    
                    # Flatten the variable arrays and compute the average of the variable
                    # array for points in each of the height bins.
                
                    vartor = vartor.flatten()
                    varbtor = varbtor.flatten()    
                    
                    # Now compute the average var in each 100-m height bin
                    # Thanks to http://stackoverflow.com/questions/6163334/binning-data-in-python-with-scipy-numpy 
                    if(numtortrajcs > 0):
                        var_avg = N.array([vartor[height_indices == j].mean() for j in range(1, len(height_bins))])
                        varb_avg = N.array([varbtor[heightb_indices == j].mean() for j in range(1, len(height_bins))])
                        sc = axtor.scatter(vartor,zvalstor,c='k',edgecolors='none',s=2)
                        axtor.plot(var_avg,height_middle,c='r',lw=2)
                        axtor.axvline(x=0.)
                        axtor.set_xlim(varlim[0],varlim[1])
                        axtor.set_ylim(0.0,5000.0)
                        figtor.savefig(imagescatterdir+runname+'_'+var+'_'+treftimestring+'_height_scatter_torpoints.png',dpi=300)
                        plt.close(figtor)
                    
                        sc = axbtor.scatter(varbtor,zvalsbtor,c='k',edgecolors='none',s=2)
                        axbtor.plot(varb_avg,height_middle,c='r',lw=2)
                        axbtor.axvline(x=0.)
                        axbtor.set_xlim(varlim[0],varlim[1])
                        axbtor.set_ylim(0.0,2000.0)
                        figbtor.savefig(imagescatterdir+runname+'_'+var+'_'+treftimestring+'_height_scatter_btorpoints.png',dpi=300)
                        plt.close(figbtor)
                    else:
                        var_avg = N.nan
                        varb_avg = N.nan
                    
                    # Save averages as function of height for a later timeseries plot
                    var_avg_TS[i,t,:] = var_avg
                    varb_avg_TS[i,t,:] = varb_avg
        if(plot_wforce):
            for i,var in zip(xrange(len(wforce_varnames)),wforce_varnames):
                figtor = plt.figure()
                axtor = figtor.add_subplot(111)
                plt.title(var+' vs. height AGL within TLV for reference time '+treftimestring+' s')
                figbtor = plt.figure()
                axbtor = figbtor.add_subplot(111)
                plt.title(var+' vs. height AGL before TLV for reference time '+treftimestring+' s')
                #figbtor = plt.figure()
                #axbtor = figbtor.add_subplot(111) 
                print "Min/max of "+var+" for trajectories for reference time "+treftimestring
                print wforceintrp[i,...].min(),wforceintrp[i,...].max()
                                
                #norm = colors.Normalize(varintrp[i,...].min(),varintrp[i,...].max())
                norm = colors.Normalize(wforce_lim3d[0],wforce_lim3d[1])
                
                #scalarMap = cm.ScalarMappable(norm=norm,cmap=cm.RdBu_r)
                
                if(masked_tornado and usetorreftime):
                    
                    wforcetor = N.ma.masked_array(wforceintrp[i,...],mask=torrefmaskfull)                   
                    wforcetor = wforcetor[:,~tormask]                    
                    wforcebtor = N.ma.masked_array(wforceintrp[i,...],mask=~torrefmaskfull)                    
                    wforcebtor = wforcebtor[:,~tormask]
                    
                    # Flatten the wforcing arrays and compute the average of the wforcing
                    # array for points in each of the height bins.
                
                    wforcetor = wforcetor.flatten()
                    wforcebtor = wforcebtor.flatten()    
                    
                    # Now compute the average wforcing in each 100-m height bin
                    # Thanks to http://stackoverflow.com/questions/6163334/binning-data-in-python-with-scipy-numpy 
                    if(numtortrajcs > 0):
                        wforce_avg = N.array([wforcetor[height_indices == j].mean() for j in range(1, len(height_bins))])
                        wforceb_avg = N.array([wforcebtor[heightb_indices == j].mean() for j in range(1, len(height_bins))])
                        sc = axtor.scatter(wforcetor,zvalstor,c='k',edgecolors='none',s=2)
                        axtor.plot(wforce_avg,height_middle,c='r',lw=2)
                        axtor.axvline(x=0.)
                        axtor.set_xlim(wforce_lim3d[0],wforce_lim3d[1])
                        axtor.set_ylim(0.0,5000.0)
                        figtor.savefig(imagescatterdir+runname+'_'+var+'_'+treftimestring+'_height_scatter_torpoints.png',dpi=300)
                        plt.close(figtor)
                        
                        sc = axbtor.scatter(wforcebtor,zvalsbtor,c='k',edgecolors='none',s=2)
                        axbtor.plot(wforceb_avg,height_middle,c='r',lw=2)
                        axbtor.axvline(x=0.)
                        axbtor.set_xlim(wforce_lim3d[0],wforce_lim3d[1])
                        axbtor.set_ylim(0.0,2000.0)
                        figbtor.savefig(imagescatterdir+runname+'_'+var+'_'+treftimestring+'_height_scatter_btorpoints.png',dpi=300)
                        plt.close(figbtor)
                    else:
                        wforce_avg = N.nan
                        wforceb_avg = N.nan
                    
                    # Save averages as function of height for a later timeseries plot
                    wforce_avg_TS[i,t,:] = wforce_avg
                    wforceb_avg_TS[i,t,:] = wforceb_avg
        if(plot_DSD):
            for i,var,varlim in zip(xrange(len(DSD_varnames)),DSD_varnames,DSDlims3D):
                figtor = plt.figure()
                axtor = figtor.add_subplot(111)
                plt.title(var+' vs. height AGL within TLV for reference time '+treftimestring+' s')
                figbtor = plt.figure()
                axbtor = figbtor.add_subplot(111) 
                plt.title(var+' vs. height AGL before TLV for reference time '+treftimestring+' s')
                print "Min/max of "+var+" for trajectories for reference time "+treftimestring
                print DSDintrp[i,...].min(),DSDintrp[i,...].max()
                                
                #norm = colors.Normalize(varintrp[i,...].min(),varintrp[i,...].max())
                norm = colors.Normalize(varlim[0],varlim[1])
                
                #scalarMap = cm.ScalarMappable(norm=norm,cmap=cm.RdBu_r)
                
                if(masked_tornado and usetorreftime):
                    
                    DSDtor = N.ma.masked_array(DSDintrp[i,...],mask=torrefmaskfull)
                    DSDtor = DSDtor[:,~tormask]
                    DSDbtor = N.ma.masked_array(DSDintrp[i,...],mask=~torrefmaskfull)
                    DSDbtor = DSDbtor[:,~tormask]
                    
                    # Flatten the DSD arrays and compute the average of the DSD
                    # array for points in each of the height bins.
                
                    DSDtor = DSDtor.flatten()
                    DSDbtor = DSDbtor.flatten()    
                    
                    # Now compute the average DSD variable in each 100-m height bin
                    # Thanks to http://stackoverflow.com/questions/6163334/binning-data-in-python-with-scipy-numpy 
                    if(numtortrajcs > 0):
                        DSD_avg = N.array([DSDtor[height_indices == j].mean() for j in range(1, len(height_bins))])
                        DSDb_avg = N.array([DSDbtor[heightb_indices == j].mean() for j in range(1, len(height_bins))])
                        sc = axtor.scatter(DSDtor,zvalstor,c='k',edgecolors='none',s=2)
                        axtor.plot(DSD_avg,height_middle,c='r',lw=2)
                        axtor.axvline(x=0.)
                        axtor.set_xlim(varlim[0],varlim[1])
                        axtor.set_ylim(0.0,5000.0)
                        figtor.savefig(imagescatterdir+runname+'_'+var+'_'+treftimestring+'_height_scatter_torpoints.png',dpi=300)
                        plt.close(figtor)
                    
                        sc = axbtor.scatter(DSDbtor,zvalsbtor,c='k',edgecolors='none',s=2)
                        axbtor.plot(DSDb_avg,height_middle,c='r',lw=2)
                        axbtor.axvline(x=0.)
                        axbtor.set_xlim(varlim[0],varlim[1])
                        axbtor.set_ylim(0.0,2000.0)
                        figbtor.savefig(imagescatterdir+runname+'_'+var+'_'+treftimestring+'_height_scatter_btorpoints.png',dpi=300)
                        plt.close(figbtor)
                    else:
                        DSD_avg = N.nan
                        DSDb_avg = N.nan
                    # Save averages as function of height for a later timeseries plot
                    DSD_avg_TS[i,t,:] = DSD_avg
                    DSDb_avg_TS[i,t,:] = DSDb_avg
# Plot averages as function of height (from scatterplots) as a timeseries, if desired

if(plot_scatter and plot_TS):
    if(plot_var):
        if(saveTS):
            npzfilename = npzsavedir+'/'+runname+'_varTH.npz'
            print "Saving time-height data for averaged trajectory variables"
            savevars = {}
            savevars['npointsTH'] = npoints_TS
            savevars['npointsbTH'] = npointsb_TS
            savevars['varTH'] = var_avg_TS
            savevars['varbTH'] = varb_avg_TS
            N.savez(npzfilename,**savevars)
        elif(loadTS):
            npzfilename = npzsavedir+'/'+runname+'_varTH.npz'
            print "Loading time-height data for averaged trajectory variables"
            varavgfile = N.load(npzfilename)
            npoints_TS = varavgfile['npointsTH']
            npointsb_TS = varavgfile['npointsbTH']
            var_avg_TS = varavgfile['varTH']
            varb_avg_TS = varavgfile['varbTH']
            #print var_avg_TS
            
        # Now mask out invalid entries (i.e. for times where no TLV is present)
        npoints_TS = N.ma.masked_invalid(npoints_TS)
        npointsb_TS = N.ma.masked_invalid(npointsb_TS)
        var_avg_TS = N.ma.masked_invalid(var_avg_TS)
        varb_avg_TS = N.ma.masked_invalid(varb_avg_TS)
        
        npoints_heightsums = npoints_TS.sum(axis=1)[:,N.newaxis]
        npointsb_heightsums = npointsb_TS.sum(axis=1)[:,N.newaxis]
        nfrac_TS = npoints_TS/npoints_heightsums
        nfracb_TS = npointsb_TS/npointsb_heightsums
        
        figtor = plt.figure()
        axtor = figtor.add_subplot(111)
        cmap = arps.blue_red1
        plt.title('Fraction of terminal trajectory points vs. height within TLV')
        plot = axtor.pcolormesh(tref_list,height_middle,nfrac_TS.T,vmin=0.0,vmax=0.25,cmap=cmap,edgecolors='None',antialiased=False,rasterized=True)
        divider = make_axes_locatable(axtor)
        cax = divider.append_axes("top",size="5%",pad=0.05)
        cax.xaxis.set_ticks_position('top')
        plt.colorbar(plot,orientation='horizontal',cax=cax)
        axtor.set_xlim(tref_list[0],tref_list[-1])
        axtor.set_ylim(0.0,5000.0)
        figtor.savefig(imagesavedir+runname+'npoints_TH_torpoints.png',dpi=300)
        
        figbtor = plt.figure()
        axbtor = figbtor.add_subplot(111)
        cmap = arps.blue_red1
        plt.title('Fraction of initial trajectory points vs. height AGL before TLV')
        plot = axbtor.pcolormesh(tref_list,height_middle,nfracb_TS.T,vmin=0.0,vmax=0.25,cmap=cmap,edgecolors='None',antialiased=False,rasterized=True)
        divider = make_axes_locatable(axbtor)
        cax = divider.append_axes("top",size="5%",pad=0.05)
        cax.xaxis.set_ticks_position('top')
        plt.colorbar(plot,orientation='horizontal',cax=cax)
        axbtor.set_xlim(tref_list[0],tref_list[-1])
        axbtor.set_ylim(0.0,5000.0)
        figbtor.savefig(imagesavedir+runname+'npoints_TH_btorpoints.png',dpi=300)
        
        for i,var,varlim in zip(xrange(len(varnames)),varnames,varlims3D):
            figtor = plt.figure()
            axtor = figtor.add_subplot(111)
            cmap = arps.blue_red1
            plt.title('Mean '+var+' vs. height AGL within TLV')
            plot = axtor.pcolormesh(tref_list,height_middle,var_avg_TS[i,...].T,vmin=varlim[0],vmax=varlim[1],cmap=cmap,edgecolors='None',antialiased=False,rasterized=True)
            divider = make_axes_locatable(axtor)
            cax = divider.append_axes("top",size="5%",pad=0.05)
            cax.xaxis.set_ticks_position('top')
            plt.colorbar(plot,orientation='horizontal',cax=cax)
            axtor.set_xlim(tref_list[0],tref_list[-1])
            axtor.set_ylim(0.0,5000.0)
            figtor.savefig(imagesavedir+runname+'_'+var+'_TH_torpoints.png',dpi=300)
            
            figbtor = plt.figure()
            axbtor = figbtor.add_subplot(111)
            cmap = arps.blue_red1
            plt.title('Mean '+var+' vs. height AGL before TLV')
            plot = axbtor.pcolormesh(tref_list,height_middle,varb_avg_TS[i,...].T,vmin=varlim[0],vmax=varlim[1],cmap=cmap,edgecolors='None',antialiased=False,rasterized=True)
            divider = make_axes_locatable(axbtor)
            cax = divider.append_axes("top",size="5%",pad=0.05)
            cax.xaxis.set_ticks_position('top')
            plt.colorbar(plot,orientation='horizontal',cax=cax)
            axbtor.set_xlim(tref_list[0],tref_list[-1])
            axbtor.set_ylim(0.0,5000.0)
            figbtor.savefig(imagesavedir+runname+'_'+var+'_TH_btorpoints.png',dpi=300)
            #figtor.savefig('/Users/ddawson/Dropbox/may0399_real_data_paper/temp/'+runname+'_'+var+'_'+treftimestring+'_height_TS_torpoints.png',dpi=300)
    if(plot_wforce):
        if(saveTS):
            npzfilename = npzsavedir+'/'+runname+'_wforceTH.npz'
            print "Saving time-height data for averaged trajectory w-forcing variables"
            savevars = {}
            savevars['wforceTH'] = wforce_avg_TS
            savevars['wforcebTH'] = wforceb_avg_TS
            N.savez(npzfilename,**savevars)
        elif(loadTS):
            npzfilename = npzsavedir+'/'+runname+'_wforceTH.npz'
            print "Loading time-height data for averaged trajectory w-forcing variables"
            varavgfile = N.load(npzfilename)
        
            wforce_avg_TS = varavgfile['wforceTH']
            wforceb_avg_TS = varavgfile['wforcebTH']
        # Now mask out invalid entries (i.e. for times where no TLV is present)
        wforce_avg_TS = N.ma.masked_invalid(wforce_avg_TS)
        wforceb_avg_TS = N.ma.masked_invalid(wforceb_avg_TS)
        for i,var in zip(xrange(len(wforce_varnames)),wforce_varnames):
            figtor = plt.figure()
            axtor = figtor.add_subplot(111)
            cmap = arps.blue_red1
            plt.title('Mean '+var+' vs. height AGL within TLV')
            plot = axtor.pcolormesh(tref_list,height_middle,wforce_avg_TS[i,...].T,vmin=wforce_lim3d[0],vmax=wforce_lim3d[1],cmap=cmap,edgecolors='None',antialiased=False,rasterized=True)
            divider = make_axes_locatable(axtor)
            cax = divider.append_axes("top",size="5%",pad=0.05)
            cax.xaxis.set_ticks_position('top')
            plt.colorbar(plot,orientation='horizontal',cax=cax)
            axtor.set_xlim(tref_list[0],tref_list[-1])
            axtor.set_ylim(0.0,5000.0)
            figtor.savefig(imagesavedir+runname+'_'+var+'_TH_torpoints.png',dpi=300)
            #figtor.savefig('/Users/ddawson/Dropbox/may0399_real_data_paper/temp/'+runname+'_'+var+'_'+treftimestring+'_height_TS_torpoints.png',dpi=300)
    
            figbtor = plt.figure()
            axbtor = figbtor.add_subplot(111)
            cmap = arps.blue_red1
            plt.title('Mean '+var+' vs. height AGL before TLV')
            plot = axbtor.pcolormesh(tref_list,height_middle,wforceb_avg_TS[i,...].T,vmin=wforce_lim3d[0],vmax=wforce_lim3d[1],cmap=cmap,edgecolors='None',antialiased=False,rasterized=True)
            divider = make_axes_locatable(axbtor)
            cax = divider.append_axes("top",size="5%",pad=0.05)
            cax.xaxis.set_ticks_position('top')
            plt.colorbar(plot,orientation='horizontal',cax=cax)
            axbtor.set_xlim(tref_list[0],tref_list[-1])
            axbtor.set_ylim(0.0,5000.0)
            figbtor.savefig(imagesavedir+runname+'_'+var+'_TH_btorpoints.png',dpi=300)
            #figtor.savefig('/Users/ddawson/Dropbox/may0399_real_data_paper/temp/'+runname+'_'+var+'_'+treftimestring+'_height_TS_torpoints.png',dpi=300)
    if(plot_DSD):
        if(saveTS):
            npzfilename = npzsavedir+'/'+runname+'_DSDTH.npz'
            print "Saving time-height data for averaged trajectory DSD variables"
            savevars = {}
            savevars['DSDTH'] = DSD_avg_TS
            savevars['DSDbTH'] = DSDb_avg_TS
            N.savez(npzfilename,**savevars)
        elif(loadTS):
            npzfilename = npzsavedir+'/'+runname+'_DSDTH.npz'
            print "Loading time-height data for averaged trajectory variables"
            varavgfile = N.load(npzfilename)
        
            DSD_avg_TS = varavgfile['DSDTH']
            DSDb_avg_TS = varavgfile['DSDbTH']
        # Now mask out invalid entries (i.e. for times where no TLV is present)
        DSD_avg_TS = N.ma.masked_invalid(DSD_avg_TS)
        DSDb_avg_TS = N.ma.masked_invalid(DSDb_avg_TS)
        for i,var,varlim in zip(xrange(len(DSD_varnames)),DSD_varnames,DSDlims3D):
            figtor = plt.figure()
            axtor = figtor.add_subplot(111)
            cmap = arps.blue_red1
            plt.title('Mean '+var+' vs. height AGL within TLV')
            plot = axtor.pcolormesh(tref_list,height_middle,DSD_avg_TS[i,...].T,vmin=varlim[0],vmax=varlim[1],cmap=cmap,edgecolors='None',antialiased=False,rasterized=True)
            divider = make_axes_locatable(axtor)
            cax = divider.append_axes("top",size="5%",pad=0.05)
            cax.xaxis.set_ticks_position('top')
            plt.colorbar(plot,orientation='horizontal',cax=cax)
            axtor.set_xlim(tref_list[0],tref_list[-1])
            axtor.set_ylim(0.0,5000.0)
            figtor.savefig(imagesavedir+runname+'_'+var+'_TH_torpoints.png',dpi=300)
            #figtor.savefig('/Users/ddawson/Dropbox/may0399_real_data_paper/temp/'+runname+'_'+var+'_'+treftimestring+'_height_TS_torpoints.png',dpi=300)
            
            figbtor = plt.figure()
            axbtor = figbtor.add_subplot(111)
            cmap = arps.blue_red1
            plt.title('Mean '+var+' vs. height AGL before TLV')
            plot = axbtor.pcolormesh(tref_list,height_middle,DSDb_avg_TS[i,...].T,vmin=varlim[0],vmax=varlim[1],cmap=cmap,edgecolors='None',antialiased=False,rasterized=True)
            divider = make_axes_locatable(axbtor)
            cax = divider.append_axes("top",size="5%",pad=0.05)
            cax.xaxis.set_ticks_position('top')
            plt.colorbar(plot,orientation='horizontal',cax=cax)
            axbtor.set_xlim(tref_list[0],tref_list[-1])
            axbtor.set_ylim(0.0,5000.0)
            figbtor.savefig(imagesavedir+runname+'_'+var+'_TH_btorpoints.png',dpi=300)
            #figtor.savefig('/Users/ddawson/Dropbox/may0399_real_data_paper/temp/'+runname+'_'+var+'_'+treftimestring+'_height_TS_torpoints.png',dpi=300)    

# Determine whether to save/load and plot timeseries

if(plot_TS and saveTS):
    if(plot_var):
        npzfilename = npzsavedir+'/'+runname+'_varTS.npz'
        print "Saving time series for averaged trajectory variables"
        savevars = {}
        savevars['zmaxbTS'] = zmaxbTS
        savevars['zmaxaTS'] = zmaxaTS
        savevars['varminbTS'] = varminbTS
        savevars['varminaTS'] = varminaTS
        savevars['varmaxbTS'] = varmaxbTS
        savevars['varmaxaTS'] = varmaxaTS
        savevars['varavgbTS'] = varavgbTS
        savevars['varavgaTS'] = varavgaTS
        savevars['numtortrajcsTS'] = numtortrajcsTS
        savevars['varminattorTS'] = varminattorTS
        savevars['varavgattorTS'] = varavgattorTS
        savevars['varmaxattorTS'] = varmaxattorTS
        savevars['varminbtorTS'] = varminbtorTS
        savevars['varmaxbtorTS'] = varmaxbtorTS
        savevars['varavgbtorTS'] = varavgbtorTS
        savevars['varminatorTS'] = varminatorTS
        savevars['varmaxatorTS'] = varmaxatorTS
        savevars['varavgatorTS'] = varavgatorTS
        
        N.savez(npzfilename,**savevars)
    
    if(plot_wforce):
        npzfilename = npzsavedir+'/'+runname+'_wforceTS.npz'
        print "Saving time series for averaged trajectory vertical velocity forcing"
        savevars = {}
        savevars['wforceminbTS'] = wforceminbTS
        savevars['wforceminaTS'] = wforceminaTS
        savevars['wforcemaxbTS'] = wforcemaxbTS
        savevars['wforcemaxaTS'] = wforcemaxaTS
        savevars['wforceavgbTS'] = wforceavgbTS
        savevars['wforceavgaTS'] = wforceavgaTS
        savevars['wforceminattorTS'] = wforceminattorTS
        savevars['wforceavgattorTS'] = wforceavgattorTS
        savevars['wforcemaxattorTS'] = wforcemaxattorTS
        savevars['wforceminattorTS'] = wforceminattorTS
        savevars['wforceminbtorTS'] = wforceminbtorTS
        savevars['wforceavgbtorTS'] = wforceavgbtorTS
        savevars['wforcemaxbtorTS'] = wforcemaxbtorTS
        savevars['wforceminatorTS'] = wforceminatorTS
        savevars['wforceavgatorTS'] = wforceavgatorTS
        savevars['wforcemaxatorTS'] = wforcemaxatorTS
        
        N.savez(npzfilename,**savevars)
        
    if(plot_DSD):
        npzfilename = npzsavedir+'/'+runname+'_DSDTS.npz'
        print "Saving time series for averaged trajectory DSD variables"
        savevars = {}
        savevars['qrmaxbTS'] = qrmaxbTS
        savevars['qcmaxbTS'] = qcmaxbTS
        savevars['qhmaxbTS'] = qhmaxbTS
        savevars['evapqcmaxbTS'] = evapqcmaxbTS
        savevars['evapqrmaxbTS'] = evapqrmaxbTS
        savevars['meltqhmaxbTS'] = meltqhmaxbTS
        savevars['DmrmaxbTS'] = DmrmaxbTS
        savevars['DmhmaxbTS'] = DmhmaxbTS
        savevars['alpharmaxbTS'] = alpharmaxbTS
        savevars['alphahmaxbTS'] = alphahmaxbTS
        savevars['qravgbTS'] = qravgbTS
        savevars['qcavgbTS'] = qcavgbTS
        savevars['qhavgbTS'] = qhavgbTS
        savevars['evapqcavgbTS'] = evapqcavgbTS
        savevars['evapqravgbTS'] = evapqravgbTS
        savevars['meltqhavgbTS'] = meltqhavgbTS
        savevars['DmravgbTS'] = DmravgbTS
        savevars['DmhavgbTS'] = DmhavgbTS
        savevars['alpharavgbTS'] = alpharavgbTS
        savevars['alphahavgbTS'] = alphahavgbTS
        savevars['DmrminbTS'] = DmrminbTS
        savevars['DmhminbTS'] = DmhminbTS
        savevars['alpharminbTS'] = alpharminbTS
        savevars['alphahminbTS'] = alphahminbTS
        
        savevars['qrmaxbtorTS'] = qrmaxbtorTS
        savevars['qravgbtorTS'] = qravgbtorTS
        savevars['qcmaxbtorTS'] = qcmaxbtorTS
        savevars['qcavgbtorTS'] = qcavgbtorTS
        savevars['qhmaxbtorTS'] = qhmaxbtorTS
        savevars['qhavgbtorTS'] = qhavgbtorTS
        savevars['evapqcmaxbtorTS'] = evapqcmaxbtorTS
        savevars['evapqcavgbtorTS'] = evapqcavgbtorTS
        savevars['evapqrmaxbtorTS'] = evapqrmaxbtorTS
        savevars['evapqravgbtorTS'] = evapqravgbtorTS
        savevars['meltqhmaxbtorTS'] = meltqhmaxbtorTS
        savevars['meltqhavgbtorTS'] = meltqhavgbtorTS
        savevars['DmrmaxbtorTS'] = DmrmaxbtorTS
        savevars['DmravgbtorTS'] = DmravgbtorTS
        savevars['DmrminbtorTS'] = DmrminbtorTS
        savevars['DmhmaxbtorTS'] = DmhmaxbtorTS
        savevars['DmhavgbtorTS'] = DmhavgbtorTS
        savevars['DmhminbtorTS'] = DmhminbtorTS
        savevars['alpharmaxbtorTS'] = alpharmaxbtorTS
        savevars['alpharavgbtorTS'] = alpharavgbtorTS
        savevars['alpharminbtorTS'] = alpharminbtorTS
        savevars['alphahmaxbtorTS'] = alphahmaxbtorTS
        savevars['alphahavgbtorTS'] = alphahavgbtorTS
        savevars['alphahminbtorTS'] = alphahminbtorTS

        N.savez(npzfilename,**savevars)        

if(loadTS and plot_TS):
    if(plot_var):
        npzfilename = npzsavedir+'/'+runname+'_varTS.npz'
        print "Loading time series for averaged trajectory variables"
        varavgfile = N.load(npzfilename)
        
        zmaxbTS = varavgfile['zmaxbTS']
        zmaxaTS = varavgfile['zmaxaTS']
        varminbTS = varavgfile['varminbTS']
        varminaTS = varavgfile['varminaTS']
        varmaxbTS = varavgfile['varmaxbTS']
        varmaxaTS = varavgfile['varmaxaTS']
        varavgbTS = varavgfile['varavgbTS']
        varavgaTS = varavgfile['varavgaTS']
        numtortrajcsTS = varavgfile['numtortrajcsTS']
        varminattorTS = varavgfile['varminattorTS']
        varavgattorTS = varavgfile['varavgattorTS']
        varmaxattorTS = varavgfile['varmaxattorTS']
        varminbtorTS = varavgfile['varminbtorTS']
        varmaxbtorTS = varavgfile['varmaxbtorTS']
        varavgbtorTS = varavgfile['varavgbtorTS']
        varminatorTS = varavgfile['varminatorTS']
        varmaxatorTS = varavgfile['varmaxatorTS']
        varavgatorTS = varavgfile['varavgatorTS']
        
    if(plot_wforce):
        npzfilename = npzsavedir+'/'+runname+'_wforceTS.npz'
        print "Loading time series for averaged trajectory vertical velocity forcing"
        wforceavgfile = N.load(npzfilename)
        
        wforceminbTS = wforceavgfile['wforceminbTS']
        wforceminaTS = wforceavgfile['wforceminaTS']
        wforcemaxbTS = wforceavgfile['wforcemaxbTS']
        wforcemaxaTS = wforceavgfile['wforcemaxaTS']
        wforceavgbTS = wforceavgfile['wforceavgbTS']
        wforceavgaTS = wforceavgfile['wforceavgaTS']
        wforceminbtorTS = wforceavgfile['wforceminbtorTS']
        wforceavgbtorTS = wforceavgfile['wforceavgbtorTS']
        wforcemaxbtorTS = wforceavgfile['wforcemaxbtorTS']
        wforceminatorTS = wforceavgfile['wforceminatorTS']
        wforceavgatorTS = wforceavgfile['wforceavgatorTS']
        wforcemaxatorTS = wforceavgfile['wforcemaxatorTS']
        
    if(plot_DSD):
        npzfilename = npzsavedir+'/'+runname+'_DSDTS.npz'
        print "Loading time series for averaged trajectory DSD variables"
        DSDavgfile = N.load(npzfilename)
        
        qrmaxbTS = DSDavgfile['qrmaxbTS']
        qcmaxbTS = DSDavgfile['qcmaxbTS']
        qhmaxbTS = DSDavgfile['qhmaxbTS']
        evapqcmaxbTS = DSDavgfile['evapqcmaxbTS']
        evapqrmaxbTS = DSDavgfile['evapqrmaxbTS']
        meltqhmaxbTS = DSDavgfile['meltqhmaxbTS']
        DmrmaxbTS = DSDavgfile['DmrmaxbTS']
        DmhmaxbTS = DSDavgfile['DmhmaxbTS']
        alpharmaxbTS = DSDavgfile['alpharmaxbTS']
        alphahmaxbTS = DSDavgfile['alphahmaxbTS']
        qravgbTS = DSDavgfile['qravgbTS']
        qcavgbTS = DSDavgfile['qcavgbTS']
        qhavgbTS = DSDavgfile['qhavgbTS']
        evapqcavgbTS = DSDavgfile['evapqcavgbTS']
        evapqravgbTS = DSDavgfile['evapqravgbTS']
        meltqhavgbTS = DSDavgfile['meltqhavgbTS']
        DmravgbTS = DSDavgfile['DmravgbTS']
        DmhavgbTS = DSDavgfile['DmhavgbTS']
        alpharavgbTS = DSDavgfile['alpharavgbTS']
        alphahavgbTS = DSDavgfile['alphahavgbTS']
        DmrminbTS = DSDavgfile['DmrminbTS']
        DmhminbTS = DSDavgfile['DmhminbTS']
        alpharminbTS = DSDavgfile['alpharminbTS']
        alphahminbTS = DSDavgfile['alphahminbTS']
        
        qrmaxbtorTS = DSDavgfile['qrmaxbtorTS']
        qravgbtorTS = DSDavgfile['qravgbtorTS']
        qcmaxbtorTS = DSDavgfile['qcmaxbtorTS']
        qcavgbtorTS = DSDavgfile['qcavgbtorTS']
        qhmaxbtorTS = DSDavgfile['qhmaxbtorTS']
        qhavgbtorTS = DSDavgfile['qhavgbtorTS']
        evapqcmaxbtorTS = DSDavgfile['evapqcmaxbtorTS']
        evapqcavgbtorTS = DSDavgfile['evapqcavgbtorTS']
        evapqrmaxbtorTS = DSDavgfile['evapqrmaxbtorTS']
        evapqravgbtorTS = DSDavgfile['evapqravgbtorTS']
        meltqhmaxbtorTS = DSDavgfile['meltqhmaxbtorTS']
        meltqhavgbtorTS = DSDavgfile['meltqhavgbtorTS']
        DmrmaxbtorTS = DSDavgfile['DmrmaxbtorTS']
        DmravgbtorTS = DSDavgfile['DmravgbtorTS']
        DmrminbtorTS = DSDavgfile['DmrminbtorTS']
        DmhmaxbtorTS = DSDavgfile['DmhmaxbtorTS']
        DmhavgbtorTS = DSDavgfile['DmhavgbtorTS']
        DmhminbtorTS = DSDavgfile['DmhminbtorTS']
        alpharmaxbtorTS = DSDavgfile['alpharmaxbtorTS']
        alpharavgbtorTS = DSDavgfile['alpharavgbtorTS']
        alpharminbtorTS = DSDavgfile['alpharminbtorTS']
        alphahmaxbtorTS = DSDavgfile['alphahmaxbtorTS']
        alphahavgbtorTS = DSDavgfile['alphahavgbtorTS']
        alphahminbtorTS = DSDavgfile['alphahminbtorTS']
        
if(plot_TS):

    if(plot_var):

        # Number of trajectories entering tornado over subsequent 5 min
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(tref_list,numtortrajcsTS)
        plt.title('Number of trajectories entering tornado over ensuing 5 min')
    
        # Fraction of trajectories in 12 km x 12 km square (2401 total) that enter tornado
        fractortrajcsTS = numtortrajcsTS/2401.0
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(tref_list,fractortrajcsTS)
        plt.title('Fraction of trajectories entering tornado over ensuing 5 min')

        # Maximum height (trajectory average) time series before reference time
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(tref_list,zmaxbTS)
        plt.title('Maximum height (past 15 min), trajectory average')

        for i,var in enumerate(varnames):
            print "Creating time series plots for variable "+var
    
            # Minimum value (trajectory average) time series before reference time
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(tref_list,varminbTS[:,i])
            plt.title('Minimum '+var+' (past 15 min), trajectory average')
            # Maximum value before reference time
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(tref_list,varmaxbTS[:,i])
            plt.title('Maximum '+var+' (past 15 min), trajectory average')
            # Average value before reference time
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(tref_list,varavgbTS[:,i])
            plt.title('Average '+var+' (past 15 min), trajectory average')
        
            if(usetorreftime):
                # Minimum value (over trajectories) of variable at point entering tornado
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.plot(tref_list,varminattorTS[:,i])
                plt.title('Minimum '+var+' at TLV')
                # Maximum value (over trajectories) of variable at point entering tornado
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.plot(tref_list,varmaxattorTS[:,i])
                plt.title('Maximum '+var+' at TLV')
                # Average value (over trajectories) of variable at point entering tornado
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.plot(tref_list,varavgattorTS[:,i])
                plt.title('Average '+var+' (past 15 min), trajectory average')
            
                # Minimum value (trajectory average) of variable before tornado
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.plot(tref_list,varminbtorTS[:,i])
                plt.title('Minimum '+var+' before TLV, trajectory average')
                # Maximum value (trajectory average) of variable before tornado
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.plot(tref_list,varmaxbtorTS[:,i])
                plt.title('Maximum '+var+' before TLV, trajectory average')
                # Average value (trajectory average) of variable before tornado
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.plot(tref_list,varavgbtorTS[:,i])
                plt.title('Average '+var+' before TLV, trajectory average')
            
                # Minimum value (trajectory average) of variable after tornado
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.plot(tref_list,varminatorTS[:,i])
                plt.title('Minimum '+var+' after TLV, trajectory average')
                # Maximum value (trajectory average) of variable after tornado
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.plot(tref_list,varmaxatorTS[:,i])
                plt.title('Maximum '+var+' after TLV, trajectory average')
                # Average value (trajectory average) of variable after tornado
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.plot(tref_list,varavgatorTS[:,i])
                plt.title('Average '+var+' after TLV, trajectory average')

    if(plot_wforce):
    
        for i,var in enumerate(wforce_varnames):
            print "Creating time series plots for variable "+var
        
            # Minimum value (trajectory average) time series before reference time
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(tref_list,wforceminbTS[:,i])
            plt.title('Minimum '+var+' (past 15 min), trajectory average')
            # Maximum value before reference time
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(tref_list,wforcemaxbTS[:,i])
            plt.title('Maximum '+var+' (past 15 min), trajectory average')
            # Average value before reference time
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(tref_list,wforceavgbTS[:,i])
            plt.title('Average '+var+' (past 15 min), trajectory average')

    if(plot_DSD):

        #Maximum qr (trajectory average) time series before reference time
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(tref_list,qrmaxbTS)
        plt.title('Maximum qr (past 15 min), trajectory average')
        #Maximum qc (trajectory average) time series before reference time
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(tref_list,qcmaxbTS)
        plt.title('Maximum qc (past 15 min), trajectory average')
        #Maximum qh (trajectory average) time series before reference time
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(tref_list,qhmaxbTS)
        plt.title('Maximum qh (past 15 min), trajectory average')
        #Maximum evapqc (trajectory average) time series before reference time
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(tref_list,evapqcmaxbTS)
        plt.title('Maximum evapqc (past 15 min), trajectory average')
        #Maximum evapqr (trajectory average) time series before reference time
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(tref_list,evapqrmaxbTS)
        plt.title('Maximum evapqr (past 15 min), trajectory average')
        #Maximum meltqh (trajectory average) time series before reference time
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(tref_list,meltqhmaxbTS)
        plt.title('Maximum meltqh (past 15 min), trajectory average')
        #Maximum Dmr (trajectory average) time series before reference time
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(tref_list,DmrmaxbTS)
        plt.title('Maximum Dmr (past 15 min), trajectory average')
        #Maximum Dmh (trajectory average) time series before reference time
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(tref_list,DmhmaxbTS)
        plt.title('Maximum Dmh (past 15 min), trajectory average')
        #Maximum alphar (trajectory average) time series before reference time
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(tref_list,alpharmaxbTS)
        plt.title('Maximum alphar (past 15 min), trajectory average')
        #Maximum alphah (trajectory average) time series before reference time
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(tref_list,alphahmaxbTS)
        plt.title('Maximum alphah (past 15 min), trajectory average')
        #Average qr (trajectory average) time series before reference time
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(tref_list,qravgbTS)
        plt.title('Average qr (past 15 min), trajectory average')
        #Average qc (trajectory average) time series before reference time
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(tref_list,qcavgbTS)
        plt.title('Average qc (past 15 min), trajectory average')
        #Average qh (trajectory average) time series before reference time
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(tref_list,qhavgbTS)
        plt.title('Average qh (past 15 min), trajectory average')
        #Average evapqr (trajectory average) time series before reference time
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(tref_list,evapqravgbTS)
        plt.title('Average evapqr (past 15 min), trajectory average')
        #Average evapqc (trajectory average) time series before reference time
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(tref_list,evapqcavgbTS)
        plt.title('Average evapqc (past 15 min), trajectory average')
        #Average meltqh (trajectory average) time series before reference time
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(tref_list,meltqhavgbTS)
        plt.title('Average meltqh (past 15 min), trajectory average')
        #Average Dmr (trajectory average) time series before reference time
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(tref_list,DmravgbTS)
        plt.title('Average Dmr (past 15 min), trajectory average')
        #Average Dmh (trajectory average) time series before reference time
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(tref_list,DmhavgbTS)
        plt.title('Average Dmh (past 15 min), trajectory average')
        #Average alphar (trajectory average) time series before reference time
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(tref_list,alpharavgbTS)
        plt.title('Average alphar (past 15 min), trajectory average')
        #Average alphah (trajectory average) time series before reference time
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(tref_list,alphahavgbTS)
        plt.title('Average alphah (past 15 min), trajectory average')
        #Minimum Dmr (trajectory average) time series before reference time
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(tref_list,DmrminbTS)
        plt.title('Minimum Dmr (past 15 min), trajectory average')
        #Minimum Dmh (trajectory average) time series before reference time
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(tref_list,DmhminbTS)
        plt.title('Minimum Dmh (past 15 min), trajectory average')
        #Minimum alphar (trajectory average) time series before reference time
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(tref_list,alpharminbTS)
        plt.title('Minimum alphar (past 15 min), trajectory average')
        #Minimum alphah (trajectory average) time series before reference time
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(tref_list,alphahminbTS)
        plt.title('Minimum alphah (past 15 min), trajectory average')
    
#plt.show()
    
    
    
    
    
    
    
