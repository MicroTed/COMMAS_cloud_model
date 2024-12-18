#-----------------------------------------------------------------------
#
# Makefile for Fortran codes
#
#-----------------------------------------------------------------------
# Need to define SUFFIXES for fortran stuff mainly

.SUFFIXES: .F .F90 .f90 .f .c

#-----------------------------------------------------------------------
#
# Trick:  If you want to rebuild a fortran module correctly after
#         changing the source, you need to put an explicit compile
#         action after the $*.o definition.  For some OS's, the 
#         make will not refer back to the automatic rules, it tries
#         to deal with the $*.mod file directly.  Since the module
#         files are derived from both $*.f90 and $*.F90  source files,
#         I did not see a way to make a rule to automatically invoke
#         the correct compile line.  Given the nature of Make, there
#         undoubtably is a way, but my universe has a finite life time.



#-----------------------------------------------------------------------
#
# Define microphysics packages - if you want special compile options for
# each, define them using a compiler flag like $(OPT_MICRO#) in the top
# level Makefile and pass them in.

MICRO1 = kessler
MICRO2 = lfoice
MICRO3 = warmlfo
MICRO5 = ice10.nci.03g
SRC = ../src
SRCTAK = ../takbin
# this allows use of a build directory at the same level as 'src'

#-----------------------------------------------------------------------
#
INCLUDE  = -I$(SRC)/include $(NETCDF_INC)

LIB_OBJS = \
           balloon.o             \
           binio.o               \
           box.o                 \
           buoy.o                \
           cline_module.o        \
           coriolis.o            \
           commasmpi_module.o    \
           cputime_module.o      \
           elec_module.o         \
           elec_driver.o         \
           file_module.o         \
           filter.o              \
           force_module.o        \
           grid_module.o         \
           gridio_module.o       \
           index_module.o        \
           initsubs.o            \
           mathutil.o            \
           micro_module.o        \
           namelist_module.o     \
           param_module.o        \
           print.o               \
           radardbz.o            \
           sndparams.o           \
           sam2v5d.titles.o      \
           smlstep.o             \
           stats.o               \
           string_module.o       \
           tick_module.o         \
           util.o                \
           raymond_lowpass.o     \
           traj_module.o         \
           vis5d_module.o        \
           v5dp.o                \

MIC_OBJS =  advectrk.o           \
           advectcrw.o           \
           boxadvect.o           \
           boxfallx.o            \
           crowfall.01ax.o       \
           mix.o                 \
           micro_driver.o        \
           module_mp_nssl_2mom.o \
           ddfallnx.o            \
           qvexcess.o            \
           rainconv.o            \
           rimconv.o             \
           ice10_setvt.o         \
           radflx.o              \
           sfcphy.o              \
           setvt_zvd.o           \
           solver.o              \
           tak_modules.o         \
           tdlbc.o               \
           vegebm.o              \
           vegpbl.o              \
           zieg3.driver.o        \
           ziegfall.o            \
           $(MICRO1).o           \
           $(MICRO2).o           \
           $(MICRO3).o           \
           icezvd_dr.o           \
           nssl_2mom_gs.o        \
           nucond.o              \
           micro_MY_driver.o     \
           my3mom_main_mod.o     \
           my3mom_fncs_mod.o     \
           my2mom_fncs_mod.o     \
           my3mom_sedi_mod.o     \
           module_mp_thompson.o  \
           module_mp_thompson411.o  \
           module_mp_thompson381wof.o  \
           module_mp_morr_two_moment.o \
           module_mp_radar.o

#
#           $(MICRO5).o           \
#            warmzieg01.o          \
#            warmzieg.driver.o     \
#            wenofall.o            \
           
#           icezvd_gs.o \
#           $(MICRO4).o           \
#           lfo3.driver.o         \

# HCM_OBJS =  sam.08.hcm.o         \
#             hcmdriver.o           \
#             solverhcm.o           


TAKCOMP =
#TAKCOMP = $(FPPDEF)-DTAKON

# TAK_OBJS =                       \
#            takdriver.o           \
#            takinterface.o        \
#            takinit.o             \
#            takcondsn.o           \
#            takcollsn.o           \
#            bincracr.o            \
#            takcrim.o             \
#            takcmelt.o            \
#            ionatttak.o           \
#            solvertak.o           \
#            ionsteptak.o          \
#            nssl_2mom_gs_tak.o    

#           icezvd_gs_tak.o



ELEC_OBJS = \
           boxmgsetup.o          \
           efield.o              \
           efieldstag.o          \
           ionstep.o             \
           ionatt.o              \
           ionattz.o             \
           iondrift.o            \
           mud02.o               \
           mud02d.o              \
           mudcom.o              \
           mud2sp.o              \
           mud3sp.o              \
           mszlight.o            \
           saund.02b.o           \
           strzsor.o             
#            mudlight.o            \
#           spark.ion.2.6.o       \
#           spark.ion.2.6.mpi.o   \

IOTEST_OBJS = \
           file_module.o         \
           gridio_module.o       \
           grid_module.o         \
           commasmpi_module.o    \
           cline_module.o        \
           param_module.o



ifeq ($(ELECON), 1)
MODEL_OBJS = $(LIB_OBJS) $(MIC_OBJS) $(ELEC_OBJS) 
else
MODEL_OBJS = $(LIB_OBJS) $(MIC_OBJS) 
endif

$(info    ELECON is $(ELECON))

PLOT_OBJS = basp1p.o           \
            basp2p.o           \
            basp3px.o           \
            plot_module.o      \
            global_module.o dualpara.o
            
#commas_init:  commas_init.o libcommas.a
#	$(LDR) $(FREE) $(SMP) -o commas_init commas_init.o  libcommas.a $(NETCDF_LIB)

commas: libcommas.a commas.o # $(MODEL_OBJS)
	$(LDR) $(SMP) $(IPO) $(FREE) -o commas commas.o libcommas.a  $(NETCDF_LIB)

# commashcm: libcommas.a commas.o $(HCM_OBJS)
# 	$(LDR) $(SMP) $(IPO) $(FREE) -o commashcm commas.o $(HCM_OBJS) libcommas.a  $(NETCDF_LIB)

# commastak: libcommas.a commas.o $(TAK_OBJS)
# 	$(LDR) $(SMP) $(IPO) $(FREE) -o commastak commas.o $(TAK_OBJS)  libcommas.a $(NETCDF_LIB)

commas_iotest: commas_iotest.o $(IOTEST_OBJS)
	$(LDR) $(SMP) $(IPO) $(FREE) -o commas_iotest commas_iotest.o $(IOTEST_OBJS)  $(NETCDF_LIB)

com2v5d: com2v5d.o libcommas.a
	$(LDR) $(SMP) $(FREE) -o com2v5d com2v5d.o libcommas.a  $(NETCDF_LIB)

comtraj: comtraj.o libcommas.a
	$(LDR) $(SMP) $(FREE) -o comtraj comtraj.o libcommas.a  $(NETCDF_LIB)

# comspectra: comspectra.o libcommas.a
# 	$(LDR) $(SMP) $(FREE) -o comspectra comspectra.o libcommas.a  $(NETCDF_LIB)

com2nc4: com2nc4.o libcommas.a
	$(LDR) $(SMP) $(FREE) -o com2nc4 com2nc4.o libcommas.a  $(NETCDF_LIB)

complot: complot.o libcommas.a $(PLOT_OBJS)
	$(NCARGF) $(SMP) $(FREE) -o complot complot.o $(PLOT_OBJS) libcommas.a $(NETCDF_LIB) 
#	$(NCARGF) $(SMP) $(FREE) -o complot complot.o libcommas.a $(NETCDF_LIB) $(PLOT_OBJS) $(NETCDF_LIB)

com2grads: com2grads.o libcommas.a
	$(LDR) $(SMP) $(FREE) -o com2grads com2grads.o libcommas.a $(NETCDF_LIB)

maxmin: maxmin.o libcommas.a
	$(LDR) $(SMP) $(FREE) -o maxmin maxmin.o libcommas.a  $(NETCDF_LIB)

MEMFLAG =
MEMUSED = 
#MEMUSED = ../src/getmemused.o
#MEMFLAG = -DGETMEM
combine: combine.o libcommas.a $(MEMUSED)
	$(LDR) $(SMP) $(FREE) -o combine combine.o libcommas.a $(MEMUSED) $(NETCDF_LIB)

com2trmm:  trmmdatdrive.o libcommas.a com2trmm.o
	$(LDR) $(SMP) $(FREE) -o com2trmm com2trmm.o trmmdatdrive.o libcommas.a $(NETCDF_LIB)

libcommas.a: $(MODEL_OBJS) 
	rm -f libcommas.a
	$(AR) libcommas.a $(MODEL_OBJS)
	ranlib libcommas.a


clean:
	rm -f *.o *.mod libcommas.a commas commas_init com2v5d *mod.f90 com2nc4 complot

#-----------------------------------------------------------------------
# Compilation rules  -- These don't seem to work with the SRC prefixes so need a rule for each

#.c.o:
#	$(CC) $(CFLAGS) $(INCLUDE) -c $(SRC)/$*.c

#.f90.o:
#	$(F90) $(FREE) $(SMP) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.f90

#.f.o:
#	$(F77) $(FIXED) $(SMP) $(OPT_F77) $(INCLUDE) -c $(SRC)/$*.f

#.F.o:
#	$(FPP) $(FIXED) $(SMP) $(DFLAGS) $(OPT_F77) $(INCLUDE) -c $(SRC)/$*.F

#.F90.o:
#	$(FPP) $(FREE) $(SMP) $(DFLAGS) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

#-----------------------------------------------------------------------
#
#
#-----------------------------------------------------------------------


commas.o:             $(SRC)/commas.F90 grid_module.o gridio_module.o namelist_module.o \
                      file_module.o cline_module.o micro_module.o param_module.o \
                      vis5d_module.o cputime_module.o index_module.o commasmpi_module.o \
                      tick_module.o
	$(LDR) $(FREE) $(SMP) $(DFLAGS) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90
#	$(FPP) $(FREE) $(SMP) $(DFLAGS) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

com2v5d.o:            $(SRC)/com2v5d.F90 grid_module.o gridio_module.o namelist_module.o \
                      file_module.o cline_module.o micro_module.o param_module.o \
                      vis5d_module.o cputime_module.o
	$(FPP) $(FREE) $(SMP) $(DFLAGS) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

comtraj.o:            $(SRC)/comtraj.F90 grid_module.o gridio_module.o namelist_module.o \
                      file_module.o cline_module.o micro_module.o param_module.o \
                      traj_module.o cputime_module.o
	$(FPP) $(FREE) $(SMP) $(DFLAGS) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

maxmin.o:             $(SRC)/maxmin.F90 grid_module.o gridio_module.o namelist_module.o \
                      file_module.o cline_module.o micro_module.o param_module.o \
                      vis5d_module.o cputime_module.o
	$(FPP) $(FREE) $(SMP) $(DFLAGS) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

comspectra.o:            $(SRC)/comspectra.F90 grid_module.o gridio_module.o namelist_module.o \
                      file_module.o cline_module.o micro_module.o param_module.o \
                      vis5d_module.o cputime_module.o
	$(FPP) $(FREE) $(SMP) $(DFLAGS) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

com2nc4.o:            $(SRC)/com2nc4.F90 grid_module.o gridio_module.o namelist_module.o \
                      file_module.o cline_module.o micro_module.o param_module.o \
                      vis5d_module.o cputime_module.o index_module.o commasmpi_module.o
	$(FPP) $(FREE) $(SMP) $(DFLAGS) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

combine.o:            $(SRC)/combine.F90 grid_module.o gridio_module.o namelist_module.o \
                      file_module.o cline_module.o micro_module.o param_module.o \
                      cputime_module.o index_module.o commasmpi_module.o
	$(FPP) $(FREE) $(SMP) $(DFLAGS) $(OPT_F90) $(INCLUDE) $(MEMFLAG) -c $(SRC)/$*.F90

complot.o:            $(SRC)/complot.F90 grid_module.o gridio_module.o namelist_module.o \
                      file_module.o cline_module.o micro_module.o param_module.o \
                      cputime_module.o plot_module.o
	$(FPP) $(FREE) $(SMP) $(DFLAGS) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

com2grads.o:     $(SRC)/com2grads.F90 grid_module.o gridio_module.o namelist_module.o \
                      file_module.o cline_module.o micro_module.o param_module.o
	$(FPP) $(FREE) $(SMP) $(DFLAGS) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90


#	$(FPP) $(FREE)  $(DFLAGS) $(INCLUDE) -O0 -g -c $*.F90

#commas_init.o:        $(SRC)/commas_init.F90 ens_module.o gridio_module.o namelist_module.o \
#                      file_module.o cline_module.o micro_module.o param_module.o commasmpi_module.o
#	$(FPP) $(FREE) $(SMP) $(DFLAGS) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

advectrk.o:           $(SRC)/advectrk.F90 grid_module.o param_module.o cputime_module.o commasmpi_module.o initsubs.o
	$(F90) $(FREE) $(SMP) $(DFLAGS) $(OPT_F90) $(OPT_ADV) $(INCLUDE) -c $(SRC)/$*.F90

advectcrw.o:          $(SRC)/advectcrw.F90 grid_module.o param_module.o cputime_module.o commasmpi_module.o
	$(FPP) $(FREE) $(SMP) $(DFLAGS) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

balloon.o:            $(SRC)/balloon.F90 grid_module.o index_module.o elec_module.o file_module.o 
	$(FPP) $(FREE) $(SMP) $(DFLAGS) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

binio.o:              $(SRC)/binio.c $(SRC)/include/binio.h
	$(CC) $(CFLAGS) $(INCLUDE) -c $(SRC)/$*.c

box.o:                $(SRC)/box.F90 grid_module.o commasmpi_module.o initsubs.o param_module.o
	$(FPP) $(FREE) $(SMP) $(DFLAGS) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

boxfallx.o:           $(SRC)/boxfallx.F90 commasmpi_module.o micro_module.o
	$(FPP) $(FREE) $(SMP) $(DFLAGS) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

boxmgsetup.o:         $(SRC)/boxmgsetup.F commasmpi_module.o elec_module.o cputime_module.o grid_module.o
	$(FPP) $(FIXED) $(SMP) $(DFLAGS) $(OPT_F77) $(INCLUDE) -c $(SRC)/$*.F

boxadvect.o:          $(SRC)/boxadvect.F commasmpi_module.o
	$(FPP) $(FIXED) $(SMP) $(DFLAGS) $(OPT_F77) $(INCLUDE) -c $(SRC)/$*.F

basp1p.o:             $(SRC)/basp1p.F plot_module.o grid_module.o param_module.o micro_module.o index_module.o
	$(FPP) $(FIXED) $(INCLUDE) $(OPT_F77) -c $(SRC)/$*.F

basp2p.o:             $(SRC)/basp2p.F plot_module.o
	$(FPP) $(FIXED) $(INCLUDE) $(OPT_F77) -c $(SRC)/$*.F

basp3px.o:            $(SRC)/basp3px.F plot_module.o grid_module.o micro_module.o index_module.o
	$(FPP) $(FIXED) $(INCLUDE) $(OPT_F77) -c $(SRC)/$*.F
#basp3px.o:            $(SRC)/basp3px.F90 plot_module.o grid_module.o micro_module.o index_module.o
#	$(FPP) $(FREE) $(INCLUDE) $(OPT_F90) -c $(SRC)/$*.F90

global_module.o:            $(SRC)/global_module.f90
	$(FPP) $(FREE) $(INCLUDE) $(OPT_F90) -c $(SRC)/$*.f90

dualpara.o:            $(SRC)/dualpara.F90 global_module.o
	$(FPP) $(FREE) $(INCLUDE) $(OPT_F90) -DCOMMAS -c $(SRC)/$*.F90

raymond_lowpass.o:             $(SRC)/raymond_lowpass.f
	$(FPP) $(FIXED) $(INCLUDE) $(OPT_F77) -c $(SRC)/$*.f

buoy.o:               $(SRC)/buoy.F90 grid_module.o param_module.o commasmpi_module.o
	$(FPP) $(FREE) $(DFLAGS) $(OPT_F90) $(OPT_BUOY) $(INCLUDE) -c $(SRC)/$*.F90

cline_module.o:       $(SRC)/cline_module.F90 
	$(FPP) $(FREE) $(DFLAGS) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

commasmpi_module.o:   $(SRC)/commasmpi_module.F90
	$(FPP) $(FREE) $(DFLAGS) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

coriolis.o:           $(SRC)/coriolis.F90  param_module.o commasmpi_module.o
	$(FPP) $(FREE) $(DFLAGS)  $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

cputime_module.o:     $(SRC)/cputime_module.F90 param_module.o
	$(FPP) $(FREE) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

crowfall.01ax.o:      $(SRC)/crowfall.01ax.F commasmpi_module.o
	$(FPP) $(FIXED) $(SMP) $(DFLAGS) $(OPT_F77) $(INCLUDE) -c $(SRC)/$*.F

efield.o:             $(SRC)/efield.F grid_module.o index_module.o commasmpi_module.o
	$(FPP) $(FIXED) $(SMP) $(DFLAGS) $(OPT_F77) $(INCLUDE) -c $(SRC)/$*.F

efieldstag.o:         $(SRC)/efieldstag.F grid_module.o index_module.o commasmpi_module.o
	$(FPP) $(FIXED) $(SMP) $(DFLAGS) $(OPT_F77) $(INCLUDE) -c $(SRC)/$*.F

elec_module.o:        $(SRC)/elec_module.F90
	$(FPP) $(FREE) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

elec_driver.o:        $(SRC)/elec_driver.F90 index_module.o file_module.o vis5d_module.o \
                      grid_module.o micro_module.o cputime_module.o elec_module.o \
                      commasmpi_module.o namelist_module.o traj_module.o
	$(FPP) $(FREE) $(OPT_F90) $(INCLUDE)  $(DFLAGS) -c $(SRC)/$*.F90

file_module.o:        $(SRC)/file_module.F90 grid_module.o commasmpi_module.o
	$(FPP) $(FREE) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

filter.o:             $(SRC)/filter.F90 grid_module.o param_module.o commasmpi_module.o
	$(FPP) $(FREE) $(SMP) $(DFLAGS) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

force_module.o:       $(SRC)/force_module.F90 param_module.o commasmpi_module.o tak_modules.o \
                      index_module.o micro_module.o initsubs.o
	$(FPP) $(FREE) $(DFLAGS) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

grid_module.o:        $(SRC)/grid_module.F90 $(SRC)/grid_stat.F90 $(SRC)/grid_variable.F90 $(SRC)/grid_attribute.F90 \
                      commasmpi_module.o
	$(FPP) $(FREE)  $(DFLAGS) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

gridio_module.o:      $(SRC)/gridio_module.F90 string_module.o  grid_module.o file_module.o commasmpi_module.o \
                      $(SRC)/grid_create.F90 $(SRC)/grid_netcdf.F90 $(SRC)/grid_pnetcdf.F90 $(SRC)/grid_hdf5.F90 \
                      $(SRC)/grid_binary.F90 index_module.o param_module.o commasmpi_module.o micro_module.o \
                      elec_module.o namelist_module.o force_module.o
	$(FPP) $(FREE)  $(DFLAGS) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

hcmdriver.o:    $(SRC)/micro_driver.F90 index_module.o file_module.o vis5d_module.o \
                      grid_module.o micro_module.o cputime_module.o elec_module.o commasmpi_module.o  elec_driver.o
	$(FPP) $(OPT_F90) $(FIXED) $(INCLUDE) -o $*.o $(DFLAGS) $(FPPDEF)-DHCMON  $(FPPDEF)-DCHGELEC=0 -c $(SRC)/micro_driver.F90

takdriver.o:    $(SRC)/micro_driver.F90 index_module.o file_module.o vis5d_module.o \
                      grid_module.o micro_module.o cputime_module.o elec_module.o commasmpi_module.o  elec_driver.o
	$(FPP) $(OPT_F90) $(FREE) $(INCLUDE) -o $*.o $(DFLAGS) $(FPPDEF)-DTAKON  $(FPPDEF)-DCHGELEC=0 -c $(SRC)/micro_driver.F90

takinterface.o: $(SRCTAK)/takinterface.F90 tak_modules.o index_module.o micro_module.o grid_module.o elec_module.o \
                      commasmpi_module.o cputime_module.o traj_module.o
	$(FPP) $(OPT_F90) $(FREE) $(INCLUDE) $(DFLAGS) -c $(SRCTAK)/$*.F90

takinit.o: $(SRCTAK)/takinit.F90 tak_modules.o index_module.o micro_module.o
	$(FPP) $(OPT_F90) $(FREE)  $(DFLAGS) $(INCLUDE) -c $(SRCTAK)/$*.F90

takcondsn.o: $(SRCTAK)/takcondsn.F90 tak_modules.o index_module.o micro_module.o traj_module.o
	$(FPP) $(OPT_F90) $(FREE)  $(DFLAGS) $(INCLUDE) -c $(SRCTAK)/$*.F90

takcollsn.o: $(SRCTAK)/takcollsn.F90 tak_modules.o index_module.o micro_module.o
	$(FPP) $(OPT_F90) $(FREE)  $(DFLAGS) $(INCLUDE) -c $(SRCTAK)/$*.F90

bincracr.o: $(SRCTAK)/bincracr.F90 tak_modules.o index_module.o micro_module.o
	$(FPP) $(OPT_F90) $(FREE)  $(DFLAGS) $(INCLUDE) -c $(SRCTAK)/$*.F90

takajsmnt.o: $(SRC)/takajsmnt.F90 tak_modules.o index_module.o 
	$(FPP) $(OPT_F90) $(FREE)  $(DFLAGS) $(INCLUDE) -c $(SRCTAK)/$*.F90

takcrim.o: $(SRCTAK)/takcrim.F90 tak_modules.o index_module.o micro_module.o elec_module.o grid_module.o traj_module.o
	$(FPP) $(OPT_F90) $(FREE)  $(DFLAGS) $(INCLUDE) -c $(SRCTAK)/$*.F90

takcmelt.o: $(SRCTAK)/takcmelt.F90 tak_modules.o index_module.o micro_module.o traj_module.o
	$(FPP) $(OPT_F90) $(FREE) $(DFLAGS) $(INCLUDE) -c $(SRCTAK)/$*.F90


tak_modules.o: $(SRC)/tak_modules.F90 $(SRC)/tak_modules_include.F90 index_module.o micro_module.o
	$(FPP) $(OPT_F90) $(FREE) $(DFLAGS) $(INCLUDE) -c $(SRC)/$*.F90

initsubs.o:           $(SRC)/initsubs.F90 micro_module.o file_module.o param_module.o grid_module.o \
                      elec_module.o commasmpi_module.o index_module.o
	$(FPP) $(FREE) $(SMP) $(DFLAGS) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

ice10.nci.03g.o:      $(SRC)/ice10.nci.03g.F micro_module.o $(SRC)/include/sam.def.h index_module.o \
                      grid_module.o elec_module.o
	$(F77) $(FIXED) $(MOPT5) $(DFLAGS) $(IPO) $(INCLUDE) -c $(SRC)/$*.F

icezvd_dr.o:          $(SRC)/icezvd_dr.F90 index_module.o micro_module.o grid_module.o elec_module.o \
                      commasmpi_module.o cputime_module.o traj_module.o param_module.o
	$(FPP) $(FREE) $(MOPT4) $(DFLAGS) -DCOMMAS  $(IPO) $(INCLUDE) -c $(SRC)/$*.F90


nssl_2mom_gs.o:       $(SRC)/nssl_2mom_gs.F90 index_module.o micro_module.o grid_module.o elec_module.o \
                      commasmpi_module.o cputime_module.o traj_module.o tak_modules.o namelist_module.o
	rm -f $*.tmp.F90
	cpp -nostdinc -C -P -traditional-cpp $(DFLAGS) -DCOMMAS $(SRC)/nssl_2mom_gs.F90 > $*.tmp.F90
	$(FPP) $(FREE) $(MOPT4) $(DFLAGS) $(FPPDEF)-DCOMMAS $(IPO) $(INCLUDE) -c $*.tmp.F90 -o $*.o
#	$(FPP) $(FIXED) $(MOPT4) $(DFLAGS) -DCOMMAS $(IPO) $(INCLUDE) -c $(SRC)/$*.F90

nssl_2mom_gs_tak.o:       $(SRC)/nssl_2mom_gs.F90 index_module.o micro_module.o grid_module.o elec_module.o \
                      commasmpi_module.o cputime_module.o traj_module.o tak_modules.o
	rm -f $*.tmp.F90
	cpp -nostdinc -C -P -traditional-cpp $(DFLAGS) $(FPPDEF)-DTAKON $(FPPDEF)-DCOMMAS  $(SRC)/nssl_2mom_gs.F90 > $*.tmp.F90
	$(FPP) $(FREE) $(MOPT4) $(DFLAGS) $(FPPDEF)-DTAKON $(FPPDEF)-DCOMMAS $(IPO) $(INCLUDE) -c $*.tmp.F90 -o $*.o

nucond.o:       $(SRC)/nucond.F90 index_module.o micro_module.o grid_module.o traj_module.o
	rm -f $*.tmp.F90
	cpp -nostdinc -C -P -traditional-cpp $(DFLAGS) -DCOMMAS $(SRC)/nucond.F90 > $*.tmp.F90
	$(FPP) $(FREE) $(MOPT4) $(DFLAGS) -DCOMMAS $(IPO) $(INCLUDE) -c $*.tmp.F90 -o $*.o

icezvd_gs_tak.o:          $(SRC)/icezvd_gs.F90 index_module.o micro_module.o grid_module.o elec_module.o \
                      commasmpi_module.o cputime_module.o traj_module.o tak_modules.o
	$(FPP) $(FREE) $(MOPT4) $(DFLAGS) $(IPO) $(INCLUDE) $(FPPDEF)-DTAKON -c -o $*.o $(SRC)/icezvd_gs.F90

sam.08.hcm.o:          $(SRCTAK)/sam.08.hcm.F90 index_module.o micro_module.o grid_module.o elec_module.o \
                      commasmpi_module.o cputime_module.o traj_module.o
	$(FPP) $(FIXED) $(MOPT4) $(DFLAGS) $(IPO) $(INCLUDE) -c $(SRCTAK)/$*.F90

index_module.o:       $(SRC)/index_module.F90
	$(FPP) $(FREE) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

ionstep.o:            $(SRC)/ionstep.F90 grid_module.o index_module.o elec_module.o commasmpi_module.o \
                      cputime_module.o micro_module.o
	$(FPP) $(FREE) $(SMP) $(DFLAGS) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

ionsteptak.o:         $(SRC)/ionstep.F90 grid_module.o index_module.o elec_module.o commasmpi_module.o cputime_module.o
	$(FPP) $(FREE) $(SMP) $(DFLAGS)  $(FPPDEF)-DTAKON $(OPT_F77) $(INCLUDE)  -o $*.o -c $(SRC)/ionstep.F90

ionatt.o:             $(SRC)/ionatt.F index_module.o micro_module.o commasmpi_module.o
	$(FPP) $(FIXED) $(SMP) $(DFLAGS) $(OPT_F77) $(INCLUDE) -c $(SRC)/$*.F

ionattz.o:            $(SRC)/ionattz.F90 elec_module.o index_module.o micro_module.o commasmpi_module.o
	$(FPP) $(FREE) $(SMP)  $(DFLAGS)  $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

ionatttak.o:          $(SRCTAK)/ionatttak.F90 index_module.o micro_module.o commasmpi_module.o tak_modules.o
	$(FPP) $(FREE) $(SMP)  $(DFLAGS)  $(OPT_F90) $(INCLUDE) -c $(SRCTAK)/$*.F90

iondrift.o:           $(SRC)/iondrift.F90 grid_module.o index_module.o commasmpi_module.o param_module.o
	$(FPP) $(FREE) $(SMP) $(DFLAGS) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

kessler.o:            $(SRC)/kessler.F90 param_module.o commasmpi_module.o traj_module.o
	$(FPP) $(FREE) $(SMP) $(DFLAGS) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

lfoice.o:             $(SRC)/lfoice.F param_module.o index_module.o commasmpi_module.o traj_module.o
	$(FPP) $(FIXED) $(SMP) $(DFLAGS) $(OPT_F77) $(INCLUDE) -c $(SRC)/$*.F

mathutil.o:           $(SRC)/mathutil.F90 micro_module.o param_module.o index_module.o
	$(FPP) $(FREE) $(SMP) $(OPT_F90) $(IPO) $(INCLUDE) -c $(SRC)/$*.F90

micro_module.o:       $(SRC)/micro_module.F90 param_module.o commasmpi_module.o file_module.o index_module.o elec_module.o
	$(FPP) $(FREE) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

mix.o:                $(SRC)/mix.F90 param_module.o grid_module.o commasmpi_module.o force_module.o micro_module.o
	$(FPP) $(FREE) $(SMP) $(DFLAGS) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

mud02.o:              $(SRC)/mud02.F elec_module.o grid_module.o index_module.o
	$(FPP) $(FIXED) $(SMP) $(DFLAGS) $(OPT_F77) $(INCLUDE) -c $(SRC)/$*.F

mud02d.o:             $(SRC)/mud02d.F elec_module.o grid_module.o index_module.o
	$(FPP) $(FIXED) $(SMP) $(DFLAGS) $(OPT_F77) $(INCLUDE) -c $(SRC)/$*.F

mudcom.o:             $(SRC)/mudcom.F
	$(FPP) $(FIXED) $(SMP) $(SOROPT) -c $(SRC)/$*.F

mud2sp.o:             $(SRC)/mud2sp.F
	$(FPP) $(FIXED) $(SMP) $(SOROPT) -c $(SRC)/$*.F

mud3sp.o:             $(SRC)/mud3sp.F
	$(FPP) $(FIXED) $(SMP) $(SOROPT) -c $(SRC)/$*.F

mudlight.o:           $(SRC)/mudlight.F90 grid_module.o index_module.o elec_module.o
	$(FPP) $(FREE) $(SMP) $(DFLAGS) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

mszlight.o:           $(SRC)/mszlight.F90 grid_module.o index_module.o elec_module.o
	$(FPP) $(FREE) $(SMP) $(DFLAGS) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

namelist_module.o:    $(SRC)/namelist_module.F90 string_module.o param_module.o micro_module.o \
                      elec_module.o force_module.o balloon.o initsubs.o traj_module.o \
                      commasmpi_module.o box.o vis5d_module.o
	$(FPP) $(FREE) $(DFLAGS) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

param_module.o:       $(SRC)/param_module.F90
	$(FPP) $(FREE) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

plot_module.o:     $(SRC)/plot_module.F90 global_module.o dualpara.o
	$(FPP) $(FREE) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

print.o:              $(SRC)/print.F90 grid_module.o file_module.o param_module.o micro_module.o \
                      commasmpi_module.o
	$(FPP) $(FREE) $(DFLAGS) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

radardbz.o:           $(SRC)/radardbz.F90 index_module.o commasmpi_module.o micro_module.o
	$(FPP) $(FREE) $(SMP) $(DFLAGS) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

rimdens.o:            $(SRC)/rimdens.F
	$(FPP) $(FIXED) $(SMP) $(DFLAGS) $(OPT_F77) $(INCLUDE) -c $(SRC)/$*.F

micro_driver.o: $(SRC)/micro_driver.F90 $(SRC)/include/sam.def.h \
                      vis5d_module.o file_module.o grid_module.o micro_module.o \
                      index_module.o cputime_module.o elec_module.o commasmpi_module.o elec_driver.o
	$(FPP) $(FREE) $(SMP) $(DFLAGS) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

micro_MY_driver.o: $(SRC)/micro_MY_driver.F90 file_module.o grid_module.o micro_module.o \
                         index_module.o cputime_module.o my3mom_fncs_mod.o my2mom_fncs_mod.o \
                         my3mom_sedi_mod.o my3mom_main_mod.o commasmpi_module.o
	$(FPP) $(FREE) $(SMP) $(DFLAGS) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

my3mom_fncs_mod.o: $(SRC)/my3mom_fncs_mod.F90
	$(FPP) $(FREE) $(SMP) $(DFLAGS) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

my2mom_fncs_mod.o: $(SRC)/my2mom_fncs_mod.F90
	$(FPP) $(FREE) $(SMP) $(DFLAGS) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

my3mom_sedi_mod.o: $(SRC)/my3mom_sedi_mod.F90 my2mom_fncs_mod.o
	$(FPP) $(FREE) $(SMP) $(DFLAGS) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

my3mom_main_mod.o: $(SRC)/my3mom_main_mod.F90 my3mom_sedi_mod.o
	$(FPP) $(FREE) $(SMP) $(DFLAGS) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

ddfallnx.o:   $(SRC)/ddfallnx.F $(SRC)/include/sam.def.h index_module.o
	$(FPP) $(FIXED) $(SMP) $(DFLAGS) $(OPT_F77) $(INCLUDE) -c $(SRC)/$*.F

ice10_setvt.o:          $(SRC)/ice10_setvt.F index_module.o commasmpi_module.o
	$(FPP) $(FIXED) $(SMP) $(DFLAGS) $(OPT_F77) $(IPO) $(INCLUDE) -c $(SRC)/$*.F

qvexcess.o:     $(SRC)/qvexcess.F90 micro_module.o param_module.o
	$(FPP) $(FREE) $(SMP) $(DFLAGS) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

rainconv.o:       $(SRC)/rainconv.F90
	$(FPP) $(FREE) $(SMP) $(OPT_F90) $(IPO) $(INCLUDE) -c $(SRC)/$*.F90

radflx.o:  $(SRC)/radflx.F90 param_module.o sfcphy.o namelist_module.o index_module.o
	$(FPP) $(FREE) $(DFLAGS) $(OPT_F90)  $(INCLUDE) -c $(SRC)/$*.F90

rimconv.o:        $(SRC)/rimconv.F
	$(FPP) $(FIXED) $(SMP) $(OPT_F77) $(IPO) $(INCLUDE) -c $(SRC)/$*.F

sam2v5d.titles.o:     $(SRC)/sam2v5d.titles.F
	$(FPP) $(FIXED) $(SMP) $(DFLAGS) $(OPT_F77) $(INCLUDE) -c $(SRC)/$*.F

saund.02b.o:          $(SRC)/saund.02b.F90
	$(FPP) $(FREE) $(SMP) $(DFLAGS) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

setvt_zvd.o:      $(SRC)/setvt_zvd.F90 index_module.o micro_module.o
	$(FPP) $(FREE) $(SMP) $(DFLAGS) $(OPT_F90) $(IPO) $(INCLUDE) -c $(SRC)/$*.F90

sfcphy.o:  $(SRC)/sfcphy.F90 param_module.o grid_module.o
	$(FPP) $(FREE) $(DFLAGS) $(OPT_F90)  $(INCLUDE) -c $(SRC)/$*.F90

sndparams.o:          $(SRC)/sndparams.F 
	$(FPP) $(FIXED) $(SMP) $(DFLAGS) $(OPT_F77) $(INCLUDE) -c $(SRC)/$*.F

smlstep.o:            $(SRC)/smlstep.F90 param_module.o commasmpi_module.o cputime_module.o
	$(FPP) $(FREE) $(SMP) $(DFLAGS) $(OPT_F90) $(OPT_ADV)  $(INCLUDE) -c $(SRC)/$*.F90

solver.o:             $(SRC)/solver.F90 param_module.o grid_module.o cputime_module.o micro_module.o \
                      index_module.o force_module.o traj_module.o commasmpi_module.o \
                      module_mp_nssl_2mom.o module_mp_thompson.o initsubs.o zieg3.driver.o \
                      module_mp_morr_two_moment.o module_mp_thompson411.o module_mp_thompson381wof.o
	$(FPP) $(FREE) $(DFLAGS) $(OPT_F90) $(OPT_SLV) $(GRDOPT)  $(INCLUDE) -c $(SRC)/$*.F90

solverhcm.o:          $(SRC)/solver.F90 param_module.o grid_module.o cputime_module.o micro_module.o \
                      index_module.o force_module.o traj_module.o initsubs.o commasmpi_module.o
	$(FPP) $(FREE) $(DFLAGS) -DHCMON $(OPT_F90) $(OPT_SLV) $(GRDOPT)  $(INCLUDE) -o $*.o -c $(SRC)/solver.F90

solvertak.o:          $(SRC)/solver.F90 param_module.o grid_module.o cputime_module.o micro_module.o \
                      initsubs.o index_module.o force_module.o traj_module.o commasmpi_module.o tak_modules.o
	$(FPP) $(FREE) $(DFLAGS) $(FPPDEF)-DTAKON $(OPT_F90) $(OPT_SLV) $(GRDOPT)  $(INCLUDE) -o $*.o -c $(SRC)/solver.F90

solverrk2.o:          $(SRC)/solverrk2.F90 param_module.o grid_module.o cputime_module.o micro_module.o \
                      index_module.o force_module.o
	$(FPP) $(FREE) $(SMP) $(DFLAGS) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

spark.ion.2.6.o:      $(SRC)/spark.ion.2.6.F elec_module.o
	$(FPP) $(FIXED) $(SOROPT) $(SMP) $(DFLAGS)  $(INCLUDE)  -c $(SRC)/$*.F

spark.ion.2.6.mpi.o:  $(SRC)/spark.ion.2.6.mpi.F90 elec_module.o commasmpi_module.o
	$(FPP) $(FREE) $(SOROPT) $(SMP) $(DFLAGS) -c $(INCLUDE)  $(SRC)/$*.F90

#spark.ion.2.6.mpi.o:  $(SRC)/spark.ion.2.6.mpi.F elec_module.o commasmpi_module.o
#	$(FPP) $(FIXED) $(SOROPT) $(SMP) $(DFLAGS) -c $(SRC)/$*.F

stats.o:              $(SRC)/stats.F90 commasmpi_module.o grid_module.o param_module.o
	$(FPP) $(FREE) $(SMP) $(DFLAGS) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

string_module.o:      $(SRC)/string_module.F90
	$(FPP) $(FREE) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

strzsor.o:            $(SRC)/strzsor.F commasmpi_module.o elec_module.o
	$(FPP) $(FIXED) $(SMP) $(DFLAGS) $(SOROPT) $(INCLUDE)  -c $(SRC)/$*.F

tdlbc.o:             $(SRC)/tdlbc.F90 param_module.o grid_module.o cputime_module.o micro_module.o \
                       file_module.o commasmpi_module.o initsubs.o
	$(FPP) $(FREE) $(DFLAGS) $(OPT_F90) $(OPT_SLV) $(GRDOPT)  $(INCLUDE) -c $(SRC)/$*.F90

tick_module.o:        $(SRC)/tick_module.F90
	$(FPP) $(FREE) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

traj_module.o:       $(SRC)/traj_module.F90 param_module.o commasmpi_module.o tak_modules.o grid_module.o \
                     micro_module.o initsubs.o
	$(FPP) $(FREE) $(DFLAGS) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

util.o:               $(SRC)/util.F90
	$(FPP) $(FREE) $(SMP) $(DFLAGS) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

v5dp.o:               $(SRC)/v5dp.c  $(SRC)/include/v5d.h $(SRC)/include/binio.h $(SRC)/include/vis5d.h
	$(CC) $(CFLAGS) $(INCLUDE) -c $(SRC)/$*.c

vegebm.o:  $(SRC)/vegebm.F90 sfcphy.o param_module.o
	$(FPP) $(FREE) $(DFLAGS) $(OPT_F90)  $(INCLUDE) -c $(SRC)/$*.F90

vegpbl.o:  $(SRC)/vegpbl.F90 sfcphy.o param_module.o namelist_module.o
	$(FPP) $(FREE) $(DFLAGS) $(OPT_F90)  $(INCLUDE) -c $(SRC)/$*.F90

vis5d_module.o:       $(SRC)/vis5d_module.F90 $(SRC)/include/v5df90.h grid_module.o param_module.o \
                      micro_module.o index_module.o commasmpi_module.o tak_modules.o
	$(F90) $(FREE) $(DFLAGS) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

warmlfo.o:            $(SRC)/warmlfo.F param_module.o micro_module.o commasmpi_module.o traj_module.o
	$(FPP) $(FIXED) $(SMP) $(DFLAGS) $(OPT_F77) $(INCLUDE) -c $(SRC)/$*.F

warmzieg01.o:         $(SRC)/warmzieg01.F90 index_module.o micro_module.o commasmpi_module.o cputime_module.o \
                      traj_module.o
	$(FPP) $(FREE) $(MOPT4) $(DFLAGS) $(IPO) $(INCLUDE) -c $(SRC)/$*.F90

warmzieg.driver.o:    $(SRC)/micro_driver.F90 index_module.o file_module.o vis5d_module.o \
                      grid_module.o micro_module.o cputime_module.o elec_module.o commasmpi_module.o  elec_driver.o
	$(FPP) $(OPT_F90) $(FREE) $(INCLUDE) -o $*.o $(DFLAGS) $(FPPDEF)-DWARM -c $(SRC)/micro_driver.F90

module_mp_thompson.o:  $(SRC)/module_mp_thompson.F module_mp_radar.o
	$(FPP) $(FREE) $(OPT_F90) $(DFLAGS) -DCOMMAS $(IPO) $(INCLUDE) -c $(SRC)/$*.F
#	$(FPP) $(FREE) -O0 -g -C $(DFLAGS) -DCM1 -DCOMMAS $(IPO) $(INCLUDE) -c $(SRC)/$*.F

module_mp_thompson381wof.o:  $(SRC)/module_mp_thompson381wof.F module_mp_radar.o
	$(FPP) $(FREE) $(OPT_F90) $(DFLAGS) -DCM1 -DCOMMAS $(IPO) $(INCLUDE) -c $(SRC)/$*.F
#	$(FPP) $(FREE) -O0 -g -C $(DFLAGS) -DCM1 -DCOMMAS $(IPO) $(INCLUDE) -c $(SRC)/$*.F

module_mp_thompson411.o:  $(SRC)/module_mp_thompson411.F module_mp_radar.o
	$(FPP) $(FREE) $(OPT_F90) $(DFLAGS) -DCM1 -DCOMMAS $(IPO) $(INCLUDE) -c $(SRC)/$*.F
#	$(FPP) $(FREE) -O0 -g -C $(DFLAGS) -DCM1 -DCOMMAS $(IPO) $(INCLUDE) -c $(SRC)/$*.F

module_mp_morr_two_moment.o:  $(SRC)/module_mp_morr_two_moment.F module_mp_radar.o micro_module.o
	$(FPP) $(FREE) $(OPT_F90) $(DFLAGS) -DCM1 -DCOMMAS $(IPO) $(INCLUDE) -c $(SRC)/$*.F
#	$(FPP) $(FREE) -O0 -g -C $(DFLAGS) -DCM1 -DCOMMAS $(IPO) $(INCLUDE) -c $(SRC)/$*.F

module_mp_radar.o:  $(SRC)/module_mp_radar.F
	$(FPP) $(FREE) $(MOPT4) $(DFLAGS) -DCM1 -DCOMMAS $(IPO) $(INCLUDE) -c $(SRC)/$*.F

module_mp_nssl_2mom.o:  $(SRC)/module_mp_nssl_2mom.F namelist_module.o
	$(FPP) $(FREE) $(MOPT4) $(DFLAGS) -DCM1 -DCOMMAS $(IPO) $(INCLUDE) -c $(SRC)/$*.F

# module_mp_nssl_2mom.o:  $(SRC)/module_mp_nssl_2mom.master.F $(SRC)/saund.02b.F90 $(SRC)/nssl_2mom_gs.F90 $(SRC)/nucond.F90 namelist_module.o
# 	rm -f $*.F
# 	cpp -nostdinc -C -P -traditional-cpp $(DFLAGS) -DCM1 -DCOMMAS $(SRC)/module_mp_nssl_2mom.master.F > $*.F
# 	$(FPP) $(FREE) $(MOPT4) $(DFLAGS) -DCM1 -DCOMMAS $(IPO) $(INCLUDE) -c $*.F
# #	$(FPP) $(FREE) $(MOPT4) $(DFLAGS) -DCM1 -DCOMMAS $(IPO) $(INCLUDE) -c $(SRC)/$*.F

wenofall.o:           $(SRC)/wenofall.F
	$(FPP) $(FIXED) $(SMP) $(DFLAGS) $(OPT_F77) $(INCLUDE) -c $(SRC)/$*.F

zieg3.driver.o:       $(SRC)/micro_driver.F90 index_module.o file_module.o vis5d_module.o \
                      grid_module.o micro_module.o cputime_module.o elec_module.o \
                      commasmpi_module.o elec_driver.o
	$(FPP) $(OPT_F90) $(FREE) $(INCLUDE) -o $*.o $(DFLAGS) $(FPPDEF)-DICE3 -c $(SRC)/micro_driver.F90

lfo3.driver.o:       $(SRC)/micro_driver.F90 index_module.o file_module.o vis5d_module.o \
                      grid_module.o micro_module.o cputime_module.o elec_module.o \
                      commasmpi_module.o elec_driver.o
	$(FPP) $(OPT_F90) $(FIXED) $(INCLUDE) -o $*.o $(DFLAGS) $(FPPDEF)-DICE3 $(FPPDEF)-DLFO3 -c $(SRC)/micro_driver.F90

ziegfall.o:           $(SRC)/ziegfall.F90 index_module.o commasmpi_module.o cputime_module.o micro_module.o
	$(FPP) $(FREE) $(SMP) $(DFLAGS) $(OPT_F90) $(INCLUDE) -c $(SRC)/$*.F90

