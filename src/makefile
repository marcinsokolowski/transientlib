include flags.mak

LIBNAME=transient.a

HEADERS=writeimage.h  \

OBJECTS=cfgfile.o mystring.o cexcp.o myenv.o myfile.o basestring.o \
        myparser.o mytrace.o mydate.o random.o mypoints.o mystrtable.o\
 		  mykeytab.o mysafekeytab.o basestructs.o cmncfg.o mybits.o \
        paramtab.o tab2Ddesc.o myutil.o \
		  myfract.o basedefines.o mycmnglobals.o \
		  myinifile.o baseanal.o myhisto.o templ.o \
		  mylock.o myfastlongtab.o \
		  mydate2.o mymatrix.o \
		  cmn_tmpl.o laplace_info.o ccdrowcollist.o myframeslist.o \
		  myhisto2D.o mynet.o ccd_matrix.o Astroangle.o \
		  fits_file_base.o out_tmpl.o fitsout_trace.o \
         fitslib_globals.o ccd_globals.o ccd_trace.o cfg.o \
			calcrot.o ccd_hardware_defines.o ccd_fits_header_defs.o myfits.o\
			mathfunc.o myprogress.o ccd_util.o ccd.o ccd_analyse.o ccd_pipeline.o \
			ccd_report.o mypixellist.o ccd_defines.o ccd_single.o ccd_log.o ccd_datares.o \
			ccd_procstate.o AstroCCD.o myvalcounter.o ccd_eventlog.o ccd_starspy.o \
			mygraphcut.o ccd_state.o \
			ccd_common_struct.o ccd_fastcalc.o sat_info.o gendistr.o \
			ccd_image_creator.o ccd_tmpl.o ccd_samples.o ccd_asastransform.o \
			cfgtab.o ccd_ccdsystem.o ccd_pipeline2.o ccd_runstat.o \
			ccd_common_defines.o satcone.o AstroBody.o sat_interface.o \
			ccd_dataconst.o ccd_photometry.o mcinput.o ccd_corr_file.o \
			ccd_starcat.o ccd_starcat_cache.o

C_OBJECTS = astutil.o ricecomp.o asas_fitsio.o predict.o asas_astrometry.o \
				cat_new.o cat_old.o starcat.o
			
INCLUDES= -I/opt/pi/ext/include/ -I./ `root-config --cflags`


default : $(OBJECTS) $(C_OBJECTS)
	time -p gcc -v -shared -Wl,-soname,lib$(LIBNAME) -o sh$(LIBNAME) $(OBJECTS) $(C_OBJECTS)

install : all 
	mkdir -p /opt/pi/dev/pisys/daq/ndir/bin/
	mkdir -p /opt/pi/dev/pisys/daq/ndir/slib/
	mkdir -p $(NDIR)/slib/
	mkdir -p $(NDIR)/bin/
	cp *.a $(NDIR)/slib/
	cp find_flashes/find_flashes $(NDIR)/bin/
	cp sattest/sattest $(NDIR)/bin/
	chmod +x *!
	cp *! $(NDIR)/bin/

all   : lib
	cd find_flashes && make all
	cd sattest && make all

lib   : sh$(LIBNAME)
	cp sh$(LIBNAME) lib$(LIBNAME)
#	cp sh$(LIBNAME) lib$(LIBNAME)

sh$(LIBNAME) : $(IDL_GEN_FILES) $(OBJECTS) $(C_OBJECTS) $(CORBA_OBJECTS)
	ar sru $(LIBNAME) $^
	time -p gcc -v -shared -Wl,-soname,lib$(LIBNAME) -o sh$(LIBNAME) $(OBJECTS) $(CORBA_OBJECTS) $(C_OBJECTS) # `root-config --libs`
#	ld -o sh$(LIBNAME) $^ -shared
#	ranlib $(LIBNAME)

$(OBJECTS) : %.o :  %.cpp
	time -p $(CPP_COMP) $(INCLUDES) $(CMN_INCLUDES) $(PREPROC) $(CCFLAGS) -c $<

$(C_OBJECTS) : %.o :  %.c
	time -p $(C_COMP) $(CCFLAGS) -c $<
#	time -p $(C_COMP) $(INCLUDES) $(CMN_INCLUDES) $(PREPROC) $(CCFLAGS) -c $<



.PHONY : clean
clean :
	rm -f $(LIBNAME) $(OBJECTS) $(C_OBJECTS) sh$(LIBNAME) $(CORBA_OBJECTS) $(IDL_GEN_FILES) *.a
	cd find_flashes && make clean

