
ifndef $(MACH) 
	#MACH=MYWORKSTATION
	#MACH=MAC
	MACH=mpi
endif

include make.$(MACH)

######################## DIRECTORIES ###########################

SRCDIR       = ../source/
CVODE_SRCDIR = ../source/cvode/
PROG         = ../run/odt.x

######################## COMPILATION ###########################

includes = -I$(CANT_INC)
libs     = -L$(CANT_LIB)

ifeq ($(CANTERA), 1.7)
    D1 = -DDOCANTERA
    ext_libs =  -luser -loneD -lzeroD -ltransport -lcantera -lcvode \
                -lctmath -ltpx -lctf2c -lconverters -lctcxx 
endif
ifeq ($(CANTERA), 1.8)
    D1 = -DDOCANTERA -DCANTERA18
    ext_libs = -luser -loneD -lzeroD -lequil -lkinetics -ltransport -lthermo \
               -lctnumerics -lctmath -ltpx -lctspectra -lconverters -lctbase -lctcxx \
               -lcvode -lctf2c -lctf2c -lm 
endif
ifneq ($(MACH), MAC)
    ext_libs += -lctlapack -lctblas
endif

ifeq ($(CANTERA), 2.1)
    D1 = -DDOCANTERA -DCANTERA21
    ext_libs = 
    ext_libs = -lcantera
endif


ifneq ($(CHEMISTRY), USUAL)
    D2 = -DPROBLEMSPECIFICRR
endif

ifeq ($(MTWIST), YES)
    D3 = -DDOMTRAND
endif
ifeq ($(PARALLEL), YES)
    D4 = -DDOMPI
endif
ifeq ($(OLDMESHER), YES)
    D5 = -DOLDMESHER
endif

DIRECTIVES = $(D1) $(D2) $(D3) $(D4) $(D5)

COMPCC   = $(CC_COMPILER_NAME) $(includes) $(options) $(DIRECTIVES) -c 
COMPC    = $(C_COMPILER_NAME)  $(includes) $(options) $(DIRECTIVES) -c 

COMPDEP  = $(CC_COMPILER_NAME) $(includes) $(options) $(DIRECTIVES)  -MM -MT
COMPDEPC = $(C_COMPILER_NAME)  $(includes) $(options) $(DIRECTIVES)  -MM -MT

LINK     = $(CC_COMPILER_NAME) $(options) $(libs)

######################## SOURCE FILES ##########################

SRC  =  $(SRCDIR)main.cc                        \
        $(SRCDIR)odtline.cc                     \
        $(SRCDIR)anyline.cc                     \
        $(SRCDIR)odtParam.cc                    \
        $(SRCDIR)inputFile.cc                   \
        $(SRCDIR)odtSolver.cc                   \
        $(SRCDIR)eddy.cc                        \
        $(SRCDIR)dumper.cc                      \
        $(SRCDIR)diffuser.cc                    \
        $(SRCDIR)batchRXR.cc                    \
        $(SRCDIR)streams.cc                     \
        $(SRCDIR)processor.cc                   \
        $(SRCDIR)cantera_shell_functions.cc     \
        $(SRCDIR)radiation.cc                   \
        $(SRCDIR)particles.cc                   \
        $(SRCDIR)MOM_soot.cc                    \
        $(SRCDIR)MOM_default.cc                 \
        $(SRCDIR)ETA_default.cc                 \
        $(SRCDIR)ETA_tableLookup.cc             \
        $(SRCDIR)table.cc                       \
        $(SRCDIR)TableReader.cc                 \
        $(SRCDIR)stats.cc                       \
        $(SRCDIR)mom_Aq.cc                      \
        $(SRCDIR)ETA_Aq.cc                      \
        $(SRCDIR)pdgen.cc                      \
	$(SRCDIR)subdomain.cc 			\
	$(SRCDIR)linegroup.cc			\
	$(SRCDIR)tripletmap.cc			\
	$(SRCDIR)eddy_functions.cc		\
	$(SRCDIR)walltime.cc			\
     	$(SRCDIR)stopevent.cc 			\
      	$(SRCDIR)sdgroup.cc 			\
      	$(SRCDIR)point.cc 			\

CVSRC = $(CVODE_SRCDIR)sundials_smalldense.c    \
        $(CVODE_SRCDIR)sundials_nvector.c       \
        $(CVODE_SRCDIR)sundials_math.c          \
        $(CVODE_SRCDIR)sundials_dense.c         \
        $(CVODE_SRCDIR)nvector_serial.c         \
        $(CVODE_SRCDIR)cvode_io.c               \
        $(CVODE_SRCDIR)cvode_dense.c            \
        $(CVODE_SRCDIR)cvode.c                  

#-------Add in the desired mesher

ifeq ($(OLDMESHER), YES)
    SRC += $(SRCDIR)adaptMesh_old.cc
else
    SRC += $(SRCDIR)adaptMesh.cc
endif

#-------Add your rxn rate file here: make another ifeq block

ifeq ($(CHEMISTRY), ONESTEP)
    SRC += $(SRCDIR)oneStepRR.cc
endif
ifeq ($(CHEMISTRY), C2H4RED)
    SRC += $(SRCDIR)c2h4RedRR.cc
endif


######################## OBJECT FILES ##########################

OBJS   = $(SRC:.cc=.o)
CVOBJS = $(CVSRC:.c=.o)

######################## BUILD TARGETS ##########################

default : $(PROG)

all : depend default doxygen

$(PROG) : $(OBJS) $(CVOBJS)	
	$(LINK) $(OBJS) $(CVOBJS) -o $@ $(ext_libs)

%.o: %.cc 
	$(COMPCC) $< -o $@
%.o: %.c
	$(COMPC)  $< -o $@

clean :
	rm -f $(OBJS) $(CVOBJS) ./*.mod ./.depend
	rm -f ../run/odt.x
	rm -f ../run/*.log ../run/*.xml
	rm -rf ../doc/doxygen/html/*
	rm -f ../doc/doxygen/latex/*
	rm -f ../doc/doxygen/ODTmanual.pdf
	rm -f ../doc/doxygen/html/pdf/*

doxygen :
	./doDoxy.sh
	#./makePdfManual.sh
	
#doxygen Doxyfile > /dev/null 2>&1
    

####################### OBJECT DEPENDENCIES #####################

DEPENDS   = $(SRC:.cc=.d)
CVDEPENDS = $(CVSRC:.c=.cvd)

depend :  $(DEPENDS) $(CVDEPENDS)
	cat   $(DEPENDS) $(CVDEPENDS) > .depend
	rm -f $(DEPENDS) $(CVDEPENDS)

%.d :
	$(COMPDEP) $*.o $*.cc > $*.d
%.cvd :
	$(COMPDEPC) $*.o $*.c > $*.cvd

#ifeq ($(wildcard .depend), .depend)
-include .depend
#endif
