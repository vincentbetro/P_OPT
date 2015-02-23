
include ../Makefile.options

LIBRARY_PATH += -L../SMOOTH

LIBS += -lgeometry -lsmooth -lUtility -loctree -lPolymesh -lsgio
#LIBS += -lgeometry -lsmooth -lUtility -loctree -lpcgns -lhdf5 -lmetis -lcgns -lPolymesh -lmesh_io -lpcgns -lcgns -lsgio
#LIBS += -lgeometry -lsmooth -lUtility -loctree -lcgns -lmetis /usr/local/hdf5-1.6.4-icc/lib/libhdf5.a

LDFLAGS = $(LIBRARY_PATH) $(LIBS)

OPTC_SRCS = P_OPT.cpp \

OPTSRCS = $(OPTC_SRCS)
OPTOBJECTS = $(OPTC_SRCS:.cpp=.o)

all: P_OPT

P_OPT: smooth $(OPTOBJECTS)
	$(LINK) -o $@ $(OPTOBJECTS) $(LDFLAGS) ../../CONV/create_hybrid_maps.o
	cp $@ ../$@.$(EXE_SUFFIX)
	cp $@ ../../$@.$(EXE_SUFFIX)

smooth:
	cd ../SMOOTH; $(MAKE)

clean:
	/bin/rm -f *.o
	/bin/rm -f P_OPT

include ../suffix.rules
