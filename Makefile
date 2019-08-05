CXX=g++
CXXFLAGS = -W -Wall -Wno-write-strings -Wno-extra -Wno-unused-parameter -Wno-unknown-pragmas -DAMS_ACQT_INTERFACE -D_PGTRACK_ -Wno-unused-variable -Wno-unused-function -g

# AMS Global Environment
CVMFS_AMS_OFFLINE = /cvmfs/ams.cern.ch/Offline

############################# CERN libraries ##################################
CERNLIBS = -lgeant321 -lpacklib -lmathlib -lkernlib
CERN_LEVEL = 2005.gcc64

ifndef CERNDIR
CERNDIR = $(CVMFS_AMS_OFFLINE)/CERN/$(CERN_LEVEL)
endif

CERNSRCDIR = $(CERNDIR)

ifndef AMSLIB
AMSLIB = /$(CVMFS_AMS_OFFLINE)/lib/linux/gcc64
endif

ifndef NAGDIR
NAGDIR = $(CVMFS_AMS_OFFLINE)/CERN/NagLib
endif
######################### End of CERN library settings ########################

##################### AMS Offline Software Related Includes ###################
INCLUDES = -I${ROOTSYS}/include -I${AMSWD}/include -I./include
NTUPLE_PG = $(AMSWD)/lib/linuxx8664gcc5.34/ntuple_slc6_PG.so
################## End of AMS Offline Software related includes ###############

############################# ROOT Related Settings ###########################
ROOTLIBS = $(shell root-config --libs) -lRIO -lNet -lNetx -lMinuit -lTMVA -lMLP -lXMLIO -lTreePlayer
########################## End of ROOT related settings #######################


TARGET = anomalous

all : $(TARGET)

$(TARGET) : anomalous.o
	$(CXX) $(CXXFLAGS) $(ROOTLIBS) $(INCLUDES) -Wno-extra -o $@ $^ $(NTUPLE_PG)

anomalous.o : anomalous.cxx
	$(CXX) $(CXXFLAGS) $(INCLUDES) -Wno-extra -c $^ -o $@

clean :
	rm -v anomalous.o
	rm -v anomalous
# DO NOT DELETE
#anomalous.o: anomalous.cxx
