# Package information
LIBRARY = HOTVR
OBJDIR  = obj
DEPDIR  = $(OBJDIR)/dep
SRCDIR  = src
INCDIR  = include

INCLUDES += -I$(SFRAME_DIR)/NtupleWriter/include
INCLUDES += -I$(SFRAME_DIR)/NtupleWriter
INCLUDES += -I$(SFRAME_DIR)/SFrameTools/include
INCLUDES += -I$(SFRAME_DIR)/SFrameTools
INCLUDES += -I$(SFRAME_DIR)/SFrameTools/JetMETObjects/interface
INCLUDES += -I$(SFRAME_DIR)/SFrameAnalysis/include
INCLUDES += -I$(SFRAME_DIR)/SFrameAnalysis/
INCLUDES += -I$(SFRAME_DIR)/core
INCLUDES += -I$(SFRAME_DIR)/core/include

# configure FastJet
INCLUDES += -I$(FASTJETDIR)/include
INCLUDES += -I$(FASTJETDIR)/../include
INCLUDES += -I$(FASTJETDIR)/include/contribs/RecursiveTools


USERLDFLAGS += $(shell root-config --libs) -lMinuit
USERLDFLAGS +=  /nfs/dust/cms/user/tlapsien/fastjet-install/lib/libRecursiveTools.a
USERLDFLAGS +=  /nfs/dust/cms/user/tlapsien/fastjet-install/lib/libJetsWithoutJets.a
USERLDFLAGS +=  /nfs/dust/cms/user/tlapsien/fastjet-install/lib/libEnergyCorrelator.a
USERLDFLAGS +=  /nfs/dust/cms/user/tlapsien/fastjet-install/lib/libVariableR.a	
USERLDFLAGS +=  /nfs/dust/cms/user/tlapsien/fastjet-install/lib/libClusteringVetoPlugin.a
USERLDFLAGS +=  /nfs/dust/cms/user/tlapsien/fastjet-install/lib/libHHTopTagger.a
USERLDFLAGS +=  /nfs/dust/cms/user/tlapsien/fastjet-install/lib/libHOTVR.a
#USERLDFLAGS += /nfs/dust/cms/user/trappeu/fastjet-install/lib/libHHTopTagger.a
USERLDFLAGS +=  /nfs/dust/cms/user/tlapsien/fastjet-install/lib/libfastjetplugins.a

#INCLUDES += -I$(LHAPDFDIR)/include
INCLUDES += -I/afs/cern.ch/sw/lcg/external/MCGenerators/lhapdf/5.8.8/x86_64-slc5-gcc46-opt/include


# Include the generic compilation rules
include $(SFRAME_DIR)/Makefile.common
