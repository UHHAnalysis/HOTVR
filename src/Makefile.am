
SUBDIRS = .

ACLOCAL_AMFLAGS=-I m4
AM_CPPFLAGS = -O4 -Wno-long-long @fjcxxflags@
AM_LDFLAGS = -O4
AM_LDADD = @fjlibs@

lib_LTLIBRARIES = libDeconstruction.la
libDeconstruction_la_CPPFLAGS = $(AM_CPPFLAGS)
libDeconstruction_la_SOURCES = AnalysisParameters.cxx ISRModel.cxx Model.cxx Storage.cxx TopGluonModel.cxx \
                               WModel.cxx BackgroundModel.cxx Exception.cxx JetInfo.cxx Parameters.cxx StoredCalculations.cxx \
                               BModel.cxx HBBModel.cxx \
                               StoredJet.cxx TopModel.cxx Deconstruct.cxx Helper.cxx \
                               Message.cxx


libDeconstruction_la_LIBADD = $(AM_LDADD)
libDeconstruction_la_LDFLAGS = $(AM_LDFLAGS) -module -avoid-version -shared


