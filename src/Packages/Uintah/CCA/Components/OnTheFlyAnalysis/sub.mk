# Makefile fragment for this subdirectory

include $(SCIRUN_SCRIPTS)/smallso_prologue.mk

SRCDIR   := Packages/Uintah/CCA/Components/OnTheFlyAnalysis

SRCS     += \
	$(SRCDIR)/AnalysisModuleFactory.cc \
	$(SRCDIR)/AnalysisModule.cc	\
	$(SRCDIR)/lineExtract.cc	\
	$(SRCDIR)/flatPlate_heatFlux.cc	\
	$(SRCDIR)/pointExtract.cc	\
	$(SRCDIR)/containerExtract.cc	\
	$(SRCDIR)/vorticity.cc

PSELIBS := \
	Core/Geometry                    \
	Core/Exceptions                  \
	Core/Util                        \
	Core/OS                          \
	Packages/Uintah/CCA/Ports        \
	Packages/Uintah/Core/Disclosure  \
	Packages/Uintah/Core/Exceptions  \
	Packages/Uintah/Core/Grid        \
	Packages/Uintah/Core/Parallel    \
	Packages/Uintah/Core/Labels      \
	Packages/Uintah/Core/Parallel    \
	Packages/Uintah/Core/Util        \
	Packages/Uintah/Core/GeometryPiece \
	Packages/Uintah/Core/ProblemSpec 

LIBS := $(XML_LIBRARY) 

include $(SCIRUN_SCRIPTS)/smallso_epilogue.mk

