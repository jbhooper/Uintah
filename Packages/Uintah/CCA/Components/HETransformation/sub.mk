# Makefile fragment for this subdirectory

include $(SRCTOP)/scripts/smallso_prologue.mk

SRCDIR   := Packages/Uintah/CCA/Components/HETransformation

SRCS     += $(SRCDIR)/NullBurn.cc $(SRCDIR)/SimpleBurn.cc \
	$(SRCDIR)/BurnFactory.cc  $(SRCDIR)/Burn.cc

PSELIBS	:= Packages/Uintah/Core/ProblemSpec \
	Packages/Uintah/Core/Exceptions

LIBS := $(XML_LIBRARY) -lm

include $(SRCTOP)/scripts/smallso_epilogue.mk