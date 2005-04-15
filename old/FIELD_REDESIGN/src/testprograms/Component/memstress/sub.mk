#
# Makefile fragment for this subdirectory
# $Id$
#

SRCDIR := testprograms/Component/memstress

ifeq ($(LARGESOS),yes)
PSELIBS := SCICore
else
PSELIBS := Component/CIA Component/PIDL SCICore/Thread \
	SCICore/Exceptions SCICore/globus_threads
endif
LIBS := $(GLOBUS_LIBS) -lglobus_nexus -lglobus_dc -lglobus_common -lglobus_io

PROGRAM := $(SRCDIR)/memstress
SRCS := $(SRCDIR)/memstress.cc $(SRCDIR)/memstress_sidl.cc
GENHDRS := $(SRCDIR)/memstress_sidl.h

include $(SRCTOP)/scripts/program.mk

#
# $Log$
# Revision 1.4  2000/03/23 11:05:26  sparker
# Added *_sidl.h files to GENHDRS so that they will get built in time
#
# Revision 1.3  2000/03/21 06:13:36  sparker
# Added pattern rule for .sidl files
# Compile component testprograms
#
# Revision 1.2  2000/03/20 19:39:20  sparker
# Added VPATH support
#
# Revision 1.1  2000/03/17 09:31:06  sparker
# New makefile scheme: sub.mk instead of Makefile.in
# Use XML-based files for module repository
# Plus many other changes to make these two things work
#
#
