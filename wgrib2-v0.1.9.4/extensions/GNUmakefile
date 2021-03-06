#
# Driver Makefile
#

.SILENT:

SHELL = /bin/bash
PYTHON = python

ARCH := $(shell uname -s)
MACH := $(shell uname -m)
ifeq ($(ARCH),Darwin)
   OS := MacOSX
else
ifeq ($(ARCH),FreeBSD)
   OS := FreeBSD
else
   OS := $(shell uname -o)
endif
endif

bintype=# On windows this will be set to c, nc, nc4, dods, hdf

# Default installation directory depends on whether this is a standalone
# module or not
ifneq ($(wildcard ./include),)
   bindir := ./bin
else
   bindir := ../bin
endif
gexdir = $(bindir)/gex/$(bintype)

ifeq ($(OS),Cygwin)
	bintype=
endif

export bintype bindir gexdir

SUBDIRS = ams bjt fish gxyat hello ipc re shape shfilt gsf mf lats orb

all % :
	@ t=$@; argv="$(SUBDIRS)" ;\
          for d in $$argv; do                    \
            ( cd $$d                            ;\
              echo -n "Making $$t for extension <$$d> ... ";\
              $(MAKE) -e  $$t ;\
              echo ok ) \
          done

# Note: it may be better to have a top ut.py for this
check : install
	$(PYTHON) ut.all.py

install : force
	@/bin/mkdir -p $(bindir) $(gexdir)
	@ t=$@; argv="$(SUBDIRS)" ;\
          for d in $$argv; do                    \
            ( cd $$d                            ;\
              echo -n "Making $$t for extension <$$d> ... ";\
              $(MAKE) -e  $$t ;\
              echo ok ) \
          done
	@/bin/rm -rf $(gexdir)/udxt
	@/bin/cat $(gexdir)/*.udxt > $(gexdir)/udxt 

force : ;
