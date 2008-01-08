## Process this file with automake to produce Makefile.in
#
# Author: niklas@slac.stanford.edu
#

if DOCBUILD
DOCDIR = doc
else
DOCDIR = 
endif

if PYTHONBUILD
PYTHONDIR = python
else
PYTHONDIR = 
endif

SUBDIRS = cparamlib $(DOCDIR) $(PYTHONDIR)

EXTRA_DIST = examples/Makefile examples/params.c examples/examples.c

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = cparamlib.pc

# extra clean target for Sun CC compiler and generated distribution tarball
clean-local:
	-rm -rf SunWS_cache
	-rm -rf *.tar.gz

# extra uninstall for documentation dir
uninstall-local:
	-rm -rf $(docdir)
