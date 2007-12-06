## Process this file with automake to produce Makefile.in
#
# Author: niklas@slac.stanford.edu
#

if PYTHONBUILD
PYTHONDIR = python
else
PYTHONDIR = 
endif

SUBDIRS = cparamlib $(PYTHONDIR)

EXTRA_DIST = examples/Makefile examples/params.c examples/test.c

# define where docs are installed
docdir = $(datadir)/doc/$(PACKAGE)

#
# Targets
#

# extra clean target for the generated distribution tarball
clean-local:
	-rm -rf *.tar.gz

# extra uninstall for
# - documentation dir
# - Sun CC compiler cache dir
uninstall-local:
	-rm -rf SunWS_cache
	-rm -rf $(docdir)
