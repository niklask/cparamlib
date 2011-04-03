## Process this file with automake to produce Makefile.in
#
# Author: niklas@slac.stanford.edu
#
# Doxygen targets adopted from HippoDraw by Paul Kunz
#

INSTALLDIR = $(DESTDIR)$(datadir)/doc/$(PACKAGE)/html

EXTRA_DIST = Doxyfile.in header.html footer.html mainpage.dox install.dox tutorial.dox gamma_ind2.png gamma_pt_dist.png

all: docs

docs:
	if test -d html; then (cd html; rm -rf *.html ) fi
	doxygen Doxyfile

clean-local:
	if test -d html; then (rm -rf html ) fi

# will be used when doing "make install"
install-data-local: docs
	$(INSTALL) -d $(INSTALLDIR)
	$(INSTALL_DATA) html/* $(INSTALLDIR)

uninstall-local:
	rm -rf $(INSTALLDIR)

