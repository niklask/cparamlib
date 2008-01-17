## Process this file with automake to produce Makefile.in
#
# Author: niklas@slac.stanford.edu
#
# Doxygen targets adopted from HippoDraw by Paul Kunz
#

INSTALLDIR = $(DESTDIR)$(datadir)/doc/$(PACKAGE)/html

EXTRA_DIST = Doxyfile.in header.html footer.html mainpage.dox install.dox tutorial.dox gamma_ind2.png

# cleaned and installed
GEN_SUFFIX = html gif png css

all: docs

docs:
	if test -d html; then (cd html; rm -rf *.html ) fi
	doxygen Doxyfile

clean-local:
	if test -d html; then \
	  list='$(GEN_SUFFIX)'; for suffix in $$list; do \
	    (cd html; rm -rf *.$$suffix ) \
	  done \
	fi

# will be used when doing "make install"
install-data-local: docs
	$(INSTALL) -d $(INSTALLDIR)
	list='$(GEN_SUFFIX)'; for suffix in $$list; do \
	  (cd html; cp *.$$suffix $(INSTALLDIR)/. ) \
	done

# will be used when doing "make install-docs"
install-docs: docs
	$(INSTALL) -d $(INSTALLDIR)

uninstall-local:
	rm -rf $(INSTALLDIR)

