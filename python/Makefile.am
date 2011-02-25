## Process this file with automake to produce Makefile.in
#
# Author: niklas@slac.stanford.edu
#

BUILT_SOURCES = $(srcdir)/cparamlib_wrap.c
SWIG_SOURCES = cparamlib.i

EXTRA_DIST = ParamModel.py __init__.py

$(srcdir)/cparamlib_wrap.c cparamlib.py: $(SWIG_SOURCES)
	$(SWIG) -python -I$(top_srcdir)/cparamlib -o $@ $<

pkgpyexec_PYTHON = ParamModel.py __init__.py cparamlib.py
pkgpyexec_LTLIBRARIES = _cparamlib.la
_cparamlib_la_SOURCES = $(BUILT_SOURCES) $(SWIG_SOURCES)
nodist__cparamlib_la_SOURCES = $(BUILT_SOURCES)
_cparamlib_la_CPPFLAGS = -I$(PYTHON_INCLUDE) -I$(top_srcdir)/cparamlib
_cparamlib_la_LDFLAGS = -module
_cparamlib_la_LIBADD = $(top_builddir)/cparamlib/libcparamlib.la

# make sure to delete generated code when with 'make clean'
clean-local:
	-rm -rf SunWS_cache
	-rm -rf cparamlib_wrap.c cparamlib.py
