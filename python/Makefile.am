## Process this file with automake to produce Makefile.in
#
# Author: niklas@slac.stanford.edu
#

BUILT_SOURCES = $(srcdir)/cparamlib_wrap.c
SWIG_SOURCES = cparamlib.i

EXTRA_DIST = ParamModel.py __init__.py

$(srcdir)/cparamlib_wrap.c cparamlib.py: $(SWIG_SOURCES)
	$(SWIG) -python -I$(top_srcdir)/cparamlib -o $@ $<

if PYTHONBUILD
pkgpyexec_PYTHON = ParamModel.py __init__.py cparamlib.py
pkgpyexec_LTLIBRARIES = _cparamlib.la
_cparamlib_la_SOURCES = $(srcdir)/cparamlib_wrap.c $(SWIG_SOURCES)
_cparamlib_la_CPPFLAGS = -I$(PYTHON_INCLUDE) -I$(top_srcdir)/cparamlib
_cparamlib_la_LDFLAGS = -module
_cparamlib_la_LIBADD = $(top_builddir)/cparamlib/libcparamlib.la
endif
