## Process this file with automake to produce Makefile.in
#
# Author: niklas@slac.stanford.edu
#

BUILT_SOURCES = $(srcdir)/cparamlib_wrap.c
SWIG_SOURCES = pycparam.i

EXTRA_DIST = ParamModel.py __init__.py

$(srcdir)/cparamlib_wrap.c pycparam.py: $(SWIG_SOURCES)
	$(SWIG) -python -I$(top_srcdir)/cparamlib -o $@ $<

if PYTHONBUILD
pkgpyexec_PYTHON = ParamModel.py __init__.py pycparam.py
pkgpyexec_LTLIBRARIES = _pycparam.la
_pycparam_la_SOURCES = $(srcdir)/cparamlib_wrap.c $(SWIG_SOURCES)
_pycparam_la_CPPFLAGS = -I$(PYTHON_INCLUDE) -I$(top_srcdir)/cparamlib
_pycparam_la_LDFLAGS = -version-info 1:0:0 -module
_pycparam_la_LIBADD = $(top_builddir)/cparamlib/libcparamlib.la
endif
