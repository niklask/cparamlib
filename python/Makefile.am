## Process this file with automake to produce Makefile.in
#
# Author: niklas@slac.stanford.edu
#

BUILT_SOURCES = $(srcdir)/cparammodel_wrap.c
SWIG_SOURCES = cparammodel.i

EXTRA_DIST = ParamModel.py __init__.py pycparam.py

pkgpython_PYTHON = ParamModel.py __init__.py pycparam.py
pkgpyexec_LTLIBRARIES = _pycparam.la
_pycparam_la_SOURCES = $(srcdir)/cparammodel_wrap.c $(SWIG_SOURCES)
_pycparam_la_CPPFLAGS = $(SWIG_PYTHON_CPPFLAGS) -I$(top_srcdir)/cparamlib
_pycparam_la_LDFLAGS = -module
_pycparam_la_LIBADD = $(top_srcdir)/cparamlib/libcparamlib.la

$(srcdir)/cparammodel_wrap.c : $(SWIG_SOURCES)
	$(SWIG) $(SWIG_PYTHON_OPT) -I$(top_srcdir)/cparamlib -o $@ $<
