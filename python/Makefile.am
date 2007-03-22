## Process this file with automake to produce Makefile.in
#
# Author: niklas@slac.stanford.edu
#

BUILT_SOURCES = $(srcdir)/cparamlib_wrap.c
SWIG_SOURCES = pycparam.i

EXTRA_DIST = ParamModel.py __init__.py

$(srcdir)/cparamlib_wrap.c pycparam.py: $(SWIG_SOURCES)
	$(SWIG) $(SWIG_PYTHON_OPT) -I$(top_srcdir)/cparamlib -o $@ $<

# don't know if this is the best way of doing this, but having the if statement
# makes the distribution work
if PYTHONBUILD
pkgpython_PYTHON = ParamModel.py __init__.py pycparam.py
pkgpyexec_LTLIBRARIES = _pycparam.la
_pycparam_la_SOURCES = $(srcdir)/cparamlib_wrap.c $(SWIG_SOURCES)
_pycparam_la_CPPFLAGS = $(SWIG_PYTHON_CPPFLAGS) -I$(top_srcdir)/cparamlib
_pycparam_la_LDFLAGS = -module
_pycparam_la_LIBADD = $(top_builddir)/cparamlib/libcparamlib.la
endif
