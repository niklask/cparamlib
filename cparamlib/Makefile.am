## Process this file with automake to produce Makefile.in
#
# Author: niklas@slac.stanford.edu
#

cparamlibdir = $(includedir)/cparamlib

cparamlib_HEADERS = cparammodel.h
cparamlib_SOURCE = antinue.c antinumu.c elec.c gamma.c nue.c numu.c posi.c sigma.c

library_includedir=$(includedir)/cparamlib
library_include_HEADERS = $(h_sources)

INCLUDES = -I$(top_srcdir)

lib_LTLIBRARIES= libcparamlib.la
libcparamlib_la_SOURCES= $(cparamlib_HEADERS) $(cparamlib_SOURCE)
libcparamlib_la_LDFLAGS= -version-info 1.6.0 -release 1.6
