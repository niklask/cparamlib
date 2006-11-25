dnl Process this file with autoconf to produce a configure script
dnl
dnl Author: Niklas Karlsson <niklas@slac.stanford.edu>
dnl

AC_INIT(cparamlib/cparammodel.h)
AM_CONFIG_HEADER(config.h)
AC_CONFIG_AUX_DIR(config)
AM_INIT_AUTOMAKE(cparamlib, 1_6)

dnl Checks for programs
AC_PROG_CC
dnl disable static building
dnl AM_DISABLE_STATIC
AM_PROG_LIBTOOL
AM_PATH_PYTHON(2.3)
AC_PROG_SWIG(1.3.21)
SWIG_PYTHON
AM_SANITY_CHECK

dnl Turn on -Wall if compiling with gcc
#AC_COMPILE_WARNINGS

dnl Check for gmake and use it
AC_CHECK_PROG(GMAKE, gmake, gmake, make)
AC_SUBST(GMAKE)

dnl files to be generated 
AC_OUTPUT(Makefile \
          cparamlib/Makefile \
					python/Makefile
)

AC_MSG_RESULT([
    cparamlib project now configured.
])

