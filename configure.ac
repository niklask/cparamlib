dnl Process this file with autoconf to produce a configure script
dnl
dnl	$Source: /home/nkarlsson/usr/cvsroot/cparamlib/configure.ac,v $
dnl $Author: niklas $ $Date: 2007/03/22 22:39:47 $ $Revision: 1.14 $
dnl
dnl Author: Niklas Karlsson <niklas@slac.stanford.edu>
dnl

AC_INIT(cparamlib/cparammodel.h)
AM_CONFIG_HEADER(config.h)
AC_CONFIG_AUX_DIR(config)
AM_INIT_AUTOMAKE(cparamlib, 1.7.5)

dnl Checks for programs
AC_PROG_CC
dnl disable static building
dnl AM_DISABLE_STATIC
AM_PROG_LIBTOOL

dnl Option to enable python extension build
AC_ARG_ENABLE(pythonbuild,
[  --enable-pythonbuild    build Python interface code (default = no)],
[case "${enableval}" in
  yes) pythonbuild=true ;;
  no)  pythonbuild=false ;;
  *) AC_MSG_ERROR(bad value ${enableval} for --enable-python) ;;
esac],[pythonbuild=false])

if test x$pythonbuild != xfalse; then
	AM_PATH_PYTHON(2.2)
	AC_PYTHON_DEVEL

	AC_PROG_SWIG(1.3)
	SWIG_PYTHON
fi

AM_CONDITIONAL(PYTHONBUILD, [test x$pythonbuild = xtrue])

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
cparamlib successfully configured.])
				  
if test x$pythonbuild != xfalse; then
  AC_MSG_RESULT( compile Python interface code: yes )
else
  AC_MSG_RESULT( compile Python interface code: no )
fi
