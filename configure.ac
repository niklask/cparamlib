dnl Process this file with autoconf to produce a configure script
dnl
dnl	$Source: /home/nkarlsson/usr/cvsroot/cparamlib/configure.ac,v $
dnl $Author: niklas $ $Date: 2007/12/07 00:24:11 $ $Revision: 1.24 $
dnl
dnl Author: Niklas Karlsson <niklas@slac.stanford.edu>
dnl

AC_INIT(cparamlib/cparamlib.h)
AM_CONFIG_HEADER(config.h)
AC_CONFIG_AUX_DIR(config)
AM_INIT_AUTOMAKE(cparamlib, 2.0.0)

dnl Checks for programs
AC_PROG_CC
dnl disable static building
dnl AM_DISABLE_STATIC
AM_PROG_LIBTOOL

dnl Option to enable python extension build
AC_ARG_ENABLE(python,
[  --enable-python    build Python interface code (default = no)],
[case "${enableval}" in
  yes) pythonbuild=true ;;
  no)  pythonbuild=false ;;
  *) AC_MSG_ERROR(bad value ${enableval} for --enable-python) ;;
esac],[pythonbuild=false])

dnl If --enable-python then make sure Python is installed with version >= 2.2
dnl and set variables $(PYTHON_INCLUDE) and $(PYTHON_LIB)
if test x$pythonbuild != xfalse; then
  AM_PATH_PYTHON(2.2)

	dnl The following code adopted from HippoDraw by Paul F. Kunz
	python_include_path=`$PYTHON -c "from distutils import sysconfig; print sysconfig.get_python_inc()"`
	python_lib_path=`$PYTHON -c "from distutils import sysconfig; print sysconfig.get_python_lib(standard_lib=True)"`

  AC_SUBST([PYTHON_INCLUDE],[$python_include_path] )
  AC_SUBST([PYTHON_LIB],[$python_lib_path])
fi
	
AM_CONDITIONAL(PYTHONBUILD, [test x$pythonbuild = xtrue])

dnl Make sure SWIG is installed and set variable $(SWIG)
dnl Only warn "not found" if --enable-python was used
dnl (Do not care about version number for now)
if test x$pythonbuild != xfalse; then
  AC_PATH_PROG([SWIG],[swig])
  if test -z "$SWIG" ; then
	  AC_MSG_WARN([cannot find the 'swig' executable, please install it])
      
  fi
fi

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
  AC_MSG_RESULT(  compile Python interface code: yes )
else
  AC_MSG_RESULT(  compile Python interface code: no )
fi

