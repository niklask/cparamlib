dnl Process this file with autoconf to produce a configure script
dnl
dnl	$Source: /home/nkarlsson/usr/cvsroot/cparamlib/configure.ac,v $
dnl $Author: nkarlsson $ $Date: 2009/02/13 21:27:03 $ $Revision: 1.39 $
dnl
dnl Author: Niklas Karlsson <niklas@slac.stanford.edu>
dnl

AC_INIT([cparamlib], [2.1.2])
AM_CONFIG_HEADER(config.h)
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_SRCDIR([cparamlib/cparamlib.h])
AM_INIT_AUTOMAKE
dnl Do not forget to change library interface version in cparamlib/Makefile.am!

dnl Checks for programs
AC_PROG_CC
dnl disable static building
AM_DISABLE_STATIC
AM_PROG_LIBTOOL

dnl Check for gmake and use it
AC_CHECK_PROG(GMAKE, gmake, gmake, make)
AC_SUBST(GMAKE)

dnl Check for header files stdio.h and math.h
AC_CHECK_HEADER(stdio.h, [],
                [AC_MSG_WARN(missing required header file <stdio.h>)])

AC_CHECK_HEADER(math.h, [],
                [AC_MSG_WARN(missing required header file <math.h>)])

dnl Turn on -Wall if compiling with gcc
#AC_COMPILE_WARNINGS

dnl
dnl Option to enable python extension build
dnl
AC_ARG_ENABLE(python,
AS_HELP_STRING([--enable-python],[build Pythin interface code [[default=no]]]),
[case "${enableval}" in
  yes) pythonbuild=true ;;
  no)  pythonbuild=false ;;
  *) AC_MSG_ERROR(bad value ${enableval} for --enable-python) ;;
esac],[pythonbuild=false])

dnl Make sure Python is installed with version >= 2.2 and set variables
dnl   $(PYTHON_INCLUDE)
dnl   $(PYTHON_LIB)
dnl Do this even if --disable-python because it sets variables required when
dnl doing a 'make distcheck' 
AM_PATH_PYTHON(2.2)

dnl The following code adopted from HippoDraw by Paul F. Kunz
dnl python_include_path=`$PYTHON -c "from distutils import sysconfig; print sysconfig.get_python_inc()"`
dnl python_lib_path=`$PYTHON -c "from distutils import sysconfig; print sysconfig.get_python_lib(standard_lib=True)"`
dnl Use the above two lines instead if using Python 2.x
python_include_path=`$PYTHON -c "import sysconfig; print(sysconfig.get_config_vars()[['INCLUDEPY']])"`
python_lib_path=`$PYTHON -c "import sysconfig; print(sysconfig.get_config_vars()[['LIBDEST']])"`

AC_SUBST([PYTHON_INCLUDE],[$python_include_path] )
AC_SUBST([PYTHON_LIB],[$python_lib_path])

dnl Make sure SWIG is installed and set variable $(SWIG)
dnl (do not care about version number for now)
dnl Do this even if --disable-python because it sets a variable required when
dnl doing a 'make distcheck' 
AC_PATH_PROG([SWIG],[swig])
if test -z "$SWIG" ; then
   AC_MSG_WARN([cannot find the 'swig' executable, please install it])
fi
	
AM_CONDITIONAL(PYTHONBUILD, [test x$pythonbuild = xtrue])

dnl
dnl Option to build and install documentation
dnl
AC_ARG_ENABLE(userdoc,
AS_HELP_STRING([--enable-userdoc],[generate documentation [[default=no]]]),
[case "${enableval}" in
  yes) docbuild=true ;;
  no)  docbuild=false ;;
  *)   AC_MSG_ERROR(bad value ${enableval} for --enable-userdoc) ;;
esac],[docbuild=false])

dnl Make sure Doxygen is installed and set variable $(DOXYGEN)
dnl only warn "not found" if actually making documentation
if test x$docbuild != xfalse; then
  AC_PATH_PROG([DOXYGEN],[doxygen])
  if test -z "$DOXYGEN" ; then
    AC_MSG_WARN([cannot find the 'doxygen' executable, please install it])
    docbuild=false
  fi
fi

AM_CONDITIONAL(DOCBUILD, [test x$docbuild = xtrue])

AM_SANITY_CHECK

dnl
dnl Files to be generated 
dnl
AC_OUTPUT(Makefile \
          cparamlib.pc \
          cparamlib/Makefile \
		  python/Makefile \
          doc/Makefile \
          doc/Doxyfile
)

dnl
dnl Output message about what is going to be compiled/installed
dnl
AC_MSG_RESULT([
cparamlib successfully configured.])
				  
if test x$pythonbuild != xfalse; then
  AC_MSG_RESULT([  compile Python interface code: yes ])
else
  AC_MSG_RESULT([  compile Python interface code: no ])
fi

if test x$docbuild != xfalse; then
  AC_MSG_RESULT([  generate Doxygen documentation: yes ])
else
  AC_MSG_RESULT([  generate Doxygen documentation: no ])
fi

