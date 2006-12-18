dnl Process this file with autoconf to produce a configure script
dnl
dnl	$Source: /home/nkarlsson/usr/cvsroot/cparamlib/configure.ac,v $
dnl $Author: niklas $ $Date: 2006/12/18 17:41:57 $ $Revision: 1.5 $
dnl
dnl Author: Niklas Karlsson <niklas@slac.stanford.edu>
dnl

AC_INIT(cparamlib/cparammodel.h)
AM_CONFIG_HEADER(config.h)
AC_CONFIG_AUX_DIR(config)
AM_INIT_AUTOMAKE(cparamlib, 1.7.1)

dnl Checks for programs
AC_PROG_CC
dnl disable static building
dnl AM_DISABLE_STATIC
AM_PROG_LIBTOOL

dnl
dnl Some user defined macros
dnl

dnl CPL_ARG_ENABLE(enable-option, help-string, default-value)
AC_DEFUN([CPL_ARG_ENABLE],[
  AC_ARG_ENABLE($1,
    AS_HELP_STRING([--enable-$1], [$2]m4_if($3,,[],[ @<:@default=$3@:>@])),
    [enable_[]translit($1, [-], [_])=$enableval],
    [enable_[]translit($1, [-], [_])=$3])])

dnl List --enable-python as option
CPL_ARG_ENABLE(python, [compile Python interface code], no)

dnl Test if --enable-python is set
if test "$enable_python" = "yes"; then
    ### Redefine AC_MSG_ERROR as AC_MSG_WARN around the call to
    ### AM_PATH_PYTHON, such that if python is not found, configure will
    ### not stop.
    pushdef([AC_MSG_ERROR],defn([AC_MSG_WARN]))
    AM_PATH_PYTHON(2.2)
    AC_PYTHON_DEVEL
    popdef([AC_MSG_ERROR])
fi

AC_SUBST(PYTHON_INSTDIR)

AC_SUBST(enable_python)

AM_CONDITIONAL(enable_python, [test "$enable_python" = yes])

dnl swig normally only used by python and java so far, but always do the
dnl AC_PROG_SWIG check since swig also used to generate files for distribution.
pushdef([AC_MSG_WARN],[true])
AC_PROG_SWIG(1.3.21)
SWIG_PYTHON
popdef([AC_MSG_WARN])

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
				  
AC_MSG_RESULT([  compile Python interface code: $enable_python
])

