/*
		pycparam.i

		SWIG interface file for wrapping cparamlib into python

		$Source: /home/nkarlsson/usr/cvsroot/cparamlib/python/Attic/pycparam.i,v $
		$Author: niklas $ $Date: 2006/11/25 02:29:29 $ $Revision: 1.1 $
*/

%module pycparam

%{
#include "cparammodel.h"
%}

%include "carrays.i"
%array_class(double, doubleArray)

%include "cpointer.i"
%pointer_functions(double, doublePtr)

%include "cparammodel.h"
