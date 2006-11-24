/*
		cparammodel.i

		SWIG interface file for wrapping cparamlib into python

		$Source: /home/nkarlsson/usr/cvsroot/cparamlib/python/Attic/cparammodel.i,v $
		$Author: niklas $ $Date: 2006/11/24 23:31:31 $ $Revision: 1.1 $
*/

%module cparammodel

%{
#include "cparammodel.h"
%}

%include "carrays.i"
%array_class(double, doubleArray)

%include "cpointer.i"
%pointer_functions(double, doublePtr)

%include "cparammodel.h"
