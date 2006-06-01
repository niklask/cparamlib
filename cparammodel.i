/*
		cparammodel.i

		SWIG interface file for wrapping cparamlib into python

		$Source: /home/nkarlsson/usr/cvsroot/cparamlib/Attic/cparammodel.i,v $
		$Author: niklas $ $Date: 2006/06/01 15:57:04 $ $Revision: 1.4 $
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
