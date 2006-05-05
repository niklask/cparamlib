/*
		cparammodel.i

		SWIG interface file for wrapping cparamlib into python

		$Source: /home/nkarlsson/usr/cvsroot/cparamlib/Attic/cparammodel.i,v $
		$Author: niklas $ $Date: 2006/05/05 00:13:06 $ $Revision: 1.2 $
*/

%module cparammodel

%{
#include "cparammodel.h"
%}

%include "carrays.i"
%array_class(double, doubleArray)

%include "cpointer.i"
%pointer_functions(double, doublep)

%include "cparammodel.h"
