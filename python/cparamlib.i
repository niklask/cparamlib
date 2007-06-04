/*
		cparamlib.i

		SWIG interface file for wrapping cparamlib into python

		$Source: /home/nkarlsson/usr/cvsroot/cparamlib/python/cparamlib.i,v $
		$Author: niklas $ $Date: 2007/06/04 22:22:23 $ $Revision: 1.1 $
*/

%module cparamlib

%{
#include "cparammodel.h"
%}

%include "cparammodel.h"

%include "carrays.i"
%array_class(double, doubleArray)

%include "cpointer.i"
%pointer_class(double, doublePtr)

#%include "cpointer.i"
#%pointer_class(PARAMSET, paramsetPtr)
