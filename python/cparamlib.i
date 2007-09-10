/*
		cparamlib.i

		SWIG interface file for wrapping cparamlib into python

		$Source: /home/nkarlsson/usr/cvsroot/cparamlib/python/cparamlib.i,v $
		$Author: niklas $ $Date: 2007/09/10 21:18:05 $ $Revision: 1.2 $
*/

%module cparamlib

%{
#include "cparamlib.h"
%}

%include "cparamlib.h"

%include "carrays.i"
%array_class(double, doubleArray)

%include "cpointer.i"
%pointer_class(double, doublePtr)

#%include "cpointer.i"
#%pointer_class(PARAMSET, paramsetPtr)
