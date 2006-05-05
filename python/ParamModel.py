#
# Python package implementing the parametric model described in ApJ paper
#
# $Source: /home/nkarlsson/usr/cvsroot/cparamlib/python/ParamModel.py,v $
# $Author: niklas $ $Date: 2006/05/05 20:40:50 $ $Revision: 1.2 $
#

#
# Provides a class that implements calculations of cross sections based on the
# parametric model published in ApJ
#

import sys
from numarray import *

import cparammodel

ID_GAMMA = 0
ID_ELECTRON = 1
ID_POSITRON = 2
ID_NUE = 3
ID_NUMU = 4
ID_ANTINUE = 5
ID_ANTINUMU = 6

L_MAX = [0.96, 0.96, 0.94, 0.98, 0.98, 0.94, 0.98]
W_NDL = [15.0, 20.0, 15.0, 15.0, 15.0, 20.0, 15.0]
W_NDH = [44.0, 45.0, 47.0, 42.0, 40.0, 45.0, 40.0]

# --------------------------------------------------------------------------------
#  Class definitions
# --------------------------------------------------------------------------------

# Description:
# The class ParamModel is a presentation of the parametric model for a given secondary particle and a
# set of proton kinetic energies. On instantiation, the particle species is set and tables for
# parameter values are precalculated. This is done using the SWIG wrapper for libcparammodel. Once
# the class has been instantiated, one can calculate inclusive cross sections for any set of secondary
# particles. Note that those methods also take Tp as argument and this Tp must be the same as used
# when instantiating.

#
# Class ParamModel
#
class ParamModel:
    #
    # Constructor
    #
    # Tp should be a numarray containing the proton kinetic energies (Tp) for which
    # the parameter table is calculated
    #
    def __init__(self, Tp, particle = None):
        if (particle == None):
            # default to gamma-rays
            self.particle = ID_GAMMA
        else:
            self.particle = particle
            
        self.__paramTableND = self.paramND(Tp)
        self.__paramTableDiff = self.paramDiff(Tp)
        self.__paramTableDelta = self.paramDelta(Tp)
        self.__paramTableRes = self.paramRes(Tp)

    #
    # Method reinit
    #
    # Recalculates the parameter tables and sets new particle species
    #
    def reinit(self, Tp, particle = None):
        if (particle == None):
            # default to gamma-rays
            self.particle = ID_GAMMA
        else:
            self.particle = particle
            
        self.__paramTableND = self.paramND(Tp)
        self.__paramTableDiff = self.paramDiff(Tp)
        self.__paramTableDelta = self.paramDelta(Tp)
        self.__paramTableRes = self.paramRes(Tp)

    #
    # Method paramND
    #
    # Create a table with the non-diffraction parameter values for all Tp:s
    # This is a wrapper around <particle>_param_nd(Tp, a)
    #
    def paramND(self, Tp):
        n = Tp.shape[0]
        
        paramTable = zeros((n,9), Float32)

        a = cparammodel.doubleArray(9)
        for i in range(n):
            if (self.particle == ID_GAMMA):
                cparammodel.gamma_param_nd(Tp[i], a)
            if (self.particle == ID_ELECTRON):
                cparammodel.elec_param_nd(Tp[i], a)
            if (self.particle == ID_POSITRON):
                cparammodel.posi_param_nd(Tp[i], a)
            if (self.particle == ID_NUE):
                cparammodel.nue_param_nd(Tp[i], a)
            if (self.particle == ID_NUMU):
                cparammodel.numu_param_nd(Tp[i], a)
            if (self.particle == ID_ANTINUE):
                cparammodel.nue_param_nd(Tp[i], a)
            if (self.particle == ID_ANTINUMU):
                cparammodel.numu_param_nd(Tp[i], a)

            for j in range(9):
                paramTable[i,j] = a[j]
                
        return paramTable

    #
    # Method paramDiff
    #
    # Create a table with the diffraction parameter values for all Tp:s
    # This is a wrapper around <particle>_param_diff(Tp, a)
    #
    def paramDiff(self, Tp):
        n = Tp.shape[0]
        
        paramTable = zeros((n,8), Float32)

        b = cparammodel.doubleArray(8)
        for i in range(n):
            if (self.particle == ID_GAMMA):
                cparammodel.gamma_param_diff(Tp[i], b)
            if (self.particle == ID_ELECTRON):
                cparammodel.elec_param_diff(Tp[i], b)
            if (self.particle == ID_POSITRON):
                cparammodel.posi_param_diff(Tp[i], b)
            if (self.particle == ID_NUE):
                cparammodel.nue_param_diff(Tp[i], b)
            if (self.particle == ID_NUMU):
                cparammodel.numu_param_diff(Tp[i], b)
            if (self.particle == ID_ANTINUE):
                cparammodel.nue_param_diff(Tp[i], b)
            if (self.particle == ID_ANTINUMU):
                cparammodel.numu_param_diff(Tp[i], b)

            for j in range(8):
                paramTable[i,j] = b[j]
                
        return paramTable

    #
    # Method paramDelta
    #
    # Create a table with the delta(1232) parameter values for all Tp:s
    # This is a wrapper around <particle>_param_delta(Tp, a)
    #
    def paramDelta(self, Tp):
        n = Tp.shape[0]
        
        paramTable = zeros((n,5), Float32)

        c = cparammodel.doubleArray(5)
        for i in range(n):
            if (self.particle == ID_GAMMA):
                cparammodel.gamma_param_delta(Tp[i], c)
            if (self.particle == ID_ELECTRON):
                cparammodel.elec_param_delta(Tp[i], c)
            if (self.particle == ID_POSITRON):
                cparammodel.posi_param_delta(Tp[i], c)
            if (self.particle == ID_NUE):
                cparammodel.nue_param_delta(Tp[i], c)
            if (self.particle == ID_NUMU):
                cparammodel.numu_param_delta(Tp[i], c)
            if (self.particle == ID_ANTINUE):
                cparammodel.nue_param_delta(Tp[i], c)
            if (self.particle == ID_ANTINUMU):
                cparammodel.numu_param_delta(Tp[i], c)

            for j in range(5):
                paramTable[i,j] = c[j]
                
        return paramTable
    
    #
    # Method paramRes
    #
    # Create a table with the res(1600) parameter values for all Tp:s
    # This is a wrapper around <particle>_param_res(Tp, a)
    #
    def paramRes(self, Tp):
        n = Tp.shape[0]
        
        paramTable = zeros((n,5), Float32)

        d = cparammodel.doubleArray(5)
        for i in range(n):
            if (self.particle == ID_GAMMA):
                cparammodel.gamma_param_res(Tp[i], d)
            if (self.particle == ID_ELECTRON):
                cparammodel.elec_param_res(Tp[i], d)
            if (self.particle == ID_POSITRON):
                cparammodel.posi_param_res(Tp[i], d)
            if (self.particle == ID_NUE):
                cparammodel.nue_param_res(Tp[i], d)
            if (self.particle == ID_NUMU):
                cparammodel.numu_param_res(Tp[i], d)
            if (self.particle == ID_ANTINUE):
                cparammodel.nue_param_res(Tp[i], d)
            if (self.particle == ID_ANTINUMU):
                cparammodel.numu_param_res(Tp[i], d)

            for j in range(5):
                paramTable[i,j] = d[j]
                
        return paramTable

    #
    # Method sigmaND
    #
    # Calculates the inclusive non-diffraction cross section for given E and Tp values
    # E and Tp are numarrays, returns a numarray
    #
    def sigmaND(self, E, Tp):
        n = E.shape[0]
        m = Tp.shape[0]

        sigma = zeros((n,m), Float32)
        
				# init some variables, given in table 2
        Lmin = -2.6
        Lmax = L_MAX[self.particle]*log10(Tp)
        Wl = W_NDL[self.particle]
        Wh = W_NDH[self.particle]
        
				# calculate the log of E and Tp
        x = log10(E)
        y = log10(Tp*0.001)

        a = self.__paramTableND
        
        cutoff = (1.0/(1.0 + exp(Wl*(Lmin - x))))*(1.0/(1.0 + exp(Wh*(x - Lmax))))

        #
        # calculate inclusive cross section
        #
        for i in range(m):
            s = a[i,0]*exp(-a[i,1]*(x - a[i,3] + a[i,2]*(x - a[i,3])**2)**2) + a[i,4]*exp(-a[i,5]*(x - a[i,8] + a[i,6]*(x - a[i,8])**2 + a[i,7]*(x - a[i,8])**3)**2)
            # multiply with cutoff
            s = s*cutoff
            # if any element is less than zero then set equal to zero
            choose(less(s, 0.0), (s, 0.0))
            sigma[:,i] = s

        #
        # do the renormalization
        #
        if (self.particle == ID_GAMMA):
            renorm = 3.05*exp(-107.0*((y + 3.25)/(1.0 + 8.08*(y + 3.25)))**2)
            renorm = choose(less_equal(Tp, 1.95), (1.01, renorm))
        elif (self.particle == ID_ELECTRON):
            renorm = 3.63*exp(-106*((y + 3.26)/(1.0 + 9.21*(y + 3.26)))**2) - 0.182*y - 0.175*y**2
            renorm = choose(less_equal(Tp, 15.6), (1.01, renorm))
        elif (self.particle == ID_POSITRON):
            renorm = 2.22*exp(-98.9*((y + 3.25)/(1.0 + 1.04*(y + 3.25)))**2)
            renorm = choose(less_equal(Tp, 5.52), (1.0, renorm))
        elif (self.particle == ID_NUE):
            renorm = 0.329*exp(-249*((y + 3.26)/(1.0 + 6.56*(y + 3.26)))**2) - 9.57*y - 0.229*y**2
            renorm = choose(less_equal(Tp, 7.81), (1.0, renorm))
        elif (self.particle == ID_NUMU):
            renorm = 2.23*exp(-93.4*((y + 3.25)/(1.0 + 8.38*(y + 3.25)))**2) - 0.376*y - 0.121*y**2
            renorm = choose(less_equal(Tp, 15.6), (1.0, renorm))
        elif (self.particle == ID_ANTINUME):
            renorm = 2.67*exp(-45.7*((y + 3.27)/(1.0 + 6.59*(y + 3.27)))**2) - 0.301*y - 0.208*y**2
            renorm = choose(less_equal(Tp, 15.6), (1.0, renorm))    
        elif (self.particle == ID_ANTINUMU):
            renorm = 2.56*exp(-107*((y + 3.25)/(1.0 + 8.34*(y + 3.25)))**2) - 0.385*y - 0.125*y**2
            renorm = choose(less_equal(Tp, 15.6), (1.0, renorm))
            
        for i in range(n):
            sigma[i,:] = sigma[i,:]*renorm

        return sigma
            
    #
    # Method sigmaDiff
    #
    # Calculates the inclusive diffraction cross section for given E and Tp values
    # E and Tp are numarrays, returns a numarray
    #
    def sigmaDiff(self, E, Tp):
        n = E.shape[0]
        m = Tp.shape[0]

        sigma = zeros((n,m), Float32)

				# init some variables, given in table 2
        Lmax = log10(Tp)
        Wdiff = 75.0

				# calculate the log of E and Tp
        x = log10(E)
        y = log10(Tp*0.001)

        b = self.__paramTableDiff
        
        cutoff = 1.0/(1.0 + exp(Wdiff*(x - Lmax)))

        #
        # calculate inclusive cross section
        #
        for i in range(m):
            sigma = b[i,0]*exp(-b[i,1]*((x - b[i,2])/(1.0 + b[i,3]*(x - b[i,2])))**2) + b[i,4]*exp(-b[i,5]*((x - b[i,6])/(1.0 + b[i,7]*(x - b[i,6])))**2)
            # multiply with cutoff
            s = s*cutoff
            # if any element is less than zero then set equal to zero
            choose(less(s, 0.0), (s, 0.0))
            sigma[:,i] = s

        return sigma

    #
    # Method sigmaDelta
    #
    # Calculates the inclusive delta(1232) cross section for given E and Tp values
    # E and Tp are numarrays, returns a numarray
    #
    def sigmaDelta(self, E, Tp):
        n = E.shape[0]
        m = Tp.shape[0]

        sigma = zeros((n,m), Float32)

				# init some variables, given in table 2
        Lmax = log10(Tp)
        Wdelta = 75.0

				# calculate the log of E and Tp
        x = log10(E)
        y = log10(Tp*0.001)

        c = self.__paramTableDelta
        
        cutoff = 1.0/(1.0 + exp(Wdelta*(x - Lmax)))

        #
        # calculate inclusive cross section
        #
        for i in range(m):
            sigma = c[i,0]*exp(-c[i,1]*((x - c[i,2])/(1.0 + c[i,3]*(x - c[i,2]) + c[i,4]*(x - c[i,2])**2))**2)
            # multiply with cutoff
            s = s*cutoff
            # if any element is less than zero then set equal to zero
            choose(less(s, 0.0), (s, 0.0))
            sigma[:,i] = s

        return sigma

    #
    # Method sigmaRes
    #
    # Calculates the inclusive res(1600) cross section for given E and Tp values
    # E and Tp are numarrays, returns a numarray
    #
    def sigmaDelta(self, E, Tp):
        n = E.shape[0]
        m = Tp.shape[0]

        sigma = zeros((n,m), Float32)

				# init some variables, given in table 2
        Lmax = log10(Tp)
        Wres = 75.0

				# calculate the log of E and Tp
        x = log10(E)
        y = log10(Tp*0.001)

        d = self.__paramTableRes
        
        cutoff = 1.0/(1.0 + exp(Wres*(x - Lmax)))

        #
        # calculate inclusive cross section
        #
        for i in range(m):
            sigma = d[i,0]*exp(-d[i,1]*((x - d[i,2])/(1.0 + d[i,3]*(x - d[i,2]) + d[i,4]*(x - d[i,2])**2))**2)
            # multiply with cutoff
            s = s*cutoff
            # if any element is less than zero then set equal to zero
            choose(less(s, 0.0), (s, 0.0))
            sigma[:,i] = s

        return sigma
        
if (__name__ == "__main__"):
    sys.exit(0)
