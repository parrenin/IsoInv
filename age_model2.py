"""
This is the age_model.py module, to invert isochronal layers along a radar profile.
"""

import os
import time
import sys
import random
import math as m
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LogNorm, Normalize
from matplotlib.backends.backend_pdf import PdfPages
from scipy.interpolate import interp1d  #Deprecated, use numpy.interp
from scipy.optimize import leastsq
from scipy.special import erf
import yaml
import pandas as pd

###Registration of start time
START_TIME = time.time()



#Physical constants
Kg0 = 9.828
Kg1 = -5.7*10**-3
Lf = 333.5*10**3
rhog = 917.
cg0 = 152.5
cg1 = 7.122
ggrav = 9.81
#Tf0 = 273.16          #This is the value from Cuffey and Paterson (2010)
#Tf1 = -9.8e-8          #This is the value from Cuffey and Paterson (2010)
Tf0 = 273.16-0.024         #This is the value from Catherine Ritz's thesis
Tf1 = -7.4e-8         #This is the value from Catherine Ritz's thesis

def cg(T):
#    return 2097.*np.ones(np.shape(T))    #at 0 degC
    return cg0+cg1*T

def Kg(T, D):
    """From Patterson eqs. 9.2 and 9.4"""
#    return 2.10*np.ones(np.shape(T))     #at 0 degC
#    return Kg0*np.exp(Kg1*T)
    KiT = Kg0*np.exp(Kg1*T)
    return (2.*KiT*D)/(3.-D)

def Tf(P):
    return Tf0+Tf1*P

def interp1d_stair_aver(x, y):   #TODO: deal with the case x not sorted
    """
    Interpolation of a staircase function using averaging.
    This function returns nan outside of the input abscissa range.
    """
    def f(xp):
        yp = np.empty(np.size(xp)-1)
        xmod = x[~(np.isnan(x)+np.isnan(y))]
        ymod = y[~(np.isnan(x)+np.isnan(y))]
        yint = np.cumsum(np.concatenate((np.array([0]), ymod[:-1]*(xmod[1:]-xmod[:-1]))))
        g = interp1d(xmod, yint, bounds_error=False, fill_value=np.nan)
        yp = np.where((xp[:-1] > min(xmod))*(xp[1:] < max(xmod)), (g(xp[1:])-g(xp[:-1]))/\
             (xp[1:]-xp[:-1]), np.nan) #Maybe this is suboptimal since we compute twice g(xp[i])
        return yp
    return f

def interp1d_stair_aver_withnan(x, y):   #TODO: deal with the case x not sorted
    """
    Interpolation of a staircase function using averaging.
    This function returns nan when there are all nans in one interpolation interval.
    """
    def f(xp):
        xmod = x[~(np.isnan(x)+np.isnan(y))]
        ymod = y[~(np.isnan(x)+np.isnan(y))]
        yp = np.empty(np.size(xp)-1)
        yint = np.cumsum(np.concatenate((np.array([0]), ymod[:-1]*(xmod[1:]-xmod[:-1]))))
        g = interp1d(xmod, yint, bounds_error=False, fill_value=np.nan)
        yp = np.where((xp[:-1] > min(xmod))*(xp[1:] < max(xmod)), (g(xp[1:])-g(xp[:-1]))/\
             (xp[1:]-xp[:-1]), np.nan) #Maybe this is suboptimal since we compute twice g(xp[i])
        for i in range(np.size(xp)-1):
            if np.isnan(y[np.where((x >= xp[i])*(x < xp[i+1]))]).all():
                yp[i] = np.nan
                return yp

    return f

def interp1d_lin_aver_withoutnan(x, y):
    """
    Interpolation of a linear by parts function using averaging.
    This function returns nan when there are all nans in one interpolation interval.
    FIXME: there is a problem in this routine when the x are in decreasing order.
    """
    def f(xp):
        yp = np.empty(np.size(xp)-1)
        for i in range(np.size(xp)-1):
            xmod = x[~(np.isnan(x)+np.isnan(y))]
            ymod = y[~(np.isnan(x)+np.isnan(y))]
            xmod2 = xmod[np.where((xmod > xp[i])*(xmod < xp[i+1]))]
            ymod2 = ymod[np.where((xmod > xp[i])*(xmod < xp[i+1]))]
            xmod3 = np.concatenate((np.array([xp[i]]), xmod2, np.array([xp[i+1]])))
            g = interp1d(xmod, ymod, bounds_error=False, fill_value=np.nan)
            ymod3 = np.concatenate((np.array([g(xp[i])]), ymod2, np.array([g(xp[i+1])])))
            if np.isnan(ymod3).all():
                yp[i] = np.nan
            else:
                xmod4 = xmod3[np.where(~(np.isnan(ymod3)+np.isnan(xmod3)))]
                ymod4 = ymod3[np.where(~(np.isnan(ymod3)+np.isnan(xmod3)))]
                yp[i] = np.sum((ymod4[1:]+ymod4[:-1])/2*(xmod4[1:]-xmod4[:-1]))
                yp[i] = yp[i]/(xmod4[-1]-xmod4[0])
        return yp
    return f


def interp1d_lin_aver(x, y):
    """
    Interpolation of a linear by parts function using averaging.
    This function returns nan when there are all nans in one interpolation interval.
    FIXME: there is a problem in this routine when the x are in decreasing order.
    """
    def f(xp):
        yp = np.empty(np.size(xp)-1)
        for i in range(np.size(xp)-1):
            xmod = x[~(np.isnan(x)+np.isnan(y))]
            ymod = y[~(np.isnan(x)+np.isnan(y))]
            xmod2 = xmod[np.where((xmod > xp[i])*(xmod < xp[i+1]))]
            ymod2 = ymod[np.where((xmod > xp[i])*(xmod < xp[i+1]))]
            xmod3 = np.concatenate((np.array([xp[i]]), xmod2, np.array([xp[i+1]])))
            g = interp1d(x, y, bounds_error=False, fill_value=np.nan)
            ymod3 = np.concatenate((np.array([g(xp[i])]), ymod2, np.array([g(xp[i+1])])))
            if np.isnan(ymod3).all():
                yp[i] = np.nan
            else:
                xmod4 = xmod3[np.where(~(np.isnan(ymod3)+np.isnan(xmod3)))]
                ymod4 = ymod3[np.where(~(np.isnan(ymod3)+np.isnan(xmod3)))]
                yp[i] = np.sum((ymod4[1:]+ymod4[:-1])/2*(xmod4[1:]-xmod4[:-1]))
                yp[i] = yp[i]/(xmod4[-1]-xmod4[0])
        return yp
    return f


class RadarLine:

    def __init__(self, label):
        self.label = label

    def init(self):
        self.is_bedelev = False
        self.is_trace = False
        self.nbdsz = 0
        self.calc_sigma = True
        self.invert_G0 = False
        self.settick = 'auto'
        self.interp_method = 'lin_aver'
        self.distance_unit = 'km'
        self.nbhor = 0
        self.nbiso = 0
        self.firn_correction = 14.6
        self.resolution = 1.
        self.is_EDC = False
        self.calc_isoage = False
        self.distance_EDC = 0.
        self.age_surf = -50.
        self.dzeta = 0.01
        self.tm_iter = 5
        self.Ts = 212.74
        self.p_prior = 2
        self.p_sigma = 5
        self.G0_prior = 0.051
        self.G0_sigma = 0.025
        self.EDC_line_dashed = False
        self.is_NESW = False
        self.reverse_distance = False
        self.aspect = 0.028
        self.is_legend = True
        self.min_tick = 0.
        self.max_tick = 100.
        self.delta_tick = 10.
        self.invert_thk = False
        self.accu_min = 1.5
        self.accu_max = 2.5
        self.p_min = 0.
        self.p_max = 20.
        self.G0_min = 40.
        self.G0_max = 80.
        self.melting_min = 0.
        self.melting_max = 1.
        self.age_min = 0.5
        self.age_max = 3.
        self.height_min = 0.
        self.height_max = 200.
        self.reso_min = 5.
        self.reso_max = 30.
        self.opt_method = 'MH1D'
        self.MHnbiter = 3000
        self.MHiter_adapt1 = 1000
        self.MHiter_adapt2 = -1
        self.accu_step = 0.001
        self.p_step = 0.5
        self.G0_step = 0.002
        self.thick_step = 10.

        #definition of some global parameters
#        exec(open(self.label+'../parameters-AllRadarLines.py').read())
        data = yaml.load(open(self.label+'../parameters_all_radar_lines.yml').read(),
                         Loader=yaml.FullLoader)
        if data != None:
            self.__dict__.update(data)
        filename = self.label+'parameters.yml'
        if os.path.isfile(filename):
            data = yaml.load(open(filename).read(), Loader=yaml.FullLoader)
            if data != None:
                self.__dict__.update(data)

        #Reading the radar dataset
        nbcolumns = 6+self.nbiso+self.is_bedelev+self.is_trace+self.nbhor
        print('nbcolumns:', nbcolumns)
        filename = self.label+'radar-data.txt'
        if os.path.isfile(filename):
            readarray = np.loadtxt(self.label+'radar-data.txt', usecols=range(nbcolumns),
                                   skiprows=1)
        else:
            readarray = np.loadtxt(self.label+'radar-data.dat', usecols=range(nbcolumns), 
                                   skiprows=1)
        if readarray[0, 4] > readarray[-1, 4]:
            readarray = readarray[::-1, :]
        self.LON_raw = readarray[:, 0]
        self.LAT_raw = readarray[:, 1]
        self.x_raw = readarray[:, 2]
        self.y_raw = readarray[:, 3]
        self.distance_raw = readarray[:, 4]
        if self.distance_unit == 'm':
            self.distance_raw = self.distance_raw/1000.
#        self.thk_raw = readarray[:, 5]*self.dilatation_factor
        self.thk_raw = readarray[:, 5]+self.firn_correction
        index = 6
        if self.is_bedelev:
            self.bedelev = readarray[:, index]-self.firn_correction
            index = index+1
        if self.is_trace:
            self.trace = readarray[:, index]
            index = index+1
        self.iso_raw = np.transpose(readarray[:, index:index+self.nbiso])+self.firn_correction
        index = index+self.nbiso
        self.hor_raw = np.transpose(readarray[:, index:index+self.nbhor])+self.firn_correction

        if self.distance_start == 'auto':
            self.distance_start = np.min(self.distance_raw)+self.resolution
        if self.distance_end == 'auto':
            self.distance_end = np.max(self.distance_raw)-self.resolution

        #Interpolation of the datasets
        self.distance = np.arange(self.distance_start, self.distance_end+self.resolution,
                                  self.resolution)
        print(self.distance)
        if self.interp_method == 'stair_aver':
            #TODO: the input function is not a staircase one
            f = interp1d_stair_aver(self.distance_raw, self.thk_raw)
        elif self.interp_method == 'lin_aver':
            #TODO: the input function is not a staircase one
            f = interp1d_lin_aver_withoutnan(self.distance_raw, self.thk_raw)
        else:
            print('interpolation method not recognized')
            quit()
        self.thkreal = f(np.concatenate((self.distance-self.resolution/2,
                                         np.array([self.distance[-1]+self.resolution/2]))))
        self.thk = self.thkreal+0
        self.iso = np.zeros((self.nbiso, np.size(self.distance)))
        self.hor = np.zeros((self.nbhor, np.size(self.distance)))
        self.iso_modage = np.empty_like(self.iso)
        self.iso_modage_sigma = np.empty_like(self.iso)
        self.hor_modage = np.empty_like(self.hor)
        self.iso_EDC = np.zeros(self.nbiso)



        for i in range(self.nbiso):
            if self.interp_method == 'stair_aver':
                f = interp1d_stair_aver(self.distance_raw, self.iso_raw[i, :])
            elif self.interp_method == 'lin_aver':
                f = interp1d_lin_aver(self.distance_raw, self.iso_raw[i, :])
            else:
                print('interpolation method not recognized')
                quit()
            self.iso[i, :] = f(np.concatenate((self.distance-self.resolution/2,
                                               np.array([self.distance[-1]+self.resolution/2]))))


        for i in range(self.nbhor):
            if self.interp_method == 'stair_aver':
                f = interp1d_stair_aver(self.distance_raw, self.hor_raw[i, :])
            elif self.interp_method == 'lin_aver':
                f = interp1d_lin_aver(self.distance_raw, self.hor_raw[i, :])
            else:
                print('interpolation method not recognized')
                quit()
            self.hor[i, :] = f(np.concatenate((self.distance-self.resolution/2,
                                               np.array([self.distance[-1]+self.resolution/2]))))


        f = interp1d(self.distance_raw, self.LON_raw)
        self.LON = f(self.distance)
        f = interp1d(self.distance_raw, self.LAT_raw)
        self.LAT = f(self.distance)

        self.LON_twtt = np.empty_like(self.distance)
        self.LAT_twtt = np.empty_like(self.distance)
        for j in range(np.size(self.distance)):
            self.LON_twtt[j] = self.LON_raw[np.argmin(np.absolute(self.LON_raw-self.LON[j]) +\
                               np.absolute(self.LAT_raw-self.LAT[j]))]
            self.LAT_twtt[j] = self.LAT_raw[np.argmin(np.absolute(self.LON_raw-self.LON[j]) +\
                               np.absolute(self.LAT_raw-self.LAT[j]))]



        #Reading the AICC2012 dataset, calculation of steady age and interpolation
        readarray = np.loadtxt(self.label+'../AICC2012.txt')
        self.AICC2012_depth = readarray[:, 0]
        self.AICC2012_iedepth = readarray[:, 1]
        self.AICC2012_accu = readarray[:, 2]
        self.AICC2012_age = readarray[:, 3]
        self.AICC2012_sigma = readarray[:, 4]

        self.AICC2012_averageaccu = np.sum((self.AICC2012_age[1:] -\
                                    self.AICC2012_age[:-1])*self.AICC2012_accu[:-1])/\
                                    (self.AICC2012_age[-1]-self.AICC2012_age[0])
        print('average accu: ', self.AICC2012_averageaccu)
        self.AICC2012_steadyage = np.cumsum(np.concatenate((np.array([self.AICC2012_age[0]]),\
                                  (self.AICC2012_age[1:]-self.AICC2012_age[:-1])*\
                                  self.AICC2012_accu[:-1]/self.AICC2012_averageaccu)))
        print('steady/unsteady ratio: ', self.AICC2012_steadyage[-1]/self.AICC2012_age[-1])


        if (self.is_EDC and self.calc_isoage):
            for i in range(self.nbiso):
                f = interp1d(self.distance_raw, self.iso_raw[i, :])
                self.iso_EDC[i] = f(self.distance_EDC)

            self.z_err = np.loadtxt(self.label+'z-err.txt')

            f = interp1d(self.AICC2012_depth, self.AICC2012_age)
            self.iso_age = f(self.iso_EDC)
            self.iso_age = np.transpose([self.iso_age])
            self.iso_sigma1 = (f(self.iso_EDC+self.z_err)-f(self.iso_EDC-self.z_err))/2.
            g = interp1d(self.AICC2012_depth, self.AICC2012_sigma)
            self.iso_sigma2 = g(self.iso_EDC)
            self.iso_sigma = np.sqrt(self.iso_sigma1**2+self.iso_sigma2**2)
            self.iso_sigma = np.transpose([self.iso_sigma])

    #Code to be deleted
            self.iso_accu_sigma = np.zeros((self.nbiso, 1))
            self.iso_accu_sigma[0] = self.iso_sigma[0]/(self.iso_age[0]-self.age_surf)
            self.iso_accu_sigma[1:] = np.sqrt(self.iso_sigma[1:]**2+self.iso_sigma[:-1]**2)/\
                                      (self.iso_age[1:]-self.iso_age[:-1])



            output = np.hstack((self.iso_age, self.iso_sigma, self.iso_accu_sigma))
            with open(self.label+'ages.txt', 'w') as f:
                f.write('#age (yr BP)\tsigma_age (yr BP)\tsigma_accu\n')
                np.savetxt(f, output, delimiter="\t")

#Reading ages of isochrones and their sigmas
        if os.path.isfile(self.label+'../ages.txt'):
            readarray = np.loadtxt(self.label+'../ages.txt')
        if os.path.isfile(self.label+'ages.txt'):
            readarray = np.loadtxt(self.label+'ages.txt')
        self.iso_age = np.transpose([readarray[:, 0]])
        self.iso_age = self.iso_age[0:self.nbiso]
        self.iso_sigma = np.transpose([readarray[:, 1]])
        self.iso_sigma = self.iso_sigma[0:self.nbiso]
        f = interp1d(self.AICC2012_age, self.AICC2012_steadyage)
        self.iso_steadyage = f(self.iso_age)




        self.a = self.a*np.ones(np.size(self.distance))
        self.G0 = self.G0*np.ones_like(self.distance)
#        self.mu = self.m/self.a
        self.p = self.p*np.ones(np.size(self.distance))
        self.s = self.s*np.ones(np.size(self.distance))
        self.thkie = np.empty_like(self.distance)

        self.zetagrid = np.arange(0, 1+self.dzeta, self.dzeta)
        self.zetagrid = self.zetagrid[::-1]
        self.zetagrid = np.transpose([self.zetagrid])
        self.zeta = np.ones((np.size(self.zetagrid), np.size(self.distance)))*self.zetagrid
        self.depth = np.empty_like(self.zeta)
        self.depthie = np.empty_like(self.zeta)
        self.zetaie = np.empty_like(self.zeta)
        self.D = np.empty_like(self.zeta[:-1, :])
        self.agesteady = np.zeros((np.size(self.zetagrid), np.size(self.distance)))
        self.age = np.zeros((np.size(self.zetagrid), np.size(self.distance)))
        self.age_density = np.zeros((np.size(self.zetagrid)-1, np.size(self.distance)))
        self.T = np.empty_like(self.age)
        self.T_anal = np.empty_like(self.age)
        self.Tf = np.empty_like(self.distance)
        self.Tm = np.empty_like(self.distance)
        self.alpha = np.empty_like(self.distance)
        self.dist = np.ones((np.size(self.zetagrid), np.size(self.distance)))*self.distance
        self.DeltaT = np.empty_like(self.distance)
        self.G = np.empty_like(self.distance)
        self.m = np.empty_like(self.distance)
        self.mu = np.empty_like(self.distance)
        self.omega_D = np.empty_like(self.age)
        self.omega = np.empty_like(self.age)
        self.tau = np.empty_like(self.age)
        self.uz = np.empty_like(self.age)
        self.sigma_a = np.zeros_like(self.distance)
        self.sigma_m = np.zeros_like(self.distance)
        self.sigma_p = np.zeros_like(self.distance)
        self.sigma_G0 = np.zeros_like(self.distance)
        self.sigma_age = np.zeros_like(self.age)
        self.sigma_logage = np.zeros_like(self.age)
        self.is_fusion = np.empty_like(self.distance)

        self.agebot = np.empty_like(self.distance)
        self.realagebot = np.empty_like(self.distance)
        self.agebot10kyrm = np.empty_like(self.distance)
        self.agebot15kyrm = np.empty_like(self.distance)
        self.age100m = np.empty_like(self.distance)
        self.age150m = np.empty_like(self.distance)
        self.age200m = np.empty_like(self.distance)
        self.age250m = np.empty_like(self.distance)
        self.height0dot6Myr = np.nan*np.ones_like(self.distance)
        self.height0dot8Myr = np.nan*np.ones_like(self.distance)
        self.height1Myr = np.nan*np.ones_like(self.distance)
        self.height1dot2Myr = np.nan*np.ones_like(self.distance)
        self.height1dot5Myr = np.nan*np.ones_like(self.distance)
        self.twtt0dot6Myr = np.nan*np.ones_like(self.distance)
        self.twtt0dot8Myr = np.nan*np.ones_like(self.distance)
        self.twtt1Myr = np.nan*np.ones_like(self.distance)
        self.twtt1dot2Myr = np.nan*np.ones_like(self.distance)
        self.twtt1dot5Myr = np.nan*np.ones_like(self.distance)
        self.sigmabotage = np.empty_like(self.distance)
        self.age_density1Myr = np.nan*np.ones_like(self.distance)
        self.age_density1dot2Myr = np.nan*np.ones_like(self.distance)
        self.age_density1dot5Myr = np.nan*np.ones_like(self.distance)
        self.twttBed = np.nan*np.ones_like(self.distance)
        self.agebotmin = np.empty_like(self.distance)
        self.agebotmax = np.empty_like(self.distance)
        self.amin = np.empty_like(self.distance)
        self.amax = np.empty_like(self.distance)
        self.G0min = np.empty_like(self.distance)
        self.G0max = np.empty_like(self.distance)
        self.mmin = np.empty_like(self.distance)
        self.mmax = np.empty_like(self.distance)
        self.pmin = np.empty_like(self.distance)
        self.pmax = np.empty_like(self.distance)
        self.reso1dot5Myrmin = np.empty_like(self.distance)
        self.reso1dot5Myrmax = np.empty_like(self.distance)
        self.height1dot5Myrmin = np.empty_like(self.distance)
        self.height1dot5Myrmax = np.empty_like(self.distance)



# Model function

    def model1D(self, j):

        #depth grids
        self.thkie[j] = np.interp(self.thk[j], np.concatenate((self.AICC2012_depth,\
                        np.array([self.AICC2012_depth[-1]+3000]))), np.concatenate((\
                        self.AICC2012_iedepth, np.array([self.AICC2012_iedepth[-1]+3000]))))
        self.depth[:, j] = self.thk[j]*(1-self.zeta[:, j])
        self.depthie[:, j] = np.interp(self.depth[:, j], np.concatenate((self.AICC2012_depth,\
                             np.array([self.AICC2012_depth[-1]+3000]))), np.concatenate((\
                             self.AICC2012_iedepth, np.array([self.AICC2012_iedepth[-1]+3000]))))
        self.zetaie[:, j] = (self.thkie[j]-self.depthie[:, j])/self.thkie[j]
        self.D[:, j] = (self.depthie[1:, j]-self.depthie[:-1, j])/(self.depth[1:, j]-\
                       self.depth[:-1, j])


        #Mechanical model
        self.m[j] = 0
        self.G0[j] = 0.05
        self.omega_D[:, j] = 1-(self.p[j]+2)/(self.p[j]+1)*(1-self.zetaie[:, j])+1/\
                             (self.p[j]+1)*(1-self.zetaie[:, j])**(2+self.p[j])
        #Parrenin et al. (CP, 2007a) 2.2 (2)
        self.omega[:, j] = self.s[j]*self.zetaie[:, j]+(1-self.s[j])*self.omega_D[:, j]
        self.tau[:, j] = self.omega[:, j]

        self.age_density[:, j] = np.where((self.tau[1:, j]+self.tau[:-1, j])/2 > 0,\
                                 1/self.a[j]/(self.tau[1:, j]+self.tau[:-1, j])*2, np.nan)
        self.agesteady[:, j] = np.cumsum(np.concatenate((np.array([self.age_surf]),\
                               (self.depthie[1:, j]-self.depthie[:-1, j])*\
                               self.age_density[:, j])), axis=0)


        self.age[:, j] = np.interp(self.agesteady[:, j], np.concatenate((np.array([-1000000000]),\
                         self.AICC2012_steadyage, np.array([1e9*self.AICC2012_steadyage[-1]]))),\
                         np.concatenate((np.array([self.AICC2012_age[0]]), self.AICC2012_age,\
                         np.array([1e9*self.AICC2012_age[-1]]))))

        self.iso_modage[:, j] = np.interp(self.iso[:, j], self.depth[:, j], self.age[:, j])



        return np.concatenate((np.array([self.a[j]]), np.array([self.m[j]]),\
               np.array([self.p[j]]), self.age[:, j], np.log(self.age[1:, j]-\
               self.age_surf), np.array([self.G0[j]])))

    def model1D_finish(self, j):
        
        f = interp1d(self.depth[:, j], self.age[:, j])
        self.agebot[j] = f(max(self.depth[:, j])-60)
        self.realagebot[j] = f(min(self.thk[j], self.thkreal[j]))
        self.age100m[j] = f(max(self.depth[:, j])-100)
        self.age150m[j] = f(max(self.depth[:, j])-150)
        self.age200m[j] = f(max(self.depth[:, j])-200)
        self.age250m[j] = f(max(self.depth[:, j])-250)
        self.hor_modage[:, j] = f(self.hor[:, j])

        self.agebot10kyrm[j] = np.interp(10000., self.age_density[:, j],
                         (self.age[:-1, j]+self.age[1:,j])/2)
        self.agebot15kyrm[j] = np.interp(15000., self.age_density[:, j], 
                         (self.age[:-1, j]+self.age[1:,j])/2)
        if self.agebot10kyrm[j] > self.realagebot[j]:
            self.agebot10kyrm[j] = np.nan
        if self.agebot15kyrm[j] > self.realagebot[j]:
            self.agebot15kyrm[j] = np.nan


        h2 = interp1d(self.age[:, j], self.depth[:, j])
        if self.realagebot[j] >= 1000000.:
            self.age_density1Myr[j] = np.interp(1000000., (self.age[:-1,j]+self.age[1:,j])/2,
                                self.age_density[:,j])
        else:
            self.age_density1Myr[j] = np.nan
        if self.realagebot[j] >= 1200000.:
            self.age_density1dot2Myr[j] = np.interp(1200000., (self.age[:-1,j]+self.age[1:,j])/2,
                                self.age_density[:,j])
        else:
            self.age_density1dot2Myr[j] = np.nan
        if self.realagebot[j] >= 1500000.:
            self.age_density1dot5Myr[j] = np.interp(1500000., (self.age[:-1,j]+self.age[1:,j])/2,
                                self.age_density[:,j])
        else:
            self.age_density1dot5Myr[j] = np.nan
            
        if max(self.age[:, j]) >= 600000:
            self.height0dot6Myr[j] = self.thk[j]-h2(600000)
            self.twtt0dot6Myr[j] = (h2(600000)-self.firn_correction)*100/84.248+250.
        else:
            self.height0dot6Myr[j] = np.nan
            self.twtt0dot6Myr[j] = -98765.0
        if max(self.age[:, j]) >= 800000:
            self.height0dot8Myr[j] = self.thk[j]-h2(800000)
            self.twtt0dot8Myr[j] = (h2(800000)-self.firn_correction)*100/84.248+250.
        else:
            self.height0dot8Myr[j] = np.nan
            self.twtt0dot8Myr[j] = -98765.0
        if max(self.age[:, j]) >= 1000000:
            self.height1Myr[j] = self.thk[j]-h2(1000000)
            self.twtt1Myr[j] = (h2(1000000)-self.firn_correction)*100/84.248+250.
        else:
            self.height1Myr[j] = np.nan
            self.twtt1Myr[j] = -98765.0
        if max(self.age[:, j]) >= 1200000:
            self.height1dot2Myr[j] = self.thk[j]-h2(1200000)
            self.twtt1dot2Myr[j] = (h2(1200000)-self.firn_correction)*100/84.248+250.
        else:
            self.height1dot2Myr[j] = np.nan
            self.twtt1dot2Myr[j] = -98765.0
        if max(self.age[:, j]) >= 1500000:
            self.height1dot5Myr[j] = self.thk[j]-h2(1500000)
            self.twtt1dot5Myr[j] = (h2(1500000)-self.firn_correction)*100/84.248+250.
        else:
            self.height1dot5Myr[j] = np.nan
            self.twtt1dot5Myr[j] = -98765.0
        #TODO: make a function to convert to twtt, and make an array for the different isochrones.
        self.twttBed[j] = (self.thk[j]-self.firn_correction)*100/84.248+250.


    def model(self):  #TODO: kill this or make a call to model(j)
        for j in range(np.size(self.distance)):
            self.model1D(j)

        return np.concatenate((self.a, self.m, self.p, self.age.flatten(), self.G0))

#Residuals function

    def residuals1D(self, variables1D, j):
        var = variables1D+0
        self.a[j] = var[0]
        var = np.delete(var, [0])
#        self.m = variables[np.size(self.distance):2*np.size(self.distance)]
        self.p[j] = var[0]
        if self.p[j] < -0.9:
            self.p[j] = -0.9
        var = np.delete(var, [0])
        if self.invert_thk:
            self.thk[j] = var[0]
            if self.thk[j] < self.iso[-1,j]:
                self.thk[j] = self.iso[-1,j]
            var = np.delete(var, [0])
        if self.invert_s:
            self.s[j] = var[0]
            var = np.delete(var, [0])

        self.model1D(j)
        resi = (self.iso_age.flatten()-self.iso_modage[:, j])/self.iso_sigma.flatten()
        resi = resi[np.where(~np.isnan(resi))]
        resi = np.concatenate((resi, np.array([(self.p[j]-self.p_prior)/\
               self.p_sigma])))
        return resi

    def cost_fct(self, variables1D, j):

        res = self.residuals1D(variables1D, j)
#        cost = 1.-m.exp( -np.sum(res**2)/2. )
        cost = np.sum(res**2)/2.
        return cost



    def jacobian1D(self, j):
        epsilon = np.sqrt(np.diag(self.hess1D))/10000000000.
        model0 = self.model1D(j)
        jacob = np.empty((np.size(self.variables1D), np.size(model0)))
        for i in np.arange(np.size(self.variables1D)):
            self.variables1D[i] = self.variables1D[i]+epsilon[i]
            self.residuals1D(self.variables1D, j)
            model1 = self.model1D(j)
            self.variables1D[i] = self.variables1D[i]-epsilon[i]
            self.residuals1D(self.variables1D, j)
            model2 = self.model1D(j)
            jacob[i] = (model1-model2)/2./epsilon[i]
            self.variables1D[i] = self.variables1D[i]+epsilon[i]
        self.residuals1D(self.variables1D, j)

        return jacob



    def accu_layers(self):
        self.accusteady_layer = np.zeros((self.nbiso, np.size(self.distance)))
        self.accu_layer = np.zeros((self.nbiso, np.size(self.distance)))
        for j in range(np.size(self.distance)):
            f = interp1d(self.depth[:, j], self.age[:, j])
            self.iso_modage[:, j] = f(self.iso[:, j])
            self.accusteady_layer[0, j] = self.a[j]*(self.iso_modage[0, j]-self.age_surf)/\
                                          (self.iso_age[0]-self.age_surf)
            self.accusteady_layer[1:, j] = self.a[j]*(self.iso_modage[1:, j]-\
                                           self.iso_modage[:-1, j])/(self.iso_age[1:]-\
                                           self.iso_age[:-1]).flatten()
        self.accu_layer[0, ] = self.accusteady_layer[0, :]*(self.iso_steadyage[0]-\
                               self.age_surf)/(self.iso_age[0]-self.age_surf)
        self.accu_layer[1:, ] = self.accusteady_layer[1:, ]*(self.iso_steadyage[1:]-\
                                self.iso_steadyage[:-1])/(self.iso_age[1:]-self.iso_age[:-1])

        return

    def sigma1D(self, j):
        jacob = self.jacobian1D(j)


        index = 0
        c_model = np.dot(np.transpose(jacob[:, index:index+1]),
                         np.dot(self.hess1D, jacob[:, index:index+1]))
        self.sigma_a[j] = np.sqrt(np.diag(c_model))[0]
        index = index+1
        c_model = np.dot(np.transpose(jacob[:, index:index+1]),
                         np.dot(self.hess1D, jacob[:, index:index+1]))
        self.sigma_m[j] = np.sqrt(np.diag(c_model))[0]
        index = index+1
        c_model = np.dot(np.transpose(jacob[:, index:index+1]),
                         np.dot(self.hess1D, jacob[:, index:index+1]))
        self.sigma_p[j] = np.sqrt(np.diag(c_model))[0]
        index = index+1
        c_model = np.dot(np.transpose(jacob[:, index:index+np.size(self.age[:, j])]),
                         np.dot(self.hess1D, jacob[:, index:index+np.size(self.age[:, j])]))
        self.sigma_age[:, j] = np.sqrt(np.diag(c_model))
        index = index+np.size(self.age[:, j])
        c_model = np.dot(np.transpose(jacob[:, index:index+np.size(self.age[1:, j])]),
                         np.dot(self.hess1D, jacob[:, index:index+np.size(self.age[1:, j])]))
        self.sigma_logage[1:, j] = np.sqrt(np.diag(c_model))
        self.sigma_logage[0, j] = np.nan
        index = index+np.size(self.age[1:, j])
        c_model = np.dot(np.transpose(jacob[:, index:index+1]),
                         np.dot(self.hess1D, jacob[:, index:index+1]))
        self.sigma_G0[j] = np.sqrt(np.diag(c_model))[0]

        f = interp1d(self.depth[:, j], self.sigma_age[:, j])
        self.sigmabotage[j] = f(self.thk[j]-60.)
        self.iso_modage_sigma[:, j] = f(self.iso[:, j])

        return

#Plotting the raw and interpolated radar datasets
    def data_display(self):
        fig = plt.figure('Data')
        plt.plot(self.distance_raw, self.thk_raw, label='raw bedrock', color='0.5', linewidth=2)
        plt.plot(self.distance, self.thk, label='interpolated bedrock', color='k', linewidth=2)
        for i in range(self.nbiso):
            if i == 0:
                plt.plot(self.distance_raw, self.iso_raw[i, :], color='c', label='raw isochrones')
                plt.plot(self.distance, self.iso[i, :], color='b', label='interpolated isochrones')
            else:
                plt.plot(self.distance_raw, self.iso_raw[i, :], color='c')
                plt.plot(self.distance, self.iso[i, :], color='b')
        for i in range(self.nbhor):
            if i == 0:
                plt.plot(self.distance_raw, self.hor_raw[i, :], color='y', label='raw horizons')
                plt.plot(self.distance, self.hor[i, :], color='g', label='interpolated horizons')
            elif i > 0 and i < self.nbhor-self.nbdsz:
                plt.plot(self.distance_raw, self.hor_raw[i, :], color='y')
                plt.plot(self.distance, self.hor[i, :], color='g')
            elif i == self.nbhor-self.nbdsz:
                plt.plot(self.distance_raw, self.hor_raw[i, :], color='orange', label='raw DSZ')
                plt.plot(self.distance, self.hor[i, :], color='r', label='interpolated DSZ')
            else:
                plt.plot(self.distance_raw, self.hor_raw[i, :], color='orange')
                plt.plot(self.distance, self.hor[i, :], color='r')
        if self.is_EDC:
            EDC_x = np.array([self.distance_EDC, self.distance_EDC])
            EDC_y = np.array([0., 3200.])
            if self.EDC_line_dashed == True:
                plt.plot(EDC_x, EDC_y, label='EDC ice core', color='r', linewidth=2, linestyle='--')
            else:
                plt.plot(EDC_x, EDC_y, label='EDC ice core', color='r', linewidth=2)
        if self.is_NESW:
            plt.xlabel('<NE - distance (km) - SW>')
        else:
            plt.xlabel('distance (km)')
        plt.ylabel('depth (m)')
#        plt.legend(loc=1)
        x1, x2, y1, y2 = plt.axis()
        plt.axis((x1, x2, y2, 0))
        if self.reverse_distance:
            plt.gca().invert_xaxis()
        pp = PdfPages(self.label+'Data.pdf')
        pp.savefig(plt.figure('Data'))
        pp.close()
        plt.close(fig)

#Plot of the model results

    def model_display(self):

#        fig = plt.figure('Model steady')
        fig, plotmodel = plt.subplots()
        plotmodel.set_aspect(self.aspect)
        plt.plot(self.distance, self.thkreal, label='obs. bedrock', color='k', linewidth=2)
        for i in range(self.nbiso):
            if i == 0:
                plt.plot(self.distance, self.iso[i, :], color='w', linewidth=1,
                         label='obs. isochrones')
            else:
                plt.plot(self.distance, self.iso[i, :], color='w', linewidth=1)
        for i in range(self.nbhor):
            if i == 0:
                plt.plot(self.distance, self.hor[i, :], color='0.5', linewidth=1,
                         label='obs. horizons')
            elif i > 0 and i < self.nbhor-self.nbdsz:
                plt.plot(self.distance, self.hor[i, :], color='0.5', linewidth=1)
            elif i == self.nbhor-self.nbdsz:
                plt.plot(self.distance, self.hor[i, :], color='r', linewidth=1,
                         label='obs. DSZ')
            else:
                plt.plot(self.distance, self.hor[i, :], color='r', linewidth=1)
        levels = np.arange(0, 1600, 100)
        levels_color = np.arange(0, 1500, 10)
        plt.contourf(self.dist, self.depth, self.agesteady/1000., levels_color, cmap='jet')
        plt.fill_between(self.distance, self.thkreal, self.thk, color='0.5', label='stagnant ice')
        if self.is_EDC:
            EDC_x = np.array([self.distance_EDC, self.distance_EDC])
            EDC_y = np.array([0., 3200.])
            plt.plot(EDC_x, EDC_y, label='EDC ice core', color='r', linewidth=2)
        if self.is_NESW:
            plt.xlabel('<NE - distance (km) - SW>')
        else:
            plt.xlabel('distance (km)')
        plt.ylabel('depth (m)')
#        plt.legend(loc=2)
        cb = plt.colorbar()
        cb.set_ticks(levels)
        cb.set_ticklabels(levels)
        cb.set_label('Modeled steady age (kyr)')
        x1, x2, y1, y2 = plt.axis()
        if self.max_depth == 'auto':
            self.max_depth = y2
        plt.axis((min(self.distance), max(self.distance), self.max_depth, 0))
        if self.reverse_distance:
            plt.gca().invert_xaxis()
        pp = PdfPages(self.label+'Model-steady.pdf')
        pp.savefig(fig)
        pp.close()
        plt.close(fig)


        fig, plotmodel = plt.subplots()
        plotmodel.set_aspect(self.aspect)
        plt.plot(self.distance, self.thkreal, color='k', linewidth=2, label='bed')
#        plt.legend(loc=1)
        for i in range(self.nbiso):
            if i == 0:
                plt.plot(self.distance, self.iso[i, :], color='w', linewidth=1,
                         label='obs. isochrones')
            else:
                plt.plot(self.distance, self.iso[i, :], color='w', linewidth=1)
        for i in range(self.nbhor):
            if i == 0:
                plt.plot(self.distance, self.hor[i, :], color='0.5', linewidth=1,
                         label='obs. horizons')
            elif i > 0 and i < self.nbhor-self.nbdsz:
                plt.plot(self.distance, self.hor[i, :], color='0.5', linewidth=1)
            elif i == self.nbhor-self.nbdsz:
                plt.plot(self.distance, self.hor[i, :], color='r', linewidth=1,
                         label='obs. DSZ')
            else:
                plt.plot(self.distance, self.hor[i, :], color='r', linewidth=1)
        levels = np.arange(0, 1600, 100)
        levels_color = np.arange(0, 1500, 10)
        plt.contourf(self.dist, self.depth, self.age/1000., levels_color, cmap='jet')
        plt.fill_between(self.distance, self.thkreal, self.thk, color='0.5', label='stagnant ice')
        if self.is_EDC:
            EDC_x = np.array([self.distance_EDC, self.distance_EDC])
            EDC_y = np.array([0., 3200.])
            if self.EDC_line_dashed == True:
                plt.plot(EDC_x, EDC_y, label='EDC ice core', color='r', linewidth=2, linestyle='--')
            else:
                plt.plot(EDC_x, EDC_y, label='EDC ice core', color='r', linewidth=2)
        if self.is_NESW:
            plt.xlabel('<NE - distance (km) - SW>')
        else:
            plt.xlabel('distance (km)')
        plt.ylabel('depth (m)')
        if self.is_legend:
            leg = plt.legend(loc=1)
            frame = leg.get_frame()
            frame.set_facecolor('0.75')
        cb = plt.colorbar()
        cb.set_ticks(levels)
        cb.set_ticklabels(levels)
        cb.set_label('Modeled age (kyr)')
        x1, x2, y1, y2 = plt.axis()
        plt.axis((min(self.distance), max(self.distance), self.max_depth, 0))
        if self.reverse_distance:
            plt.gca().invert_xaxis()
        if self.settick == 'manual':
            plotmodel.set_xticks(np.arange(self.min_tick, self.max_tick+1., self.delta_tick))
        pp = PdfPages(self.label+'Model.pdf')
        pp.savefig(fig)
        pp.close()
        plt.close(fig)


        fig, plotmodel = plt.subplots()
        plotmodel.set_aspect(self.aspect)
        plt.plot(self.distance, self.thkreal, color='k', linewidth=2)
        plt.fill_between(self.distance, self.thkreal, self.thk, color='0.5', label='stagnant ice')
        norm = Normalize(vmin=-5000, vmax=5000)
        for i in range(self.nbiso):
            colorscatter = self.iso_modage[i, :]-self.iso_age[i]
            if i == 0:
                sc = plt.scatter(self.distance, self.iso[i, :], c=colorscatter,
                                 label='obs. isochrones', s=7, edgecolor=None, norm=norm)
            else:
                plt.scatter(self.distance, self.iso[i, :], c=colorscatter, s=7, edgecolor=None,
                            norm=norm)
#        levels = np.arange(0, 1600000, 100000)
#        levels_color = np.arange(0, 1500000, 10000)
#        plt.contourf(self.dist, self.depth, self.age, levels_color)
        if self.is_EDC:
            EDC_x = np.array([self.distance_EDC, self.distance_EDC])
            EDC_y = np.array([0., 3200.])
            if self.EDC_line_dashed == True:
                plt.plot(EDC_x, EDC_y, label='EDC ice core', color='r', linewidth=2,
                         linestyle='--')
            else:
                plt.plot(EDC_x, EDC_y, label='EDC ice core', color='r', linewidth=2)
        if self.is_NESW:
            plt.xlabel('<NE - distance (km) - SW>')
        else:
            plt.xlabel('distance (km)')
        plt.ylabel('depth (m)')
        if self.is_legend:
            print('test')
            leg = plt.legend(loc=1)
            frame = leg.get_frame()
            frame.set_facecolor('0.75')
        cb = plt.colorbar(sc)
#        cb.set_ticks(levels)
#        cb.set_ticklabels(levels)
        cb.set_label('Age misfit (yr)')
        x1, x2, y1, y2 = plt.axis()
        plt.axis((min(self.distance), max(self.distance), self.max_depth, 0))
        if self.reverse_distance:
            plt.gca().invert_xaxis()
        if self.settick == 'manual':
            plotmodel.set_xticks(np.arange(self.min_tick, self.max_tick+1., self.delta_tick))
        pp = PdfPages(self.label+'AgeMisfit.pdf')
        pp.savefig(fig)
        pp.close()
        plt.close(fig)


        fig, plotmodelci = plt.subplots()
        plotmodelci.set_aspect(self.aspect)
        plt.plot(self.distance, self.thkreal, color='k', linewidth=2)
        plt.fill_between(self.distance, self.thkreal, self.thk, color='0.5', label='stagnant ice')
        for i in range(self.nbiso):
            if i == 0:
                plt.plot(self.distance, self.iso[i, :], color='w', label='obs. isochrones')
            else:
                plt.plot(self.distance, self.iso[i, :], color='w')
        for i in range(self.nbhor):
            if i == 0:
                plt.plot(self.distance, self.hor[i, :], color='0.5', label='obs. horizons')
            elif i > 0 and i < self.nbhor-self.nbdsz:
                plt.plot(self.distance, self.hor[i, :], color='0.5')
            elif i == self.nbhor-self.nbdsz:
                plt.plot(self.distance, self.hor[i, :], color='r', label='obs. DSZ')
            else:
                plt.plot(self.distance, self.hor[i, :], color='r')
        levels_log = np.arange(2, 6, 0.1)
        levels = np.power(10, levels_log)
        plt.contourf(self.dist[1:,:], self.depth[1:,:], self.sigma_age[1:,:], levels, norm=LogNorm())
        cb = plt.colorbar()
        cb.set_label('Modeled age confidence interval (yr)')
        levels_labels = np.array([])
        for i in np.arange(2, 6, 1):
            levels_labels = np.concatenate((levels_labels,
                                            np.array([10**i, '', '', '', '', '', '', '', ''])))
        cb.set_ticklabels(levels_labels)
        levels_ticks = np.concatenate((np.arange(100, 1000, 100),
                                       np.arange(1000, 10000, 1000),
                                       np.arange(10000, 100000, 10000),
                                       np.arange(100000, 600000, 100000)))
        cb.set_ticks(levels_ticks)
        if self.is_EDC:
            EDC_x = np.array([self.distance_EDC, self.distance_EDC])
            EDC_y = np.array([0., 3200.])
            if self.EDC_line_dashed == True:
                plt.plot(EDC_x, EDC_y, label='EDC ice core', color='r', linewidth=2, linestyle='--')
            else:
                plt.plot(EDC_x, EDC_y, label='EDC ice core', color='r', linewidth=2)
        if self.is_NESW:
            plt.xlabel('<NE - distance (km) - SW>')
        else:
            plt.xlabel('distance (km)')
        plt.ylabel('depth (m)')
        if self.is_legend:
            leg = plt.legend(loc=1)
            frame = leg.get_frame()
            frame.set_facecolor('0.75')
        x1, x2, y1, y2 = plt.axis()
        plt.axis((min(self.distance), max(self.distance), self.max_depth, 0))
        if self.reverse_distance:
            plt.gca().invert_xaxis()
        if self.settick == 'manual':
            plotmodelci.set_xticks(np.arange(self.min_tick, self.max_tick+1., self.delta_tick))
        pp = PdfPages(self.label+'Model-confidence-interval.pdf')
        pp.savefig(fig)
        pp.close()
        plt.close(fig)


        plt.figure('Thinning')
        plt.plot(self.distance, self.thkreal, label='obs. bedrock', color='k', linewidth=2)
        plt.fill_between(self.distance, self.thkreal, self.thk, color='0.5', label='stagnant ice')
        for i in range(self.nbiso):
            if i == 0:
                plt.plot(self.distance, self.iso[i, :], color='k', label='obs. isochrones')
            else:
                plt.plot(self.distance, self.iso[i, :], color='k')
        for i in range(self.nbhor):
            if i == 0:
                plt.plot(self.distance, self.hor[i, :], color='0.5', label='obs. horizons')
            elif i > 0 and i < self.nbhor-self.nbdsz:
                plt.plot(self.distance, self.hor[i, :], color='0.5')
            elif i == self.nbhor-self.nbdsz:
                plt.plot(self.distance, self.hor[i, :], color='r', label='obs. DSZ')
            else:
                plt.plot(self.distance, self.hor[i, :], color='r')
        plt.contourf(self.dist, self.depth, self.tau)
        if self.is_EDC:
            EDC_x = np.array([self.distance_EDC, self.distance_EDC])
            EDC_y = np.array([0., 3200.])
            plt.plot(EDC_x, EDC_y, label='EDC ice core', color='r', linewidth=2)
        if self.is_NESW:
            plt.xlabel('<NE - distance (km) - SW>')
        else:
            plt.xlabel('distance (km)')
        plt.ylabel('depth (m)')
        plt.legend(loc=2)
        cb = plt.colorbar()
        cb.set_label('Modeled thinning')
        x1, x2, y1, y2 = plt.axis()
        plt.axis((min(self.distance), max(self.distance), self.max_depth, 0))
        if self.reverse_distance:
            plt.gca().invert_xaxis()
        pp = PdfPages(self.label+'Thinning.pdf')
        pp.savefig(plt.figure('Thinning'))
        pp.close()
        plt.close(fig)

        fig = plt.figure('Temperature')
        plt.plot(self.distance, self.thkreal, label='obs. bedrock', color='k', linewidth=2)
        plt.fill_between(self.distance, self.thkreal, self.thk, color='0.5', label='stagnant ice')
        plt.plot(self.distance, np.where(self.is_fusion, np.nan, self.thk), color='b', linewidth=4)
        plt.contourf(self.dist, self.depth, self.T)
        if self.is_EDC:
            EDC_x = np.array([self.distance_EDC, self.distance_EDC])
            EDC_y = np.array([0., 3200.])
            plt.plot(EDC_x, EDC_y, label='EDC ice core', color='r', linewidth=2)
        if self.is_NESW:
            plt.xlabel('<NE - distance (km) - SW>')
        else:
            plt.xlabel('distance (km)')
        plt.ylabel('depth (m)')
        plt.legend(loc=2)
        cb = plt.colorbar()
        cb.set_label('Modeled temperature (K)')
        x1, x2, y1, y2 = plt.axis()
        plt.axis((min(self.distance), max(self.distance), self.max_depth, 0))
        if self.reverse_distance:
            plt.gca().invert_xaxis()
        pp = PdfPages(self.label+'Temperature.pdf')
        pp.savefig(plt.figure('Temperature'))
        pp.close()
        plt.close(fig)



#        fig = plt.figure('Accumulation history')
        lines = [list(zip(self.distance, 917*self.accu_layer[i, :])) for i in range(self.nbiso)]
        z = (self.iso_age.flatten()[1:]+self.iso_age.flatten()[:-1])/2
        z = np.concatenate((np.array([(self.age_surf+self.iso_age.flatten()[0])/2]), z))
        fig, ax = plt.subplots()
        lines = LineCollection(lines, array=z, cmap=plt.cm.rainbow, linewidths=2)
        ax.add_collection(lines)
        ax.autoscale()
        cb = fig.colorbar(lines)
        cb.set_label('Average layer age (yr)')
        if self.is_NESW:
            plt.xlabel('<NE - distance (km) - SW>')
        else:
            plt.xlabel('distance (km)')
        plt.ylabel('steady accumulation (mm-we/yr)')
        if self.reverse_distance:
            plt.gca().invert_xaxis()


        pp = PdfPages(self.label+'AccumulationHistory.pdf')
        pp.savefig(fig)
        pp.close()
        plt.close(fig)



#Plot of the parameters

    def parameters_display(self):
        fig = plt.figure('Parameters')
#        f = plt.figure('Parameters', figsize=(4, 6))

        plotpara = plt.subplot(711,
                               aspect=self.aspect/(self.accu_max-self.accu_min)*self.max_depth/7)
        plt.plot(self.distance, self.a*100, label='accumulation', color='k')
        plt.plot(self.distance, (self.amin)*100, color='k', linestyle='--')
        plt.plot(self.distance, (self.amax)*100, color='k', linestyle='--')
        plt.ylabel('accu. (cm/yr)', fontsize=10)
        plt.tick_params(axis='y', which='both', labelsize=8)
        x1, x2, y1, y2 = plt.axis()
        plt.axis((min(self.distance), max(self.distance), self.accu_min, self.accu_max))
        if self.reverse_distance:
            plt.gca().invert_xaxis()
        if self.is_EDC:
            EDC_x = np.array([self.distance_EDC, self.distance_EDC])
            EDC_y = np.array([self.accu_min, self.accu_max])
            if self.EDC_line_dashed == True:
                plt.plot(EDC_x, EDC_y, label='EDC ice core', color='r', linewidth=2, linestyle='--')
            else:
                plt.plot(EDC_x, EDC_y, label='EDC ice core', color='r', linewidth=2)
#        plt.legend()

        if self.settick == 'manual':
            plotpara.set_xticks(np.arange(self.min_tick, self.max_tick+1., self.delta_tick))
        plotpara = plt.subplot(712, aspect=self.aspect/(self.G0_max-self.G0_min)*self.max_depth/7)
        plt.plot(self.distance, self.G0*1000, label='$G_0$', color='k')
        plt.plot(self.distance, (self.G0min)*1000, color='k', linestyle='--')
        plt.plot(self.distance, (self.G0max)*1000, color='k', linestyle='--')
        plt.ylabel('$G_0$ (mW/m$^2$)', fontsize=10)
        plt.tick_params(axis='y', which='both', labelsize=8)
        x1, x2, y1, y2 = plt.axis()
        plt.axis((min(self.distance), max(self.distance), self.G0_min, self.G0_max))
        if self.reverse_distance:
            plt.gca().invert_xaxis()
        if self.is_EDC:
            EDC_x = np.array([self.distance_EDC, self.distance_EDC])
            EDC_y = np.array([self.G0_min, self.G0_max])
            if self.EDC_line_dashed == True:
                plt.plot(EDC_x, EDC_y, label='EDC ice core', color='r', linewidth=2, linestyle='--')
            else:
                plt.plot(EDC_x, EDC_y, label='EDC ice core', color='r', linewidth=2)
        plotpara.yaxis.tick_right()
        plotpara.yaxis.set_label_position('right')


        if self.settick == 'manual':
            plotpara.set_xticks(np.arange(self.min_tick, self.max_tick+1., self.delta_tick))
        plotpara = plt.subplot(713,
                               aspect=self.aspect/(self.p_max-self.p_min)*self.max_depth/7)
        plt.plot(self.distance, self.p, label='p', color='k')
        plt.plot(self.distance, self.pmin, color='k', linestyle='--')
        plt.plot(self.distance, self.pmax, color='k', linestyle='--')
        plt.ylabel('p parameter', fontsize=10)
        plt.tick_params(axis='y', which='both', labelsize=8)
        if self.is_EDC:
            EDC_x = np.array([self.distance_EDC, self.distance_EDC])
            EDC_y = np.array([self.p_min, self.p_max])
            if self.EDC_line_dashed == True:
                plt.plot(EDC_x, EDC_y, label='EDC ice core', color='r', linewidth=2, linestyle='--')
            else:
                plt.plot(EDC_x, EDC_y, label='EDC ice core', color='r', linewidth=2)
#        plotpara.set_yticks(np.log(np.arange(1., 11.)))
#        plotpara.set_yticks(np.log(np.concatenate((np.arange(1., 10.), 10.*np.arange(1., 10.)))))
#        labels = ["1", "", "", "", "", "", "", "", "", "10"]
#        plotpara.set_yticklabels(labels)
        if self.settick == 'manual':
            plotpara.set_xticks(np.arange(self.min_tick, self.max_tick+1., self.delta_tick))
        x1, x2, y1, y2 = plt.axis()
        plt.axis((min(self.distance), max(self.distance), self.p_min, self.p_max))

        if self.settick == 'manual':
            plotpara.set_xticks(np.arange(self.min_tick, self.max_tick+1., self.delta_tick))
        plotpara = plt.subplot(714,
                               aspect=self.aspect/(self.melting_max-self.melting_min)*\
                               self.max_depth/7)
        plt.plot(self.distance, self.m*1000, label='melting', color='k')
        plt.plot(self.distance, (self.mmin)*1000, color='k', linestyle='--')
        plt.plot(self.distance, (self.mmax)*1000, color='k', linestyle='--')
        plt.ylabel('melting (mm/yr)', fontsize=10)
        plt.tick_params(axis='y', which='both', labelsize=8)
        x1, x2, y1, y2 = plt.axis()
        plt.axis((min(self.distance), max(self.distance), self.melting_min, self.melting_max))
        if self.reverse_distance:
            plt.gca().invert_xaxis()
        if self.is_EDC:
            EDC_x = np.array([self.distance_EDC, self.distance_EDC])
            EDC_y = np.array([self.melting_min, self.melting_max])
            if self.EDC_line_dashed == True:
                plt.plot(EDC_x, EDC_y, label='EDC ice core', color='r', linewidth=2, linestyle='--')
            else:
                plt.plot(EDC_x, EDC_y, label='EDC ice core', color='r', linewidth=2)
        plotpara.yaxis.tick_right()
        plotpara.yaxis.set_label_position('right')

        if self.settick == 'manual':
            plotpara.set_xticks(np.arange(self.min_tick, self.max_tick+1., self.delta_tick))
        plotpara = plt.subplot(715,
                               aspect=self.aspect/(m.log(self.age_max)-m.log(self.age_min))*\
                               self.max_depth/7)
        plt.plot(self.distance, np.log(self.agebot/1000000), label='age 60 m', color='k')
        plt.plot(self.distance, np.log(self.agebotmin/1000000), color='k', linestyle='--')
        plt.plot(self.distance, np.log(self.agebotmax/1000000), color='k', linestyle='--')
        plt.ylabel('age (Myr)', fontsize=10)
        plt.tick_params(axis='y', which='both', labelsize=8)
        plotpara.set_yticks(np.log(np.concatenate((np.arange(1., 10.)/10., np.arange(1., 11.)))))
        labels = ["", "", "", "", "0.5", "", "0.7", "", "", "1", "2", "3", "4", "", "6", "", "",
                  "", "10"]
        plotpara.set_yticklabels(labels)
        x1, x2, y1, y2 = plt.axis()
        plt.axis((min(self.distance), max(self.distance), m.log(self.age_min), m.log(self.age_max)))
        if self.reverse_distance:
            plt.gca().invert_xaxis()
        if self.is_EDC:
            EDC_x = np.array([self.distance_EDC, self.distance_EDC])
            EDC_y = np.array([self.melting_min, self.melting_max])
            if self.EDC_line_dashed == True:
                plt.plot(EDC_x, EDC_y, label='EDC ice core', color='r', linewidth=2, linestyle='--')
            else:
                plt.plot(EDC_x, EDC_y, label='EDC ice core', color='r', linewidth=2)
        plotpara.yaxis.tick_left()
        plotpara.yaxis.set_label_position('left')

        if self.settick == 'manual':
            plotpara.set_xticks(np.arange(self.min_tick, self.max_tick+1., self.delta_tick))
        plotpara = plt.subplot(716, aspect=self.aspect/(self.height_max-self.height_min)*\
                               self.max_depth/7)
        plt.plot(self.distance, self.height1dot5Myr, label='height above bed', color='k')
        plt.plot(self.distance, self.height1dot5Myrmin, label='height above bed', color='k',
                 linestyle='--')
#        plt.plot(self.distance, self.height1dot5Myrmax+self.sigma_height1dot5Myr,
#                 label='height above bed', color='k', linestyle='--')
        plt.ylabel('height (m)', fontsize=10)
        plt.tick_params(axis='y', which='both', labelsize=8)
        x1, x2, y1, y2 = plt.axis()
        plt.axis((min(self.distance), max(self.distance), self.height_min, self.height_max))
        if self.reverse_distance:
            plt.gca().invert_xaxis()
        if self.is_EDC:
            EDC_x = np.array([self.distance_EDC, self.distance_EDC])
            EDC_y = np.array([self.melting_min, self.melting_max])
            if self.EDC_line_dashed == True:
                plt.plot(EDC_x, EDC_y, label='EDC ice core', color='r', linewidth=2, linestyle='--')
            else:
                plt.plot(EDC_x, EDC_y, label='EDC ice core', color='r', linewidth=2)
        plotpara.yaxis.tick_right()
        plotpara.yaxis.set_label_position('right')

        if self.settick == 'manual':
            plotpara.set_xticks(np.arange(self.min_tick, self.max_tick+1., self.delta_tick))
        plotpara = plt.subplot(717,
                               aspect=self.aspect/(self.reso_max-self.reso_min)*self.max_depth/7)
        plt.plot(self.distance, self.age_density1dot5Myr/1000, label='resolution', color='k')
        plt.plot(self.distance, (self.reso1dot5Myrmin)/1000, label='resolution', color='k',
                 linestyle='--')
        plt.plot(self.distance, (self.reso1dot5Myrmax)/1000, label='resolution', color='k',
                 linestyle='--')
        plt.ylabel('resolution (kyr/m)', fontsize=10)
        plt.tick_params(axis='y', which='both', labelsize=8)
        x1, x2, y1, y2 = plt.axis()
        plt.axis((min(self.distance), max(self.distance), self.reso_min, self.reso_max))
        if self.reverse_distance:
            plt.gca().invert_xaxis()
        if self.is_EDC:
            EDC_x = np.array([self.distance_EDC, self.distance_EDC])
            EDC_y = np.array([self.melting_min, self.melting_max])
            if self.EDC_line_dashed == True:
                plt.plot(EDC_x, EDC_y, label='EDC ice core', color='r', linewidth=2, linestyle='--')
            else:
                plt.plot(EDC_x, EDC_y, label='EDC ice core', color='r', linewidth=2)
        plotpara.yaxis.tick_left()
        plotpara.yaxis.set_label_position('left')


        if self.reverse_distance:
            plt.gca().invert_xaxis()
        if self.is_NESW:
            plt.xlabel('<NE - distance (km) - SW>')
        else:
            plt.xlabel('distance (km)')
#        plt.legend()
        fig.subplots_adjust(hspace=0)
        plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
        pp = PdfPages(self.label+'Parameters.pdf')
        pp.savefig(plt.figure('Parameters'))
        pp.close()
        plt.close(fig)

        if self.invert_G0:
            fig = plt.figure('Geothermal heat flux')
            plt.plot(self.distance, self.G0*1000, label='G0', color='k')
            plt.plot(self.distance, (self.G0-self.sigma_G0)*1000, color='k', linestyle='--')
            plt.plot(self.distance, (self.G0+self.sigma_G0)*1000, color='k', linestyle='--')
            plt.ylabel('$G_0$ (mW/m$^2$)')
            if self.reverse_distance:
                plt.gca().invert_xaxis()
#            plt.yaxis.tick_right()
#            plt.yaxis.set_label_position('right')
            pp = PdfPages(self.label+'GeothermalHeatFlux.pdf')
            pp.savefig(fig)
            pp.close()
            plt.close(fig)


    def parameters_save(self):
        output = np.vstack((self.LON, self.LAT, self.distance, self.a, self.sigma_a,
                            self.accu_layer))
        header = '#LON\tLAT\tdistance(km)\taccu(ice-m/yr)\tsigma_accu'
        header = header + '\tlayer ' + str(int(self.age_surf/1000.)) + '-' +\
                 str(int(self.iso_age[0]/1000.)) + 'kyr'
        for i in range(self.nbiso-1):
            header = header + '\tlayer ' + str(int(self.iso_age[i]/1000.)) + '-' +\
                     str(int(self.iso_age[i+1]/1000.)) + 'kyr'
        header = header + '\n'
        with open(self.label+'a.txt', 'w') as f:
            f.write(header)
            np.savetxt(f, np.transpose(output), delimiter="\t")
        output = np.vstack((self.LON, self.LAT, self.distance, self.m, self.sigma_m))
        with open(self.label+'m.txt', 'w') as f:
            f.write('#LON\tLAT\tdistance(km)\tmelting(ice-m/yr)\tsigma_melting\n')
            np.savetxt(f, np.transpose(output), delimiter="\t")
        output = np.vstack((self.LON, self.LAT, self.distance, self.p, self.sigma_p))
        with open(self.label+'p.txt', 'w') as f:
            f.write('#LON\tLAT\tdistance(km)\tp\tsigma_p\n')
            np.savetxt(f, np.transpose(output), delimiter="\t")
        output = np.vstack((self.LON, self.LAT, self.distance, self.G0, self.sigma_G0))
        with open(self.label+'G0.txt', 'w') as f:
            f.write('#LON\tLAT\tdistance(km)\tG0\tsigma_G0\n')
            np.savetxt(f, np.transpose(output), delimiter="\t")

    def bot_age_save(self):
        output = np.vstack((self.LON, self.LAT, self.distance, self.thk, self.agebot,
                            self.agebotmin, self.age100m, self.age150m, self.age200m,
                            self.age250m, self.age_density1Myr, self.age_density1dot2Myr,
                            self.age_density1dot5Myr, self.height0dot6Myr, self.height0dot8Myr,
                            self.height1Myr, self.height1dot2Myr, self.height1dot5Myr,
                            self.agebot10kyrm, self.agebot15kyrm, self.thkreal))

        with open(self.label+'agebottom.txt', 'w') as f:
            f.write('#LON\tLAT\tdistance(km)\tinverted_thickness(m)\tage60m(yr-b1950)'
                    '\tage-min(yr-b1950)'
                    '\tage100m\tage150m\tage200m\tage250\tage_density1Myr\tage_density1.2Myr\t'
                    'age_density1.5Myr\theight0.6Myr\theight0.8Myr\theight1Myr\theight1.2Myr\t'
                    'height1.5Myr'
                    '\tage-10kyrm\tage-15kyrm\treal_thickness'
                    '\n')

            np.savetxt(f, np.transpose(output), delimiter="\t")

    def hor_age_save(self):
        output = np.vstack((self.LON, self.LAT, self.distance, self.hor_modage))
        header = '#LON\tLAT\tdistance(km)'
        for i in range(self.nbhor):
            header = header+'\thor_no_'+str(i+1)
        header = header+'\n'
        with open(self.label+'agehorizons.txt', 'w') as f:
            f.write(header)
            np.savetxt(f, np.transpose(output), delimiter="\t")
        for i in range(self.nbhor):
            print('horizon no:', i+1, ', average age: ', np.nanmean(self.hor_modage[i, :]),
                  ', stdev age: ', np.nanstd(self.hor_modage[i, :]))

    def iso_age_save(self):
        output = np.vstack((self.LON, self.LAT, self.distance, self.iso_modage,
                            self.iso_modage_sigma))
        header = '#LON\tLAT\tdistance(km)'
        for i in range(self.nbiso):
            header = header+'\tiso_no_'+str(i+1)
        for i in range(self.nbiso):
            header = header+'\tiso_no_'+str(i+1)
        header = header+'\n'
        with open(self.label+'ageisochrones.txt', 'w') as f:
            f.write(header)
            np.savetxt(f, np.transpose(output), delimiter="\t")
        for i in range(self.nbiso):
            print('isochrone no:', i+1, ', average age: ', np.nanmean(self.iso_modage[i, :]),
                  ', stdev age: ', np.nanstd(self.iso_modage[i, :]))


    def twtt_save(self):
#        b = np.chararray((np.size(self.distance)), itemsize=20)
#        b[:] = RLlabel
#        output = np.vstack((self.LON_twtt, self.LAT_twtt, self.twtt1Myr))
        current_folder_path, current_folder_name = os.path.split(RL.label[:-1])
        with open(self.label+'twtt.txt', 'w') as f:
            f.write('#LON\tLAT\ttwtt-0.6Myr(lk)\ttwtt-0.8Myr(lk)\ttwtt-1Myr(lk)\ttwtt-1.2Myr(lk)'
                    '\ttwtt-1.5Myr(lk)\tself.twttBed(lk)\tLabel\n')
            for j in range(np.size(self.distance)):
                f.write('{0:11}'.format(str(self.LON_twtt[j]))+'\t'+\
                        '{0:12}'.format(str(self.LAT_twtt[j]))+'\t'+\
                        '{0:13}'.format(str(self.twtt0dot6Myr[j]))+'\t'+\
                        '{0:13}'.format(str(self.twtt0dot8Myr[j]))+'\t'+\
                        '{0:13}'.format(str(self.twtt1Myr[j]))+'\t'+\
                        '{0:13}'.format(str(self.twtt1dot2Myr[j]))+'\t'+\
                        '{0:13}'.format(str(self.twtt1dot5Myr[j]))+'\t'+\
                        '{0:13}'.format(str(self.twttBed[j]))+'\t'+current_folder_name+'\n')
#            np.savetxt(f, np.transpose(output), delimiter="\t")

    def layer_depth_save(self):
        output = np.vstack((self.LON, self.LAT, self.distance))
        output = np.vstack((output, self.iso[0, :]/2))
        output = np.vstack((output, (self.iso[:-1, :]+self.iso[1:, :])/2))
        header = '#LON\tLAT\tdistance(km)'
        for i in range(self.nbiso):
            header = header+'\tlayer_no_'+str(i+1)
        header = header+'\n'
        with open(self.label+'depthlayers.txt', 'w') as f:
            f.write(header)
            np.savetxt(f, np.transpose(output), delimiter="\t")

    def age_2D_save(self):
        # save model results as arrays 
        # 2D data
        np.savetxt(self.label + 'distance.txt',self.dist)
        np.savetxt(self.label + 'depth.txt',self.depth)
        np.savetxt(self.label + 'age.txt',self.age/1000.)
        np.savetxt(self.label + 'ages_density.txt',self.age_density)
        np.savetxt(self.label + 'thinning.txt',self.tau)
        np.savetxt(self.label + 'sigma_age.txt',self.sigma_age)
        
        # along line data
        data_line = pd.DataFrame({'Distance': self.distance[:], 
                                  'thick':self.thkreal[:],
                                  'stagnant': self.thkreal[:] - self.thk[:]})
        data_line.to_csv (self.label + 'line.csv', index = False)
        
        # on the isochrones
        np.savetxt(self.label + 'isomod.txt',self.iso_modage)       
        np.savetxt(self.label + 'iso_obs_interp.txt',self.iso)

    def EDC(self):
        f = interp1d(self.distance, self.age)
        age_EDC = f(self.distance_EDC)
        g = interp1d(self.distance, self.sigma_age)
        sigmaage_EDC = g(self.distance_EDC)
        h = interp1d(self.distance, self.depth)
        depth_EDC = h(self.distance_EDC)
        print('Thickness at EDC is:', depth_EDC[-1])
        i = interp1d(depth_EDC, age_EDC)
        age_EDC_bot = i(max(depth_EDC)-60)
        j = interp1d(depth_EDC, sigmaage_EDC)
        sigmaage_EDC_bot = j(max(depth_EDC)-60)
        print('Age at EDC at 3200 m depth is: ', age_EDC_bot, '+-', sigmaage_EDC_bot)
        f = interp1d(self.distance, self.p)
        p_EDC = f(self.distance_EDC)
        print('p parameter at EDC is: ', p_EDC)
        f = interp1d(self.distance, self.a)
        a_EDC = f(self.distance_EDC)
        print('accumulation at EDC is: ', a_EDC)
        f = interp1d(self.distance, self.m)
        m_EDC = f(self.distance_EDC)
        print('melting at EDC is: ', m_EDC)

        readarray = np.loadtxt(self.label+'temperatures_EDC.txt')
        datatempEDC_depth = readarray[:, 0]
        datatempEDC_temp = readarray[:, 1]+273.15
        f = interp1d(self.distance, self.T)
        temp_EDC = f(self.distance_EDC)
        f = interp1d(self.distance, self.T_anal)
        temp_anal_EDC = f(self.distance_EDC)
        plt.figure('Temperature at EDC')
        plt.plot(temp_EDC, depth_EDC, label='model')
        plt.plot(temp_anal_EDC, depth_EDC, label='analytical solution')
        plt.plot(datatempEDC_temp, datatempEDC_depth, label='data')
        plt.ylabel('depth (m)')
        plt.xlabel('temperature (K)')
        x1, x2, y1, y2 = plt.axis()
        plt.axis((x1, x2, y2, 0.))
        plt.legend()
        pp = PdfPages(self.label+'temperatures_EDC.pdf')
        pp.savefig(plt.figure('Temperature at EDC'))
        pp.close()

        output = np.vstack((depth_EDC, age_EDC, temp_EDC))
        header = 'depth\tage\ttemperature\n'
        with open(self.label+'EDCresults.txt', 'w') as f:
            f.write(header)
            np.savetxt(f, np.transpose(output), delimiter="\t")

#Main
RLlabel = sys.argv[1]
if RLlabel[-1] != '/':
    RLlabel = RLlabel+'/'
print('Radar line is: ', RLlabel)
print('Creation of the radar line')
RL = RadarLine(RLlabel)
print('Initialization of radar line')
RL.init()
print('Data display')
RL.data_display()
if RL.opt_method == 'leastsq':
    print('Optimization by leastsq')
    RL.variables = np.concatenate((RL.a, RL.p))
    if RL.invert_G0:
        RL.variables = np.concatenate((RL.variables, RL.G0))
    if RL.invert_thk:
        RL.variables = np.concatenate((RL.variables, RL.thk))
#        self.variables = np.concatenate((self.a, self.m, self.s))
    RL.variables, RL.hess, infodict, mesg, ier = leastsq(RL.residuals, RL.variables, full_output=1)
    print(mesg)
    RL.residuals(RL.variables)
    print('Calculation of confidence intervals')
    if not RL.calc_sigma:
        RL.hess = np.zeros((np.size(RL.variables), np.size(RL.variables)))
    RL.sigma()
elif RL.opt_method == 'none':
    RL.variables = np.concatenate((RL.a, RL.p))
    if RL.invert_G0:
        RL.variables = np.concatenate((RL.variables, RL.G0))
    if RL.invert_thk:
        RL.variables = np.concatenate((RL.variables, RL.thk))
    print('No optimization')
    RL.residuals(RL.variables)
elif RL.opt_method == 'none1D':
    print('Forward model 1D')
    for j in range(np.size(RL.distance)):
        RL.variables1D = np.array([RL.a[j], RL.p[j]])
        if RL.invert_G0:
            RL.variables1D = np.append(RL.variables1D, RL.G0[j])
        if RL.invert_thk:
            RL.variables1D = np.append(RL.variables1D, RL.thk[j])
        RL.residuals1D(RL.variables1D, j)
elif RL.opt_method == 'leastsq1D':
    print('Optimization by leastsq1D')
    for j in range(np.size(RL.distance)):
        print('index along the radar line: ', j)
        RL.variables1D = np.array([RL.a[j], RL.p[j]])
        if RL.invert_thk:
            RL.variables1D = np.append(RL.variables1D, RL.thk[j])
        if RL.invert_s:
            RL.variables1D = np.append(RL.variables1D, RL.s[j])

        RL.variables1D, RL.hess1D, infodict, mesg, ier = leastsq(RL.residuals1D, RL.variables1D,
                                                                 args=(j), full_output=1)
        RL.residuals1D(RL.variables1D, j)
        RL.model1D_finish(j)
        print(RL.variables1D)
        if not RL.calc_sigma:
            RL.hess1D = np.zeros((np.size(RL.variables1D), np.size(RL.variables1D)))
        if np.size(RL.hess1D) != 1:
            RL.sigma1D(j)
    RL.agebotmin = RL.agebot-RL.sigmabotage
    RL.agebotmax = RL.agebot+RL.sigmabotage


else:
    print(RL.opt_method, ': Optimization method not recognized.')
    quit()
print('calculating per layer accumulations')

RL.accu_layers()





print('Model display')
RL.model_display()
print('parameters display')
if RL.is_EDC:
    RL.EDC()
RL.parameters_display()
RL.parameters_save()
#RL.max_age()
RL.bot_age_save()
RL.hor_age_save()
RL.iso_age_save()
RL.twtt_save()
RL.layer_depth_save()
RL.age_2D_save()
print('Program execution time: ', time.time() - START_TIME, 'seconds')
