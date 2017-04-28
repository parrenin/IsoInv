import time
import sys
import math as m
import numpy as np
import matplotlib.pyplot as plt
import os
import random
from matplotlib.collections import LineCollection
from matplotlib.colors import LogNorm, Normalize
from scipy.interpolate import interp1d
from scipy.optimize import leastsq, basinhopping, minimize
from matplotlib.backends.backend_pdf import PdfPages
from scipy.special import erf


###Registration of start time
start_time = time.time()     



#Physical constants
Kg0=9.828
Kg1=-5.7*10**-3
Lf=333.5*10**3
rhog=917.
cg0=152.5
cg1=7.122
ggrav=9.81
#Tf0=273.16          #This is the value from Cuffey and Paterson (2010)
#Tf1=-9.8e-8          #This is the value from Cuffey and Paterson (2010)
Tf0=273.16-0.024         #This is the value from Catherine Ritz's thesis
Tf1=-7.4e-8         #This is the value from Catherine Ritz's thesis

def cg(T):
#    return 2097.*np.ones(np.shape(T))    #at 0 degC
    return cg0+cg1*T

def Kg(T, D):
    """From Patterson eqs. 9.2 and 9.4"""
#    return 2.10*np.ones(np.shape(T))     #at 0 degC
#    return Kg0*np.exp(Kg1*T)
    KiT=Kg0*np.exp(Kg1*T)
    return (2.*KiT*D)/(3.-D)

def Tf(P):
    return Tf0+Tf1*P

def interp1d_stair_aver(x, y):   #TODO: deal with the case x not sorted
    """
    Interpolation of a staircase function using averaging.
    This function returns nan outside of the input abscissa range.
    """
    def f(xp):
        yp=np.empty(np.size(xp)-1)
        xmod=x[~(np.isnan(x)+np.isnan(y))]
        ymod=y[~(np.isnan(x)+np.isnan(y))]
        yint=np.cumsum(np.concatenate((np.array([0]),ymod[:-1]*(xmod[1:]-xmod[:-1]))))
        g=interp1d(xmod,yint, bounds_error=False, fill_value=np.nan)
#        yp=np.where((xp[:-1]>min(xmod))*(xp[1:]<max(xmod)),(g(xp[1:])-g(xp[:-1]))/(xp[1:]-xp[:-1]),np.nan)     #Maybe this is suboptimal since we compute twice g(xp[i])
        yp=np.where((xp[:-1]>min(xmod))*(xp[1:]<max(xmod)),(g(xp[1:])-g(xp[:-1]))/(xp[1:]-xp[:-1]),np.nan)     #Maybe this is suboptimal since we compute twice g(xp[i])
        return yp

    return f

def interp1d_stair_aver_withnan(x, y):   #TODO: deal with the case x not sorted
    """
    Interpolation of a staircase function using averaging.
    This function returns nan when there are all nans in one interpolation interval.
    """
    def f(xp):
        xmod=x[~(np.isnan(x)+np.isnan(y))]
        ymod=y[~(np.isnan(x)+np.isnan(y))]
        yp=np.empty(np.size(xp)-1)
        yint=np.cumsum(np.concatenate((np.array([0]),ymod[:-1]*(xmod[1:]-xmod[:-1]))))
        g=interp1d(xmod,yint, bounds_error=False, fill_value=np.nan)
#        yp=np.where((xp[:-1]>min(xmod))*(xp[1:]<max(xmod)),(g(xp[1:])-g(xp[:-1]))/(xp[1:]-xp[:-1]),np.nan)     #Maybe this is suboptimal since we compute twice g(xp[i])
        yp=np.where((xp[:-1]>min(xmod))*(xp[1:]<max(xmod)),(g(xp[1:])-g(xp[:-1]))/(xp[1:]-xp[:-1]),np.nan)     #Maybe this is suboptimal since we compute twice g(xp[i])
        for i in range(np.size(xp)-1):
            if np.isnan(y[np.where((x>=xp[i])*(x<xp[i+1]))]).all():
                yp[i]=np.nan
        return yp

    return f

      
def interp1d_lin_aver_withoutnan(x,y):   #FIXME: there is a problem in this routine when the x are in decreasing order.
    """
    Interpolation of a linear by parts function using averaging.
    This function returns nan when there are all nans in one interpolation interval.
    """
    def f(xp):
        yp=np.empty(np.size(xp)-1)
        for i in range(np.size(xp)-1):
#            print i, xp[i], xp[i+1]
            xmod=x[~(np.isnan(x)+np.isnan(y))]
            ymod=y[~(np.isnan(x)+np.isnan(y))]
            xmod2=xmod[np.where((xmod>xp[i])*(xmod<xp[i+1]))]
            ymod2=ymod[np.where((xmod>xp[i])*(xmod<xp[i+1]))]
            xmod3=np.concatenate((np.array([xp[i]]),xmod2,np.array([xp[i+1]])))
            g=interp1d(xmod,ymod, bounds_error=False, fill_value=np.nan)
            ymod3=np.concatenate((np.array([g(xp[i])]),ymod2,np.array([g(xp[i+1])])))
#                print xmod3
#                print ymod3
            if np.isnan(ymod3).all():
                yp[i]=np.nan
            else:
                xmod4=xmod3[np.where(~(np.isnan(ymod3)+np.isnan(xmod3)))]
                ymod4=ymod3[np.where(~(np.isnan(ymod3)+np.isnan(xmod3)))]
#                if i==9:
#                    print xmod4,ymod4
                yp[i]=np.sum((ymod4[1:]+ymod4[:-1])/2*(xmod4[1:]-xmod4[:-1]))
                yp[i]=yp[i]/(xmod4[-1]-xmod4[0])
#                print yp[i]
        return yp
    return f


def interp1d_lin_aver(x,y):   #FIXME: there is a problem in this routine when the x are in decreasing order.
    """
    Interpolation of a linear by parts function using averaging.
    This function returns nan when there are all nans in one interpolation interval.
    """
    def f(xp):
        yp=np.empty(np.size(xp)-1)
        for i in range(np.size(xp)-1):
#            print i, xp[i], xp[i+1]
            xmod=x[~(np.isnan(x)+np.isnan(y))]
            ymod=y[~(np.isnan(x)+np.isnan(y))]
            xmod2=xmod[np.where((xmod>xp[i])*(xmod<xp[i+1]))]
            ymod2=ymod[np.where((xmod>xp[i])*(xmod<xp[i+1]))]
            xmod3=np.concatenate((np.array([xp[i]]),xmod2,np.array([xp[i+1]])))
            g=interp1d(x,y, bounds_error=False, fill_value=np.nan)
            ymod3=np.concatenate((np.array([g(xp[i])]),ymod2,np.array([g(xp[i+1])])))
#                print xmod3
#                print ymod3
            if np.isnan(ymod3).all():
                yp[i]=np.nan
            else:
                xmod4=xmod3[np.where(~(np.isnan(ymod3)+np.isnan(xmod3)))]
                ymod4=ymod3[np.where(~(np.isnan(ymod3)+np.isnan(xmod3)))]
#                if i==9:
#                    print xmod4,ymod4
                yp[i]=np.sum((ymod4[1:]+ymod4[:-1])/2*(xmod4[1:]-xmod4[:-1]))
                yp[i]=yp[i]/(xmod4[-1]-xmod4[0])
#                print yp[i]
        return yp
    return f


class RadarLine:

    def __init__(self, label):
        self.label=label

    def init(self):
        self.is_bedelev=False
        self.is_trace=False
        self.nbdsz=0
        self.calc_sigma=True
        self.invert_G0=False
        self.settick='auto'
        self.interp_method='lin_aver'
        self.distance_unit='km'
        self.nbhor=0

        #definition of some global parameters
        execfile(self.label+'../parameters-AllRadarLines.py')
        filename=self.label+'parameters.py'
        if os.path.isfile(filename):
            execfile(filename)


        #Reading the radar dataset
        nbcolumns=6+self.nbiso+self.is_bedelev+self.is_trace+self.nbhor
        print 'nbcolumns:',nbcolumns
        readarray=np.loadtxt(self.label+'radar-data.txt', usecols=range(nbcolumns))
        if readarray[0,4]>readarray[-1,4]:
            readarray=readarray[::-1,:]
        self.LON_raw=readarray[:,0]
        self.LAT_raw=readarray[:,1]
        self.x_raw=readarray[:,2]
        self.y_raw=readarray[:,3]
        self.distance_raw=readarray[:,4]
        if self.distance_unit=='m':
            self.distance_raw=self.distance_raw/1000.
#        self.thk_raw=readarray[:,5]*self.dilatation_factor
        self.thk_raw=readarray[:,5]+self.firn_correction
        index=6
        if self.is_bedelev:
            self.bedelev=readarray[:,index]-self.firn_correction
            index=index+1
        if self.is_trace:
            self.trace=readarray[:,index]
            index=index+1
        self.iso_raw=np.transpose(readarray[:,index:index+self.nbiso])+self.firn_correction
        index=index+self.nbiso
        self.hor_raw=np.transpose(readarray[:,index:index+self.nbhor])+self.firn_correction

#        #Interpolation of the radar dataset
#        if self.reverse_distance:
#            toto=self.distance_raw[-1]-self.distance_end
#            self.distance_end=self.distance_raw[-1]-self.distance_start
#            self.distance_start=toto+0.
##            print self.distance_start, self.distance_end
#            self.distance_raw=self.distance_raw[-1]-self.distance_raw
#            self.LON_raw=self.LON_raw[::-1]
#            self.LAT_raw=self.LAT_raw[::-1]
#            self.x_raw=self.x_raw[::-1]
#            self.y_raw=self.y_raw[::-1]
#            self.distance_raw=self.distance_raw[::-1]
#            self.thk_raw=self.thk_raw[::-1]
#            self.iso_raw=self.iso_raw[:,::-1]
       
#        if self.reset_distance:
#            self.distance_raw=self.distance_raw-self.distance_start


        if self.distance_start=='auto':
            self.distance_start=int(np.min(self.distance_raw)+2.99)
        if self.distance_end=='auto':
            self.distance_end=int(np.max(self.distance_raw)-2.)

        #Interpolation of the datasets
        self.distance=np.arange(self.distance_start, self.distance_end+self.resolution, self.resolution) 
#        self.distance=np.arange(90, 100+self.resolution, self.resolution) 
#        f=interp1d(self.distance_raw,self.thk_raw)  #TODO: we want to integrate here to smooth the record. Same for iso.
#        self.thk=f(self.distance)
        if self.interp_method=='stair_aver':
            f=interp1d_stair_aver(self.distance_raw,self.thk_raw)    #TODO: the input function is not a staircase one
        elif self.interp_method=='lin_aver':
            f=interp1d_lin_aver_withoutnan(self.distance_raw,self.thk_raw)    #TODO: the input function is not a staircase one
        else:
            print 'interpolation method not recognized'
            quit()
        self.thkreal=f(np.concatenate((self.distance-self.resolution/2, np.array([self.distance[-1]+self.resolution/2]))))
        self.thk=self.thkreal+0
#        print 'Interpolation of bedrock: done.'
        self.iso=np.zeros((self.nbiso,np.size(self.distance)))
        self.hor=np.zeros((self.nbhor,np.size(self.distance)))
        self.iso_modage=np.empty_like(self.iso)
        self.iso_modage_sigma=np.empty_like(self.iso)
        self.hor_modage=np.empty_like(self.hor)
        self.iso_EDC=np.zeros(self.nbiso)



        for i in range(self.nbiso):
            if self.interp_method=='stair_aver':
                f=interp1d_stair_aver(self.distance_raw,self.iso_raw[i,:])
            elif self.interp_method=='lin_aver':
                f=interp1d_lin_aver(self.distance_raw,self.iso_raw[i,:])
            else:
                print 'interpolation method not recognized'
                quit()
            self.iso[i,:]=f(np.concatenate((self.distance-self.resolution/2, np.array([self.distance[-1]+self.resolution/2]))))
#            print 'Interpolation of isochrone no: ',i,': done.'

        for i in range(self.nbhor):
            if self.interp_method=='stair_aver':
                f=interp1d_stair_aver(self.distance_raw,self.hor_raw[i,:])
            elif self.interp_method=='lin_aver':
                f=interp1d_lin_aver(self.distance_raw,self.hor_raw[i,:])
            else:
                print 'interpolation method not recognized'
                quit()
            self.hor[i,:]=f(np.concatenate((self.distance-self.resolution/2, np.array([self.distance[-1]+self.resolution/2]))))


        f=interp1d(self.distance_raw,self.LON_raw)
        self.LON=f(self.distance)
        f=interp1d(self.distance_raw, self.LAT_raw)
        self.LAT=f(self.distance)

        self.LON_twtt=np.empty_like(self.distance)
        self.LAT_twtt=np.empty_like(self.distance)
        for j in range(np.size(self.distance)):
            self.LON_twtt[j]=self.LON_raw[np.argmin(np.absolute(self.LON_raw-self.LON[j])+np.absolute(self.LAT_raw-self.LAT[j]))]
            self.LAT_twtt[j]=self.LAT_raw[np.argmin(np.absolute(self.LON_raw-self.LON[j])+np.absolute(self.LAT_raw-self.LAT[j]))]



        #Reading the AICC2012 dataset, calculation of steady age and interpolation
        readarray=np.loadtxt(self.label+'../AICC2012.txt')
        self.AICC2012_depth=readarray[:,0]
        self.AICC2012_iedepth=readarray[:,1]
        self.AICC2012_accu=readarray[:,2]
        self.AICC2012_age=readarray[:,3]
        self.AICC2012_sigma=readarray[:,4]

#        self.AICC2012_steadyage=np.cumsum(np.concatenate((np.array([self.AICC2012_age[0]]),(self.AICC2012_age[1:]-self.AICC2012_age[:-1])*self.AICC2012_accu[:-1]/(self.AICC2012_iedepth[1:]-self.AICC2012_iedepth[:-1])/self.AICC2012_accu[0]*(self.AICC2012_depth[1:]-self.AICC2012_depth[:-1]))))
        self.AICC2012_averageaccu=np.sum((self.AICC2012_age[1:]-self.AICC2012_age[:-1])*self.AICC2012_accu[:-1])/(self.AICC2012_age[-1]-self.AICC2012_age[0])
        print 'average accu: ',self.AICC2012_averageaccu
        self.AICC2012_steadyage=np.cumsum(np.concatenate((np.array([self.AICC2012_age[0]]),(self.AICC2012_age[1:]-self.AICC2012_age[:-1])*self.AICC2012_accu[:-1]/self.AICC2012_averageaccu)))
#        self.AICC2012_steadyage=self.AICC2012_steadyage*self.AICC2012_age[-1]/self.AICC2012_steadyage[-1]
        print 'steady/unsteady ratio: ', self.AICC2012_steadyage[-1]/self.AICC2012_age[-1]

#        print 'prior value on melting',self.m_EDC

        if (self.is_EDC and self.calc_isoage):
#        if self.is_EDC:
            for i in range(self.nbiso):
                f=interp1d(self.distance_raw,self.iso_raw[i,:])
                self.iso_EDC[i]=f(self.distance_EDC)

            self.z_err=np.loadtxt(self.label+'z-err.txt')

            f=interp1d(self.AICC2012_depth,self.AICC2012_age)
            self.iso_age=f(self.iso_EDC)
            self.iso_age=np.transpose([self.iso_age])
            self.iso_sigma1=(f(self.iso_EDC+self.z_err)-f(self.iso_EDC-self.z_err))/2.
            g=interp1d(self.AICC2012_depth,self.AICC2012_sigma)
            self.iso_sigma2=g(self.iso_EDC)
            self.iso_sigma=np.sqrt(self.iso_sigma1**2+self.iso_sigma2**2)
            self.iso_sigma=np.transpose([self.iso_sigma])

    #Code to be deleted
            self.iso_accu_sigma=np.zeros((self.nbiso,1))
            self.iso_accu_sigma[0]=self.iso_sigma[0]/(self.iso_age[0]-self.age_surf)
            self.iso_accu_sigma[1:]=np.sqrt(self.iso_sigma[1:]**2+self.iso_sigma[:-1]**2)/(self.iso_age[1:]-self.iso_age[:-1])


            

            output=np.hstack((self.iso_age, self.iso_sigma, self.iso_accu_sigma))
#            with open(self.label+'-ages.txt','w') as f:
#                f.write('#age (yr BP)\tsigma_age (yr BP)\n')
#                np.savetxt(f,output, delimiter="\t") 
            with open(self.label+'ages.txt','w') as f:
                f.write('#age (yr BP)\tsigma_age (yr BP)\tsigma_accu\n')
                np.savetxt(f,output, delimiter="\t")

#Reading ages of isochrones and their sigmas
        if os.path.isfile(self.label+'../ages.txt'):
            readarray=np.loadtxt(self.label+'../ages.txt')
        if os.path.isfile(self.label+'ages.txt'):
            readarray=np.loadtxt(self.label+'ages.txt')
        self.iso_age=np.transpose([readarray[:,0]])
        self.iso_age=self.iso_age[0:self.nbiso]
        self.iso_sigma=np.transpose([readarray[:,1]])
        self.iso_sigma=self.iso_sigma[0:self.nbiso]
        f=interp1d(self.AICC2012_age,self.AICC2012_steadyage)
        self.iso_steadyage=f(self.iso_age)




        self.a=self.a*np.ones(np.size(self.distance))
        self.G0=self.G0*np.ones_like(self.distance)
#        self.mu=self.m/self.a
        self.pprime=self.pprime*np.ones(np.size(self.distance))
        self.p=np.empty_like(self.pprime)
        self.s=self.s*np.ones(np.size(self.distance))
        self.thkie=np.empty_like(self.distance)

        self.zetagrid=np.arange(0,1+self.dzeta,self.dzeta)
        self.zetagrid=self.zetagrid[::-1]
        self.zetagrid=np.transpose([self.zetagrid])
        self.zeta=np.ones((np.size(self.zetagrid),np.size(self.distance)))*self.zetagrid
        self.depth=np.empty_like(self.zeta)
        self.depthie=np.empty_like(self.zeta)
        self.zetaie=np.empty_like(self.zeta)
        self.D=np.empty_like(self.zeta[:-1,:])
        self.agesteady=np.zeros((np.size(self.zetagrid),np.size(self.distance)))
        self.age=np.zeros((np.size(self.zetagrid),np.size(self.distance)))
        self.age_density=np.zeros((np.size(self.zetagrid)-1,np.size(self.distance)))
        self.T=np.empty_like(self.age)
        self.T_anal=np.empty_like(self.age)
        self.Tf=np.empty_like(self.distance)
        self.Tm=np.empty_like(self.distance)
        self.alpha=np.empty_like(self.distance)
        self.dist=np.ones((np.size(self.zetagrid),np.size(self.distance)))*self.distance
        self.DeltaT=np.empty_like(self.distance)
        self.G=np.empty_like(self.distance)
        self.m=np.empty_like(self.distance)
        self.mu=np.empty_like(self.distance)
        self.omega_D=np.empty_like(self.age)
        self.omega=np.empty_like(self.age)
        self.tau=np.empty_like(self.age)
        self.uz=np.empty_like(self.age)
        self.sigma_a=np.zeros_like(self.distance)
        self.sigma_m=np.zeros_like(self.distance)
        self.sigma_pprime=np.zeros_like(self.distance)
        self.sigma_G0=np.zeros_like(self.distance)
        self.sigma_age=np.zeros_like(self.age)
        self.sigma_logage=np.zeros_like(self.age)
        self.is_fusion=np.empty_like(self.distance)

        self.agebot=np.empty_like(self.distance)
        self.age100m=np.empty_like(self.distance)
        self.age150m=np.empty_like(self.distance)
        self.age200m=np.empty_like(self.distance)
        self.age250m=np.empty_like(self.distance)
        self.height0dot6Myr=np.nan*np.ones_like(self.distance)
        self.height0dot8Myr=np.nan*np.ones_like(self.distance)
        self.height1Myr=np.nan*np.ones_like(self.distance)
        self.height1dot2Myr=np.nan*np.ones_like(self.distance)
        self.height1dot5Myr=np.nan*np.ones_like(self.distance)
        self.twtt0dot6Myr=np.nan*np.ones_like(self.distance)
        self.twtt0dot8Myr=np.nan*np.ones_like(self.distance)
        self.twtt1Myr=np.nan*np.ones_like(self.distance)
        self.twtt1dot2Myr=np.nan*np.ones_like(self.distance)
        self.twtt1dot5Myr=np.nan*np.ones_like(self.distance)
        self.sigmabotage=np.empty_like(self.distance)
        self.agebotmin=np.empty_like(self.distance)
        self.age_density1Myr=np.nan*np.ones_like(self.distance)
        self.age_density1dot2Myr=np.nan*np.ones_like(self.distance)
        self.age_density1dot5Myr=np.nan*np.ones_like(self.distance)
        self.twttBed=np.nan*np.ones_like(self.distance)



# Model function

    def model1D(self, j):

        #depth grids
        self.thkie[j]=np.interp(self.thk[j],np.concatenate((self.AICC2012_depth,np.array([self.AICC2012_depth[-1]+3000]))),np.concatenate((self.AICC2012_iedepth,np.array([self.AICC2012_iedepth[-1]+3000]))))
#        print np.shape(self.zeta[:,j])
        self.depth[:,j]=self.thk[j]*(1-self.zeta[:,j])
        self.depthie[:,j]=np.interp(self.depth[:,j],np.concatenate((self.AICC2012_depth,np.array([self.AICC2012_depth[-1]+3000]))),np.concatenate((self.AICC2012_iedepth,np.array([self.AICC2012_iedepth[-1]+3000]))))
        self.zetaie[:,j]=(self.thkie[j]-self.depthie[:,j])/self.thkie[j]
        self.D[:,j]=(self.depthie[1:,j]-self.depthie[:-1,j])/(self.depth[1:,j]-self.depth[:-1,j])

        #Steady plug flow without melting thermal model (cf. document from Catherine)
        self.p[j]=-1+m.exp(self.pprime[j])
        self.Tf[j]=Tf(rhog*ggrav*self.thkie[j])     #FIXME: take into account temporal variations of ice thickness
        self.Tm[j]=(self.Ts+self.Tf[j])/2       #We assume first that we have melting point everywhere
        self.alpha[j]=m.sqrt(self.a[j]/365./24./3600./self.thk[j]/Kg(self.Tm[j], 1.)*rhog*cg(self.Tm[j])/2.)
        self.DeltaT[j]=self.Ts-self.Tf[j]
        self.G[j]=-self.DeltaT[j]*2*Kg(self.Tm[j], 1.)*self.alpha[j]/m.sqrt(m.pi)/erf(self.thkie[j]*self.alpha[j])
        self.is_fusion[j]=(self.G0[j]>self.G[j])
        self.G[j]=min(self.G[j], self.G0[j])
        self.m[j]=(self.G0[j]-self.G[j])*365.242*24*3600/rhog/Lf
        self.T[:,j]=self.Ts-self.G[j]*m.sqrt(m.pi)/2./Kg(self.Tm[j], 1.)/self.alpha[j]*(erf(self.alpha[j]*self.zeta[:,j]*self.thkie[j])-erf(self.alpha[j]*self.thkie[j]))
        self.T_anal[:,j]=self.T[:,j]+0.

        #Mechanical model
        self.mu[j]=self.m[j]/self.a[j]
        self.omega_D[:,j]=1-(self.p[j]+2)/(self.p[j]+1)*(1-self.zetaie[:,j])+1/(self.p[j]+1)*(1-self.zetaie[:,j])**(self.p[j]+2)	#Parrenin et al. (CP, 2007a) 2.2 (3)
        self.omega[:,j]=self.s[j]*self.zetaie[:,j]+(1-self.s[j])*self.omega_D[:,j]   #Parrenin et al. (CP, 2007a) 2.2 (2)
        self.tau[:,j]=(1-self.mu[j])*self.omega[:,j]+self.mu[j]
        self.uz[:,j]=-self.a[j]/365.242/24/3600*self.tau[:,j]

        for it in range(self.tm_iter):

            #Steady non-plug flow thermal model
            self.G=np.empty_like(self.distance)
            if self.is_fusion[j]:
                self.T[0,j]=self.Ts
                self.T[-1,j]=self.Tf[j]
                self.Btemp=np.transpose(np.zeros(np.size(self.zetagrid)-2))
                self.Atemp=1/(self.dzeta*self.thk[j])*(np.diag(Kg((self.T[1:-2,j]+self.T[2:-1,j])/2, self.D[1:-1,j]),1)+np.diag(Kg((self.T[1:-2,j]+self.T[2:-1,j])/2, self.D[1:-1,j]),-1)+np.diag(-Kg((self.T[:-2,j]+self.T[1:-1,j])/2, self.D[:-1,j])-Kg((self.T[1:-1,j]+self.T[2:,j])/2, self.D[1:,j]),0))
                self.Atemp=self.Atemp-0.5*rhog*np.transpose([(self.D[:-1,j]+self.D[1:,j])/2*self.uz[1:-1,j]*cg(self.T[1:-1,j])])*(np.diag(np.ones(np.size(self.zetagrid)-3),-1)-np.diag(np.ones(np.size(self.zetagrid)-3),1))   #TODO:We are missing a term due to heat production by deformatio
                self.Btemp[0]=-(Kg((self.T[0,j]+self.T[1,j])/2, self.D[0,j])/self.dzeta/self.thk[j]-0.5*rhog*(self.D[0,j]+self.D[1,j])/2*cg(self.T[1,j])*self.uz[1,j])*self.Ts
                self.Btemp[-1]=-(Kg((self.T[-2,j]+self.T[-1,j])/2, self.D[-1,j])/self.dzeta/self.thk[j]+0.5*rhog*(self.D[-1,j]+self.D[-2,j])/2*cg(self.T[-2,j])*self.uz[-2,j])*self.Tf[j]
                self.T[1:-1,j]=np.linalg.solve(self.Atemp,self.Btemp)
                self.G[j]=-Kg((self.T[-1,j]+self.T[-2,j])/2, self.D[-1,j])*(self.T[-2,j]-self.T[-1,j])/(self.depthie[-1,j]-self.depthie[-2,j])
                self.m[j]=(self.G0[j]-self.G[j])*365.242*24*3600/rhog/Lf
                if self.G0[j]<=self.G[j]:
                    self.is_fusion[j]=False
            else:
                self.T[0,j]=self.Ts
                self.Btemp=np.transpose(np.zeros(np.size(self.zetagrid)-1))
                self.Atemp=1/(self.dzeta*self.thk[j])*(np.diag(Kg((self.T[1:-2,j]+self.T[2:-1,j])/2, self.D[1:-1,j]),1)+np.diag(Kg((self.T[1:-2,j]+self.T[2:-1,j])/2, self.D[1:-1,j]),-1)+np.diag(-Kg((self.T[:-2,j]+self.T[1:-1,j])/2, self.D[:-1,j])-Kg((self.T[1:-1,j]+self.T[2:,j])/2, self.D[1:,j]),0))
                self.Atemp=self.Atemp-0.5*rhog*np.transpose([(self.D[:-1,j]+self.D[1:,j])/2*self.uz[1:-1,j]*cg(self.T[1:-1,j])])*(np.diag(np.ones(np.size(self.zetagrid)-3),-1)-np.diag(np.ones(np.size(self.zetagrid)-3),1))   #TODO:We are missing a term due to heat production by deformation
                vector=np.concatenate((np.zeros(np.size(self.zetagrid)-3), np.array([ Kg((self.T[-2,j]+self.T[-1,j])/2, self.D[-1,j])/self.dzeta/self.thk[j] - 0.5*rhog*(self.D[-1,j]+self.D[-2,j])/2*cg(self.T[-2,j])*self.uz[-2,j] ]) ))
                self.Atemp=np.concatenate(( self.Atemp, np.array([ vector ]).T ), axis=1)
                vector=np.concatenate(( np.zeros(np.size(self.zetagrid)-3), Kg((self.T[-2,j]+self.T[-1,j])/2, self.D[-1,j])/self.dzeta/self.thk[j] * np.array([1,-1]) ))
                self.Atemp=np.concatenate(( self.Atemp, np.array([ vector ])  ), axis=0)

                self.Btemp[0]=-(Kg((self.T[0,j]+self.T[1,j])/2, self.D[0,j])/self.dzeta/self.thk[j]-0.5*(self.D[0,j]+self.D[1,j])/2*rhog*cg(self.T[0,j])*self.uz[1,j])*self.Ts  #Is there an indice problem here?
                self.Btemp[-1]=-self.G0[j]
                self.T[1:,j]=np.linalg.solve(self.Atemp,self.Btemp)
                self.G[j]=self.G0[j]
                self.m[j]=0.
                if self.T[-1,j]>self.Tf[j]:
                    self.is_fusion[j]=True

            #Mechanical model
            self.mu[j]=self.m[j]/self.a[j]
            self.tau[:,j]=(1-self.mu[j])*self.omega[:,j]+self.mu[j]
            self.uz[:,j]=-self.a[j]/365.242/24/3600*self.tau[:,j]


        self.age_density[:,j]=np.where((self.tau[1:,j]+self.tau[:-1,j])/2>0 ,1/self.a[j]/(self.tau[1:,j]+self.tau[:-1,j])*2,np.nan)
#        toto=np.concatenate(( np.array([[ self.age_surf ]]),np.array([(self.depthie[1:,j]-self.depthie[:-1,j])*self.age_density[:,j]]) ))
        self.agesteady[:,j]=np.cumsum(np.concatenate((np.array([ self.age_surf ]),(self.depthie[1:,j]-self.depthie[:-1,j])*self.age_density[:,j] )), axis=0)

#        for j in range(np.size(self.distance)):
#            self.agesteady[:,j]=np.cumsum(np.concatenate((np.array([self.age_surf]),(self.depthie[1:,j]-self.depthie[:-1,j])/self.a[j]/(self.tau[1:,j]+self.tau[:-1,j])*2)))

#        self.agesteady[0,:]=self.age_surf
#        for i in range(1,np.size(self.zetagrid)):
#            self.agesteady[i,:]=np.where(self.tau[i,:]+self.tau[i-1,:]>0,self.agesteady[i-1,:]+(self.depthie[i,:]-self.depthie[i-1,:])/self.a/(self.tau[i,:]+self.tau[i-1,:])*2,1000000)
#            self.agesteady[i,:]=np.where(self.agesteady[i,:]<1000000,self.agesteady[i,:],1000000)

        self.age[:,j]=np.interp(self.agesteady[:,j],np.concatenate((np.array([-1000000000]),self.AICC2012_steadyage,np.array([1e9*self.AICC2012_steadyage[-1]]))),np.concatenate((np.array([self.AICC2012_age[0]]),self.AICC2012_age,np.array([1e9*self.AICC2012_age[-1]]))))

        f=interp1d(self.depth[:,j],self.age[:,j])
        self.agebot[j]=f(max(self.depth[:,j])-60)
        self.age100m[j]=f(max(self.depth[:,j])-100)
        self.age150m[j]=f(max(self.depth[:,j])-150)
        self.age200m[j]=f(max(self.depth[:,j])-200)
        self.age250m[j]=f(max(self.depth[:,j])-250)
        self.hor_modage[:,j]=f(self.hor[:,j])
        self.iso_modage[:,j]=np.interp(self.iso[:,j],self.depth[:,j],self.age[:,j])
        h1=interp1d(self.age[:-1,j],self.age_density[:,j])
        h2=interp1d(self.age[:,j],self.depth[:,j])
        if self.agebot[j]>=1000000:
            self.age_density1Myr[j]=h1(1000000)
        else:
            self.age_density1Myr[j]=np.nan
        if self.agebot[j]>=1200000:
            self.age_density1dot2Myr[j]=h1(1200000)
        else:
            self.age_density1dot2Myr[j]=np.nan
        if self.agebot[j]>=1500000:
            self.age_density1dot5Myr[j]=h1(1500000)
        else:
            self.age_density1dot5Myr[j]=np.nan
        if max(self.age[:,j])>=600000:
            self.height0dot6Myr[j]=self.thk[j]-h2(600000)
            self.twtt0dot6Myr[j]=(h2(600000)-self.firn_correction)*100/84.248+250.
        else:
            self.height0dot6Myr[j]=np.nan
            self.twtt0dot6Myr[j]=-98765.0
        if max(self.age[:,j])>=800000:
            self.height0dot8Myr[j]=self.thk[j]-h2(800000)
            self.twtt0dot8Myr[j]=(h2(800000)-self.firn_correction)*100/84.248+250.
        else:
            self.height0dot8Myr[j]=np.nan
            self.twtt0dot8Myr[j]=-98765.0
        if max(self.age[:,j])>=1000000:
            self.height1Myr[j]=self.thk[j]-h2(1000000)
            self.twtt1Myr[j]=(h2(1000000)-self.firn_correction)*100/84.248+250.
        else:
            self.height1Myr[j]=np.nan
            self.twtt1Myr[j]=-98765.0
        if max(self.age[:,j])>=1200000:
            self.height1dot2Myr[j]=self.thk[j]-h2(1200000)
            self.twtt1dot2Myr[j]=(h2(1200000)-self.firn_correction)*100/84.248+250.
        else:
            self.height1dot2Myr[j]=np.nan
            self.twtt1dot2Myr[j]=-98765.0
        if max(self.age[:,j])>=1500000:
            self.height1dot5Myr[j]=self.thk[j]-h2(1500000)
            self.twtt1dot5Myr[j]=(h2(1500000)-self.firn_correction)*100/84.248+250.
        else:
            self.height1dot5Myr[j]=np.nan
            self.twtt1dot5Myr[j]=-98765.0

        self.twttBed[j]=(self.thk[j]-self.firn_correction)*100/84.248+250.  #TODO: make a function to convert to twtt, and make an array for the different isochrones.


        return np.concatenate(( np.array([self.a[j]]),np.array([self.m[j]]),np.array([self.pprime[j]]),self.age[:,j],np.log(self.age[1:,j]-self.age_surf),np.array([self.G0[j]]) ))

    def model(self):  #TODO: kill this or make a call to model(j)
        for j in range(np.size(self.distance)):
            self.model1D(j)

        return np.concatenate((self.a,self.m,self.pprime,self.age.flatten(),self.G0))

#Residuals function

    def residuals1D(self, variables1D, j):
        var=variables1D
        self.a[j]=var[0]
        var=np.delete(var,[0])
#        print 'a: ',self.a[j]
#        self.m=variables[np.size(self.distance):2*np.size(self.distance)]
        self.pprime[j]=var[0]
        var=np.delete(var,[0])
        if self.invert_G0:
            self.G0[j]=var[0]
            var=np.delete(var,[0])
        if self.invert_thk:
            self.thk[j]=var[0]
            var=np.delete(var,[0])

        self.model1D(j)
        f=interp1d(self.depth[:,j],self.age[:,j])  #Is this useful?
#        print self.age[:,j]
#        print self.iso_modage[:,j]
#        print 'shape iso_age and iso_sigma:', np.shape(self.iso_age),np.shape(self.iso_sigma)
        resi=(self.iso_age.flatten()-self.iso_modage[:,j])/self.iso_sigma.flatten()
#        print resi
        resi=resi[np.where(~np.isnan(resi))]
        resi=np.concatenate((resi,np.array([ (self.pprime[j]-self.pprime_prior)/self.pprime_sigma ]) ))
        if self.invert_G0:
            resi=np.concatenate((resi, np.array([ (self.G0[j]-self.G0_prior)/self.G0_sigma ]) ))
#        resi=np.append(resi,(self.pprime[j]-self.pprime_prior)/self.pprime_sigma)
#        if self.invert_G0:
#            resi=np.append(resi,(self.G0[j]-self.G0_prior)/self.G0_sigma)
#        print 'Residual: ',resi
        return resi

    def cost_fct(self, variables1D, j):

        res=self.residuals1D(variables1D, j)
#        cost=1.-m.exp( -np.sum(res**2)/2. )
        cost=np.sum(res**2)/2.
        return cost


    def residuals(self, variables):
        var=variables
        self.a=var[0:np.size(self.distance)]
        var=np.delete[var,np.zeros_like(self.distance)]
#        self.m=variables[np.size(self.distance):2*np.size(self.distance)]
        self.pprime=var[0:np.size(self.distance)]
        var=np.delete[var,np.zeros_like(self.distance)]
        if self.invert_G0:
            self.G0=var[0:np.size(self.distance)]
            var=np.delete[var,np.zeros_like(self.distance)]
        if self.invert_G0:
            self.thk=var[0:np.size(self.distance)]
            var=np.delete[var,np.zeros_like(self.distance)]

        age=self.model()
        for j in range(np.size(self.distance)):
            f=interp1d(self.depth[:,j],self.age[:,j])
            self.iso_modage[:,j]=f(self.iso[:,j])
        resi=(self.iso_age-self.iso_modage)/self.iso_sigma
        iso_age_flatten=self.iso_age.flatten()
        resi=resi.flatten()
        resi=resi[np.where(~np.isnan(resi))]
        resi=np.concatenate((resi,(self.pprime-self.pprime_prior)/self.pprime_sigma))
        if self.invert_G0:
            resi=np.concatenate((resi,(self.G0-self.G0_prior)/self.G0_sigma))
        return resi


#    def jacobian1D(self, j):
#        epsilon=np.sqrt(np.diag(self.hess1D))/100000000.
#        model0=self.model1D(j)
#        jacob=np.empty((np.size(self.variables1D), np.size(model0)))
#        for i in np.arange(np.size(self.variables1D)):
#            self.variables1D[i]=self.variables1D[i]+epsilon[i]
#            self.residuals1D(self.variables1D, j)
#            model1=self.model1D(j)
#            jacob[i]=(model1-model0)/epsilon[i]
#            self.variables1D[i]=self.variables1D[i]+epsilon[i]
#        self.residuals1D(self.variables1D, j)
#       return jacob

    def jacobian1D(self, j):
        epsilon=np.sqrt(np.diag(self.hess1D))/10000000000.
        model0=self.model1D(j)
        jacob=np.empty((np.size(self.variables1D), np.size(model0)))
        for i in np.arange(np.size(self.variables1D)):
            self.variables1D[i]=self.variables1D[i]+epsilon[i]
            self.residuals1D(self.variables1D, j)
            model1=self.model1D(j)
            self.variables1D[i]=self.variables1D[i]-epsilon[i]
            self.residuals1D(self.variables1D, j)
            model2=self.model1D(j)
            jacob[i]=(model1-model2)/2./epsilon[i]
            self.variables1D[i]=self.variables1D[i]+epsilon[i]
        self.residuals1D(self.variables1D, j)

        return jacob


    def jacobian(self):
        epsilon=np.sqrt(np.diag(self.hess))/100000000.
        model0=self.model()
        jacob=np.empty((np.size(self.variables), np.size(model0)))
        for i in np.arange(np.size(self.variables)):
            self.variables[i]=self.variables[i]+epsilon[i]
            self.residuals(self.variables)
            model1=self.model()
            jacob[i]=(model1-model0)/epsilon[i]
            self.variables[i]=self.variables[i]-epsilon[i]
        model0=self.model()

        return jacob



    def accu_layers(self):
        self.accusteady_layer=np.zeros((self.nbiso,np.size(self.distance)))
        self.accu_layer=np.zeros((self.nbiso,np.size(self.distance)))
        for j in range(np.size(self.distance)):
            f=interp1d(self.depth[:,j],self.age[:,j])
            self.iso_modage[:,j]=f(self.iso[:,j])
            self.accusteady_layer[0,j]=self.a[j]*(self.iso_modage[0,j]-self.age_surf)/(self.iso_age[0]-self.age_surf)
            self.accusteady_layer[1:,j]=self.a[j]*(self.iso_modage[1:,j]-self.iso_modage[:-1,j])/(self.iso_age[1:]-self.iso_age[:-1]).flatten()
        self.accu_layer[0,]=self.accusteady_layer[0,:]*(self.iso_steadyage[0]-self.age_surf)/(self.iso_age[0]-self.age_surf)
        self.accu_layer[1:,]=self.accusteady_layer[1:,]*(self.iso_steadyage[1:]-self.iso_steadyage[:-1])/(self.iso_age[1:]-self.iso_age[:-1])

        return

    def sigma1D(self,j):
        jacob=self.jacobian1D(j)


        index=0
        c_model=np.dot(np.transpose(jacob[:,index:index+1]),np.dot(self.hess1D,jacob[:,index:index+1]))
        self.sigma_a[j]=np.sqrt(np.diag(c_model))[0]
        index=index+1
        c_model=np.dot(np.transpose(jacob[:,index:index+1]),np.dot(self.hess1D,jacob[:,index:index+1]))
        self.sigma_m[j]=np.sqrt(np.diag(c_model))[0]
        index=index+1
        c_model=np.dot(np.transpose(jacob[:,index:index+1]),np.dot(self.hess1D,jacob[:,index:index+1]))
        self.sigma_pprime[j]=np.sqrt(np.diag(c_model))[0]
        index=index+1
        c_model=np.dot(np.transpose(jacob[:,index:index+np.size(self.age[:,j])]),np.dot(self.hess1D,jacob[:,index:index+np.size(self.age[:,j])]))
        self.sigma_age[:,j]=np.sqrt(np.diag(c_model))
#        print np.size(self.sigma_age)
        index=index+np.size(self.age[:,j])
        c_model=np.dot(np.transpose(jacob[:,index:index+np.size(self.age[1:,j])]),np.dot(self.hess1D,jacob[:,index:index+np.size(self.age[1:,j])]))
        self.sigma_logage[1:,j]=np.sqrt(np.diag(c_model))
        self.sigma_logage[0,j]=np.nan
        index=index+np.size(self.age[1:,j])
        c_model=np.dot(np.transpose(jacob[:,index:index+1]),np.dot(self.hess1D,jacob[:,index:index+1]))
        self.sigma_G0[j]=np.sqrt(np.diag(c_model))[0]

        f=interp1d(self.depth[:,j],self.sigma_age[:,j])
        self.sigmabotage[j]=f(self.thk[j]-60.)
        self.iso_modage_sigma[:,j]=f(self.iso[:,j])

        return

    def sigma(self):
        jacob=self.jacobian()

        index=0
        c_model=np.dot(np.transpose(jacob[:,index:index+np.size(self.a)]),np.dot(self.hess,jacob[:,index:index+np.size(self.a)]))
        self.sigma_a=np.sqrt(np.diag(c_model))
        index=index+np.size(self.a)
        c_model=np.dot(np.transpose(jacob[:,index:index+np.size(self.m)]),np.dot(self.hess,jacob[:,index:index+np.size(self.m)]))
        self.sigma_m=np.sqrt(np.diag(c_model))
        index=index+np.size(self.m)
        c_model=np.dot(np.transpose(jacob[:,index:index+np.size(self.p)]),np.dot(self.hess,jacob[:,index:index+np.size(self.p)]))
        self.sigma_pprime=np.sqrt(np.diag(c_model))
        index=index+np.size(self.p)
        c_model=np.dot(np.transpose(jacob[:,index:index+np.size(self.age)]),np.dot(self.hess,jacob[:,index:index+np.size(self.age)]))
        self.sigma_age=np.sqrt(np.diag(c_model))
#        print np.size(self.sigma_age)
        self.sigma_age=np.reshape(self.sigma_age,(np.size(self.zetagrid),np.size(self.distance)))
        index=index+np.size(self.age)
        c_model=np.dot(np.transpose(jacob[:,index:index+np.size(self.G0)]),np.dot(self.hess,jacob[:,index:index+np.size(self.G0)]))
        self.sigma_G0=np.sqrt(np.diag(c_model))

        return

#Plotting the raw and interpolated radar datasets
    def data_display(self):
        plt.figure('Data')
        plt.plot(self.distance_raw, self.thk_raw, label='raw bedrock', color='0.5', linewidth=2)
        plt.plot(self.distance, self.thk, label='interpolated bedrock', color='k', linewidth=2)
        for i in range(self.nbiso):
            if i==0:
                plt.plot(self.distance_raw, self.iso_raw[i,:], color='c', label='raw isochrones')
                plt.plot(self.distance, self.iso[i,:], color='b', label='interpolated isochrones')
            else:
                plt.plot(self.distance_raw, self.iso_raw[i,:], color='c')
                plt.plot(self.distance, self.iso[i,:], color='b')
        for i in range(self.nbhor):
            if i==0:
                plt.plot(self.distance_raw, self.hor_raw[i,:], color='y', label='raw horizons')
                plt.plot(self.distance, self.hor[i,:], color='g', label='interpolated horizons')
            elif i>0 and i<self.nbhor-self.nbdsz:
                plt.plot(self.distance_raw, self.hor_raw[i,:], color='y')
                plt.plot(self.distance, self.hor[i,:], color='g')
            elif i==self.nbhor-self.nbdsz:
                plt.plot(self.distance_raw, self.hor_raw[i,:], color='orange', label='raw DSZ')
                plt.plot(self.distance, self.hor[i,:], color='r', label='interpolated DSZ')
            else:
                plt.plot(self.distance_raw, self.hor_raw[i,:], color='orange')
                plt.plot(self.distance, self.hor[i,:], color='r')
        if self.is_EDC:
            EDC_x=np.array([self.distance_EDC, self.distance_EDC])
            EDC_y=np.array([0., 3200.])
            if self.EDC_line_dashed==True:
                plt.plot(EDC_x, EDC_y, label='EDC ice core', color='r', linewidth=2, linestyle='--')
            else:
                plt.plot(EDC_x, EDC_y, label='EDC ice core', color='r', linewidth=2)
        if self.is_NESW:
            plt.xlabel('<NE - distance (km) - SW>')
        else:
            plt.xlabel('distance (km)')
        plt.ylabel('depth (m)')
#        plt.legend(loc=1)
        x1,x2,y1,y2 = plt.axis()
        plt.axis((x1,x2,y2,0))
        if self.reverse_distance:
            plt.gca().invert_xaxis()
        pp=PdfPages(self.label+'Data.pdf')
        pp.savefig(plt.figure('Data'))
        pp.close()

#Plot of the model results

    def model_display(self):

        fig=plt.figure('Model steady')
        plotmodel = fig.add_subplot(111, aspect=self.aspect)
        plt.plot(self.distance, self.thkreal, label='obs. bedrock', color='k', linewidth=2)
        plt.fill_between(self.distance, self.thkreal, self.thk, color='0.5', label='stagnant ice')
        for i in range(self.nbiso):
            if i==0:
                plt.plot(self.distance, self.iso[i,:], color='w', label='obs. isochrones')
            else:
                plt.plot(self.distance, self.iso[i,:], color='w')
        for i in range(self.nbhor):
            if i==0:
                plt.plot(self.distance, self.hor[i,:], color='0.5', label='obs. horizons')
            elif i>0 and i<self.nbhor-self.nbdsz:
                plt.plot(self.distance, self.hor[i,:], color='0.5')
            elif i==self.nbhor-self.nbdsz:
                plt.plot(self.distance, self.hor[i,:], color='r', label='obs. DSZ')
            else:
                plt.plot(self.distance, self.hor[i,:], color='r')
        levels=np.arange(0, 1600, 100)
        levels_color=np.arange(0, 1500, 10)
        plt.contourf(self.dist, self.depth, self.agesteady/1000., levels_color)
        if self.is_EDC:
            EDC_x=np.array([self.distance_EDC, self.distance_EDC])
            EDC_y=np.array([0., 3200.])
            plt.plot(EDC_x, EDC_y, label='EDC ice core', color='r', linewidth=2)
        if self.is_NESW:
            plt.xlabel('<NE - distance (km) - SW>')
        else:
            plt.xlabel('distance (km)')
        plt.ylabel('depth (m)')
#        plt.legend(loc=2)
        cb=plt.colorbar()
        cb.set_ticks(levels)
        cb.set_ticklabels(levels)
        cb.set_label('Modeled steady age (kyr)')
        x1,x2,y1,y2 = plt.axis()
        if self.max_depth=='auto':
            self.max_depth=y2
        plt.axis((min(self.distance),max(self.distance),self.max_depth,0))
        if self.reverse_distance:
            plt.gca().invert_xaxis()
        pp=PdfPages(self.label+'Model-steady.pdf')
        pp.savefig(plt.figure('Model steady'))
        pp.close()


        fig=plt.figure('Model')
        plotmodel = fig.add_subplot(111, aspect=self.aspect)
        plt.plot(self.distance, self.thkreal, color='k', linewidth=2, label='bed')
        plt.fill_between(self.distance, self.thkreal, self.thk, color='0.5', label='stagnant ice')
#        plt.legend(loc=1)
        for i in range(self.nbiso):
            if i==0:
                plt.plot(self.distance, self.iso[i,:], color='w', label='obs. isochrones')
            else:
                plt.plot(self.distance, self.iso[i,:], color='w')
        for i in range(self.nbhor):
            if i==0:
                plt.plot(self.distance, self.hor[i,:], color='0.5', label='obs. horizons')
            elif i>0 and i<self.nbhor-self.nbdsz:
                plt.plot(self.distance, self.hor[i,:], color='0.5')
            elif i==self.nbhor-self.nbdsz:
                plt.plot(self.distance, self.hor[i,:], color='r', label='obs. DSZ')
            else:
                plt.plot(self.distance, self.hor[i,:], color='r')
        levels=np.arange(0, 1600, 100)
        levels_color=np.arange(0, 1500, 10)
        plt.contourf(self.dist, self.depth, self.age/1000., levels_color)
        if self.is_EDC:
            EDC_x=np.array([self.distance_EDC, self.distance_EDC])
            EDC_y=np.array([0., 3200.])
            if self.EDC_line_dashed==True:
                plt.plot(EDC_x, EDC_y, label='EDC ice core', color='r', linewidth=2, linestyle='--')
            else:
                plt.plot(EDC_x, EDC_y, label='EDC ice core', color='r', linewidth=2)
        if self.is_NESW:
            plt.xlabel('<NE - distance (km) - SW>')
        else:
            plt.xlabel('distance (km)')
        plt.ylabel('depth (m)')
        if self.is_legend:
            leg=plt.legend(loc=1)
            frame=leg.get_frame()
            frame.set_facecolor('0.75')
        cb=plt.colorbar()
        cb.set_ticks(levels)
        cb.set_ticklabels(levels)
        cb.set_label('Modeled age (kyr)')
        x1,x2,y1,y2 = plt.axis()
        plt.axis((min(self.distance),max(self.distance),self.max_depth,0))
        if self.reverse_distance:
            plt.gca().invert_xaxis()
        if self.settick=='manual':
            plotmodel.set_xticks(np.arange(self.min_tick,self.max_tick+1.,self.delta_tick))
        pp=PdfPages(self.label+'Model.pdf')
        pp.savefig(plt.figure('Model'))
        pp.close()

        fig=plt.figure('Age misfit')
        plotmodel = fig.add_subplot(111, aspect=self.aspect)
        plt.plot(self.distance, self.thkreal, color='k', linewidth=2)
        plt.fill_between(self.distance, self.thkreal, self.thk, color='0.5', label='stagnant ice')
        norm=Normalize(vmin=-5000, vmax=5000)
        for i in range(self.nbiso):
            colorscatter=self.iso_modage[i,:]-self.iso_age[i]
#            print i, np.shape(colorscatter),np.shape(self.iso[i,:]),colorscatter
            if i==0:
                sc=plt.scatter(self.distance, self.iso[i,:], c=colorscatter, label='obs. isochrones', s=7, edgecolor='', norm=norm)
            else:
                plt.scatter(self.distance, self.iso[i,:], c=colorscatter, s=7, edgecolor='', norm=norm)
#        levels=np.arange(0, 1600000, 100000)
#        levels_color=np.arange(0, 1500000, 10000)
#        plt.contourf(self.dist, self.depth, self.age, levels_color)
        if self.is_EDC:
            EDC_x=np.array([self.distance_EDC, self.distance_EDC])
            EDC_y=np.array([0., 3200.])
            if self.EDC_line_dashed==True:
                plt.plot(EDC_x, EDC_y, label='EDC ice core', color='r', linewidth=2, linestyle='--')
            else:
                plt.plot(EDC_x, EDC_y, label='EDC ice core', color='r', linewidth=2)
        if self.is_NESW:
            plt.xlabel('<NE - distance (km) - SW>')
        else:
            plt.xlabel('distance (km)')
        plt.ylabel('depth (m)')
        if self.is_legend:
            print 'test'
            leg=plt.legend(loc=1)
            frame=leg.get_frame()
            frame.set_facecolor('0.75')
        cb=plt.colorbar(sc)
#        cb.set_ticks(levels)
#        cb.set_ticklabels(levels)
        cb.set_label('Age misfit (yr)')
        x1,x2,y1,y2 = plt.axis()
        plt.axis((min(self.distance),max(self.distance),self.max_depth,0))
        if self.reverse_distance:
            plt.gca().invert_xaxis()
        if self.settick=='manual':
            plotmodel.set_xticks(np.arange(self.min_tick,self.max_tick+1.,self.delta_tick))
        pp=PdfPages(self.label+'AgeMisfit.pdf')
        pp.savefig(plt.figure('Age misfit'))
        pp.close()


        fig = plt.figure('Model confidence interval')
        plotmodelci = fig.add_subplot(111, aspect=self.aspect)
        plt.plot(self.distance, self.thkreal, color='k', linewidth=2)
        plt.fill_between(self.distance, self.thkreal, self.thk, color='0.5', label='stagnant ice')
        for i in range(self.nbiso):
            if i==0:
                plt.plot(self.distance, self.iso[i,:], color='w', label='obs. isochrones')
            else:
                plt.plot(self.distance, self.iso[i,:], color='w')
        for i in range(self.nbhor):
            if i==0:
                plt.plot(self.distance, self.hor[i,:], color='0.5', label='obs. horizons')
            elif i>0 and i<self.nbhor-self.nbdsz:
                plt.plot(self.distance, self.hor[i,:], color='0.5')
            elif i==self.nbhor-self.nbdsz:
                plt.plot(self.distance, self.hor[i,:], color='r', label='obs. DSZ')
            else:
                plt.plot(self.distance, self.hor[i,:], color='r')
#        levels=np.arange(0, 200000, 20000)
        levels_log=np.arange(2, 6, 0.1)
        levels=np.power(10, levels_log)
#        levels=np.array([100, 250, 500, 1000, 2500, 5000, 10000, 25000, 50000, 100000, 250000])
        plt.contourf(self.dist, self.depth, self.sigma_age, levels, norm = LogNorm())
        cb=plt.colorbar()
        cb.set_label('Modeled age confidence interval (yr)')
#        levels_labels=np.where( np.equal(np.mod(np.arange(2,6,0.1), 1), 0) , np.power(10., np.arange(2,6,0.1)), "" ) 
        levels_labels=np.array([])
        for i in np.arange(2,6,1):
            levels_labels=np.concatenate((levels_labels, np.array([10**i, '', '', '', '', '', '', '', '']) ))
        cb.set_ticklabels(levels_labels)
        levels_ticks=np.concatenate(( np.arange(100, 1000, 100), np.arange(1000, 10000, 1000), np.arange(10000, 100000, 10000), np.arange(100000, 600000, 100000) )) 
        cb.set_ticks(levels_ticks)
        if self.is_EDC:
            EDC_x=np.array([self.distance_EDC, self.distance_EDC])
            EDC_y=np.array([0., 3200.])
            if self.EDC_line_dashed==True:
                plt.plot(EDC_x, EDC_y, label='EDC ice core', color='r', linewidth=2, linestyle='--')
            else:
                plt.plot(EDC_x, EDC_y, label='EDC ice core', color='r', linewidth=2)
        if self.is_NESW:
            plt.xlabel('<NE - distance (km) - SW>')
        else:
            plt.xlabel('distance (km)')
        plt.ylabel('depth (m)')
        if self.is_legend:
            leg=plt.legend(loc=1)
            frame=leg.get_frame()
            frame.set_facecolor('0.75')
        x1,x2,y1,y2 = plt.axis()
        plt.axis((min(self.distance),max(self.distance),self.max_depth,0))
        if self.reverse_distance:
            plt.gca().invert_xaxis()
        if self.settick=='manual':
            plotmodelci.set_xticks(np.arange(self.min_tick,self.max_tick+1.,self.delta_tick))
        pp=PdfPages(self.label+'Model-confidence-interval.pdf')
        pp.savefig(plt.figure('Model confidence interval'))
        pp.close()


        plt.figure('Thinning')
        plt.plot(self.distance, self.thkreal, label='obs. bedrock', color='k', linewidth=2)
        plt.fill_between(self.distance, self.thkreal, self.thk, color='0.5', label='stagnant ice')
        for i in range(self.nbiso):
            if i==0:
                plt.plot(self.distance, self.iso[i,:], color='k', label='obs. isochrones')
            else:
                plt.plot(self.distance, self.iso[i,:], color='k')
        for i in range(self.nbhor):
            if i==0:
                plt.plot(self.distance, self.hor[i,:], color='0.5', label='obs. horizons')
            elif i>0 and i<self.nbhor-self.nbdsz:
                plt.plot(self.distance, self.hor[i,:], color='0.5')
            elif i==self.nbhor-self.nbdsz:
                plt.plot(self.distance, self.hor[i,:], color='r', label='obs. DSZ')
            else:
                plt.plot(self.distance, self.hor[i,:], color='r')
        plt.contourf(self.dist, self.depth, self.tau)
        if self.is_EDC:
            EDC_x=np.array([self.distance_EDC, self.distance_EDC])
            EDC_y=np.array([0., 3200.])
            plt.plot(EDC_x, EDC_y, label='EDC ice core', color='r', linewidth=2)
        if self.is_NESW:
            plt.xlabel('<NE - distance (km) - SW>')
        else:
            plt.xlabel('distance (km)')
        plt.ylabel('depth (m)')
        plt.legend(loc=2)
        cb=plt.colorbar()
        cb.set_label('Modeled thinning')
        x1,x2,y1,y2 = plt.axis()
        plt.axis((min(self.distance),max(self.distance),self.max_depth,0))
        if self.reverse_distance:
            plt.gca().invert_xaxis()
        pp=PdfPages(self.label+'Thinning.pdf')
        pp.savefig(plt.figure('Thinning'))
        pp.close()

        fig=plt.figure('Temperature')
        plt.plot(self.distance, self.thkreal, label='obs. bedrock', color='k', linewidth=2)
        plt.fill_between(self.distance, self.thkreal, self.thk, color='0.5', label='stagnant ice')
        plt.plot(self.distance, np.where(self.is_fusion,np.nan,self.thk), color='b', linewidth=4)
        plt.contourf(self.dist, self.depth, self.T)
        if self.is_EDC:
            EDC_x=np.array([self.distance_EDC, self.distance_EDC])
            EDC_y=np.array([0., 3200.])
            plt.plot(EDC_x, EDC_y, label='EDC ice core', color='r', linewidth=2)
        if self.is_NESW:
            plt.xlabel('<NE - distance (km) - SW>')
        else:
            plt.xlabel('distance (km)')
        plt.ylabel('depth (m)')
        plt.legend(loc=2)
        cb=plt.colorbar()
        cb.set_label('Modeled temperature (K)')
        x1,x2,y1,y2 = plt.axis()
        plt.axis((min(self.distance),max(self.distance),self.max_depth,0))
        if self.reverse_distance:
            plt.gca().invert_xaxis()
        pp=PdfPages(self.label+'Temperature.pdf')
        pp.savefig(plt.figure('Temperature'))
        pp.close()
#        plt.close(fig)



#        fig = plt.figure('Accumulation history')
        lines=[zip(self.distance, 917*self.accu_layer[i,:]) for i in range(self.nbiso)]
        z=(self.iso_age.flatten()[1:]+self.iso_age.flatten()[:-1])/2
        z=np.concatenate(( np.array([(self.age_surf+self.iso_age.flatten()[0])/2]) , z ))
        fig, ax = plt.subplots()
        lines = LineCollection(lines, array=z, cmap=plt.cm.rainbow, linewidths=2)
        ax.add_collection(lines)
        ax.autoscale()
        cb=fig.colorbar(lines)
        cb.set_label('Average layer age (yr)')
        if self.is_NESW:
            plt.xlabel('<NE - distance (km) - SW>')
        else:
            plt.xlabel('distance (km)')
        plt.ylabel('steady accumulation (mm-we/yr)')
        if self.reverse_distance:
            plt.gca().invert_xaxis()


        pp=PdfPages(self.label+'AccumulationHistory.pdf')
        pp.savefig(fig)
        pp.close()
#        plt.close(fig)



#Plot of the parameters

    def parameters_display(self):
        f=plt.figure('Parameters')
#        f=plt.figure('Parameters', figsize=(4,6))
        plotpara = plt.subplot(311, aspect=self.aspect/(self.accu_max-self.accu_min)*self.max_depth/3)
        plt.plot(self.distance, self.a*100, label='accumulation', color='b')
        plt.plot(self.distance, (self.a-self.sigma_a)*100, color='b', linestyle='--')
        plt.plot(self.distance, (self.a+self.sigma_a)*100, color='b', linestyle='--')
        plt.ylabel('accu. (cm/yr)')
        x1,x2,y1,y2 = plt.axis()
        plt.axis((min(self.distance),max(self.distance), self.accu_min, self.accu_max))
        if self.reverse_distance:
            plt.gca().invert_xaxis()
        if self.is_EDC:
            EDC_x=np.array([self.distance_EDC, self.distance_EDC])
            EDC_y=np.array([self.accu_min, self.accu_max])
            if self.EDC_line_dashed==True:
                plt.plot(EDC_x, EDC_y, label='EDC ice core', color='r', linewidth=2, linestyle='--')
            else:
                plt.plot(EDC_x, EDC_y, label='EDC ice core', color='r', linewidth=2)
#        plt.legend()
        if self.settick=='manual':
            plotpara.set_xticks(np.arange(self.min_tick,self.max_tick+1.,self.delta_tick))
        plotpara=plt.subplot(312, aspect=self.aspect/(self.melting_max-self.melting_min)*self.max_depth/3)
        plt.plot(self.distance, self.m*1000, label='melting', color='r')
        plt.plot(self.distance, (self.m-self.sigma_m)*1000, color='r', linestyle='--')
        plt.plot(self.distance, (self.m+self.sigma_m)*1000, color='r', linestyle='--')
        plt.ylabel('melting (mm/yr)')
        x1,x2,y1,y2 = plt.axis()
        plt.axis((min(self.distance),max(self.distance),self.melting_min,self.melting_max))
        if self.reverse_distance:
            plt.gca().invert_xaxis()
        if self.is_EDC:
            EDC_x=np.array([self.distance_EDC, self.distance_EDC])
            EDC_y=np.array([self.melting_min, self.melting_max])
            if self.EDC_line_dashed==True:
                plt.plot(EDC_x, EDC_y, label='EDC ice core', color='r', linewidth=2, linestyle='--')
            else:
                plt.plot(EDC_x, EDC_y, label='EDC ice core', color='r', linewidth=2)
        plotpara.yaxis.tick_right()
        plotpara.yaxis.set_label_position('right')
        if self.settick=='manual':
            plotpara.set_xticks(np.arange(self.min_tick,self.max_tick+1.,self.delta_tick))
        plotpara=plt.subplot(313, aspect=self.aspect/(m.log(self.p_max+1)-m.log(self.p_min+1))*self.max_depth/3)
        plt.plot(self.distance, self.pprime, label='p', color='g')
        plt.plot(self.distance, self.pprime-self.sigma_pprime, color='g', linestyle='--')
        plt.plot(self.distance, self.pprime+self.sigma_pprime, color='g', linestyle='--')
        plt.ylabel('p+1 parameter')
        if self.is_EDC:
            EDC_x=np.array([self.distance_EDC, self.distance_EDC])
            EDC_y=np.array([m.log(self.p_min), m.log(self.p_max)])
            if self.EDC_line_dashed==True:
                plt.plot(EDC_x, EDC_y, label='EDC ice core', color='r', linewidth=2, linestyle='--')
            else:
                plt.plot(EDC_x, EDC_y, label='EDC ice core', color='r', linewidth=2)
#        plotpara.set_yticks(np.log(np.arange(1.,11.)))
        plotpara.set_yticks(np.log(np.concatenate((np.arange(1.,10.),10.*np.arange(1.,10.) )) ))
        labels=["1", "", "", "", "", "", "", "", "", "10"]
        plotpara.set_yticklabels(labels)
        if self.settick=='manual':
            plotpara.set_xticks(np.arange(self.min_tick,self.max_tick+1.,self.delta_tick))
        x1,x2,y1,y2 = plt.axis()
        plt.axis((min(self.distance),max(self.distance),m.log(self.p_min+1),m.log(self.p_max+1)))
        if self.reverse_distance:
            plt.gca().invert_xaxis()
        if self.is_NESW:
            plt.xlabel('<NE - distance (km) - SW>')
        else:
            plt.xlabel('distance (km)')
#        plt.legend()
        f.subplots_adjust(hspace=0)
        plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
        pp=PdfPages(self.label+'Parameters.pdf')
        pp.savefig(plt.figure('Parameters'))
        pp.close()

        if self.invert_G0:
            plt.figure('Geothermal heat flux')
            plt.plot(self.distance, self.G0*1000, label='G0', color='k')
            plt.plot(self.distance, (self.G0-self.sigma_G0)*1000, color='k', linestyle='--')
            plt.plot(self.distance, (self.G0+self.sigma_G0)*1000, color='k', linestyle='--')
            plt.ylabel('$G_0$ (mW/m$^2$)')
        if self.reverse_distance:
            plt.gca().invert_xaxis()
#            plt.yaxis.tick_right()
#            plt.yaxis.set_label_position('right')
            pp=PdfPages(self.label+'GeothermalHeatFlux.pdf')
            pp.savefig(plt.figure('Geothermal heat flux'))
            pp.close()


    def parameters_save(self):
        output=np.vstack((self.LON, self.LAT, self.distance,self.a, self.sigma_a, self.accu_layer))
        header='#LON\tLAT\tdistance(km)\taccu(ice-m/yr)\tsigma_accu'
        header=header+'\tlayer '+str(int(self.age_surf/1000.))+'-'+str(int(self.iso_age[0]/1000.))+'kyr'
        for i in range(self.nbiso-1):
            header=header+'\tlayer '+str(int(self.iso_age[i]/1000.))+'-'+str(int(self.iso_age[i+1]/1000.))+'kyr'
        header=header+'\n'
        with open(self.label+'a.txt','w') as f:
            f.write(header)
            np.savetxt(f,np.transpose(output), delimiter="\t") 
        output=np.vstack((self.LON, self.LAT, self.distance,self.m, self.sigma_m))
        with open(self.label+'m.txt','w') as f:
            f.write('#LON\tLAT\tdistance(km)\tmelting(ice-m/yr)\tsigma_melting\n')
            np.savetxt(f,np.transpose(output), delimiter="\t") 
        output=np.vstack((self.LON, self.LAT, self.distance,self.pprime, self.sigma_pprime))
        with open(self.label+'pprime.txt','w') as f:
            f.write('#LON\tLAT\tdistance(km)\tpprime\tsigma_pprime\n')
            np.savetxt(f,np.transpose(output), delimiter="\t") 
        output=np.vstack((self.LON, self.LAT, self.distance,self.G0, self.sigma_G0))
        with open(self.label+'G0.txt','w') as f:
            f.write('#LON\tLAT\tdistance(km)\tG0\tsigma_G0\n')
            np.savetxt(f,np.transpose(output), delimiter="\t") 

    def bot_age_save(self):
        output=np.vstack((self.LON,self.LAT,self.distance,self.thk,self.agebot,self.agebotmin,self.age100m,self.age150m,self.age200m,self.age250m,self.age_density1Myr,self.age_density1dot2Myr,self.age_density1dot5Myr, self.height0dot6Myr,self.height0dot8Myr,self.height1Myr,self.height1dot2Myr,self.height1dot5Myr))
        with open(self.label+'agebottom.txt','w') as f:
            f.write('#LON\tLAT\tdistance(km)\tthickness(m)\tage60m(yr-b1950)\tage-min(yr-b1950)\tage100m\tage150m\tage200m\tage250\tage_density1Myr\tage_density1.2Myr\tage_density1.5Myr\theight0.6Myr\theight0.8Myr\theight1Myr\theight1.2Myr\theight1.5Myr\n')
            np.savetxt(f,np.transpose(output), delimiter="\t") 

    def hor_age_save(self):
        output=np.vstack((self.LON, self.LAT, self.distance, self.hor_modage))
        header='#LON\tLAT\tdistance(km)'
        for i in range(self.nbhor):
            header=header+'\thor_no_'+str(i+1)
        header=header+'\n'
        with open(self.label+'agehorizons.txt','w') as f:
            f.write(header)
            np.savetxt(f,np.transpose(output), delimiter="\t") 
        for i in range(self.nbhor):
            print 'horizon no:',i+1,', average age: ',np.nanmean(self.hor_modage[i,:]),', stdev age: ',np.nanstd(self.hor_modage[i,:])

    def iso_age_save(self):
        output=np.vstack((self.LON, self.LAT, self.distance, self.iso_modage, self.iso_modage_sigma))
        header='#LON\tLAT\tdistance(km)'
        for i in range(self.nbiso):
            header=header+'\tiso_no_'+str(i+1)
        for i in range(self.nbiso):
            header=header+'\tiso_no_'+str(i+1)
        header=header+'\n'
        with open(self.label+'ageisochrones.txt','w') as f:
            f.write(header)
            np.savetxt(f,np.transpose(output), delimiter="\t") 
        for i in range(self.nbiso):
            print 'isochrone no:',i+1,', average age: ',np.nanmean(self.iso_modage[i,:]),', stdev age: ',np.nanstd(self.iso_modage[i,:])


    def twtt_save(self):
#        b=np.chararray((np.size(self.distance)), itemsize=20)
#        b[:]=RLlabel
#        print b
#        print self.LON
#        output=np.vstack((self.LON_twtt,self.LAT_twtt,self.twtt1Myr))
        current_folder_path, current_folder_name = os.path.split(RL.label[:-1])
#        print 'folder: ',RL.label,current_folder_name
        with open(self.label+'twtt.txt','w') as f:
            f.write('#LON\tLAT\ttwtt-0.6Myr(lk)\ttwtt-0.8Myr(lk)\ttwtt-1Myr(lk)\ttwtt-1.2Myr(lk)\ttwtt-1.5Myr(lk)\tself.twttBed(lk)\tLabel\n')
            for j in range(np.size(self.distance)):
                f.write('{0:11}'.format(str(self.LON_twtt[j]))+'\t'+'{0:12}'.format(str(self.LAT_twtt[j]))+'\t'+'{0:13}'.format(str(self.twtt0dot6Myr[j]))+'\t'+'{0:13}'.format(str(self.twtt0dot8Myr[j]))+'\t'+'{0:13}'.format(str(self.twtt1Myr[j]))+'\t'+'{0:13}'.format(str(self.twtt1dot2Myr[j]))+'\t'+'{0:13}'.format(str(self.twtt1dot5Myr[j]))+'\t'+'{0:13}'.format(str(self.twttBed[j]))+'\t'+current_folder_name+'\n')
#            np.savetxt(f,np.transpose(output), delimiter="\t") 


    def EDC(self):
        f=interp1d(self.distance,self.age)
        age_EDC=f(self.distance_EDC)
        g=interp1d(self.distance,self.sigma_age)
        sigmaage_EDC=g(self.distance_EDC)
        h=interp1d(self.distance,self.depth)
        depth_EDC=h(self.distance_EDC)
        print 'Thickness at EDC is:',depth_EDC[-1]
        i=interp1d(depth_EDC,age_EDC)
        age_EDC_bot=i(max(depth_EDC)-60)
        j=interp1d(depth_EDC,sigmaage_EDC)
        sigmaage_EDC_bot=j(max(depth_EDC)-60)
        print 'Age at EDC at 3200 m depth is: ',age_EDC_bot,'+-',sigmaage_EDC_bot
        f=interp1d(self.distance,self.p)
        p_EDC=f(self.distance_EDC)
        print 'p parameter at EDC is: ', p_EDC
        f=interp1d(self.distance,self.a)
        a_EDC=f(self.distance_EDC)
        print 'accumulation at EDC is: ', a_EDC
        f=interp1d(self.distance,self.m)
        m_EDC=f(self.distance_EDC)
        print 'melting at EDC is: ', m_EDC

        readarray=np.loadtxt(self.label+'temperatures_EDC.txt')
        datatempEDC_depth=readarray[:,0]
        datatempEDC_temp=readarray[:,1]+273.15
        f=interp1d(self.distance,self.T)
        temp_EDC=f(self.distance_EDC)
        f=interp1d(self.distance,self.T_anal)
        temp_anal_EDC=f(self.distance_EDC)
        plt.figure('Temperature at EDC')
        plt.plot(temp_EDC, depth_EDC, label='model')
        plt.plot(temp_anal_EDC, depth_EDC, label='analytical solution')
        plt.plot(datatempEDC_temp, datatempEDC_depth, label='data')
        plt.ylabel('depth (m)')
        plt.xlabel('temperature (K)')
        x1,x2,y1,y2 = plt.axis()
        plt.axis((x1,x2,y2,0.))
        plt.legend()
        pp=PdfPages(self.label+'temperatures_EDC.pdf')
        pp.savefig(plt.figure('Temperature at EDC'))
        pp.close()

#    def max_age(self):
#        j=np.argmax(self.age100)
#        self.agemax=max(self.agebot)
#        self.sigmaagemax=self.sigmaagebot[j]
#        self.distanceagemax=self.distance[j]
#        f=interp1d(self.distance_raw,self.LON_raw)
#        g=interp1d(self.distance_raw,self.LAT_raw)
#        self.LONagemax=f(self.distanceagemax)
#        self.LATagemax=g(self.distanceagemax)
#        print 'Maximum age is ',self.agemax,'+-',self.sigmaagemax
#        print 'It occurs at a distance of ',self.distanceagemax,' km, with coordinates ',self.LONagemax,', ',self.LATagemax        





#Main
RLlabel=sys.argv[1]
if RLlabel[-1]!='/':
    RLlabel=RLlabel+'/'
print 'Radar line is: ',RLlabel
print 'Creation of the radar line'
RL=RadarLine(RLlabel)
print 'Initialization of radar line'
RL.init()
print 'Data display'
RL.data_display()
if RL.opt_method=='leastsq':
    print 'Optimization by leastsq'
    RL.variables=np.concatenate((RL.a,RL.pprime))
    if RL.invert_G0:
        RL.variables=np.concatenate((RL.variables,RL.G0))
    if RL.invert_thk:
        RL.variables=np.concatenate((RL.variables,RL.thk))
#        self.variables=np.concatenate((self.a,self.m,self.s))
    RL.variables,RL.hess,infodict,mesg,ier=leastsq(RL.residuals, RL.variables, full_output=1)
    print mesg
    RL.residuals(RL.variables)
    print 'Calculation of confidence intervals'
    if RL.calc_sigma==False:
        RL.hess=np.zeros((np.size(RL.variables),np.size(RL.variables)))
    RL.sigma()
elif RL.opt_method=='none':
    RL.variables=np.concatenate((RL.a,RL.pprime))
    if RL.invert_G0:
        RL.variables=np.concatenate((RL.variables,RL.G0))
    if RL.invert_thk:
        RL.variables=np.concatenate((RL.variables,RL.thk))
    print 'No optimization'
    RL.residuals(RL.variables)
elif RL.opt_method=='none1D':
    print 'Forward model 1D'
    for j in range(np.size(RL.distance)):
        RL.variables1D=np.array([RL.a[j],RL.pprime[j]])
        if RL.invert_G0:
            RL.variables1D=np.append(RL.variables1D,RL.G0[j])
        if RL.invert_thk:
            RL.variables1D=np.append(RL.variables1D,RL.thk[j])
        RL.residuals1D(RL.variables1D,j)
elif RL.opt_method=='leastsq1D':
    print 'Optimization by leastsq1D'    
    for j in range(np.size(RL.distance)):
#        if RL.thk[j]<>np.nan and 
        print 'index along the radar line: ', j
        RL.variables1D=np.array([RL.a[j],RL.pprime[j]])
        if RL.invert_G0:
            RL.variables1D=np.append(RL.variables1D,RL.G0[j])
        if RL.invert_thk:
            RL.variables1D=np.append(RL.variables1D,RL.thk[j])
        RL.variables1D,RL.hess1D,infodict,mesg,ier=leastsq(RL.residuals1D, RL.variables1D, args=(j), full_output=1)
        RL.residuals1D(RL.variables1D,j)
        if RL.calc_sigma==False:
          RL.hess1D=np.zeros((np.size(RL.variables1D),np.size(RL.variables1D)))
        RL.sigma1D(j)
    RL.agebotmin=RL.agebot-RL.sigmabotage
elif RL.opt_method=='MH1D':
    print 'Optimization by MH1D'    
    for j in range(np.size(RL.distance)):
        print 'index along the radar line: ', j
        RL.variables1D=np.array([RL.a[j],RL.pprime[j]])
        if RL.invert_G0:
            RL.variables1D=np.append(RL.variables1D,RL.G0[j])
        if RL.invert_thk:
            RL.variables1D=np.append(RL.variables1D,RL.thk[j])
        RL.variables1D,RL.hess1D,infodict,mesg,ier=leastsq(RL.residuals1D, RL.variables1D, args=(j), full_output=1)
        step=RL.hess1D
        if RL.invert_G0 and RL.invert_thk and RL.variables1D[3]>RL.thkreal[j]:
            RL.variables1D[3]=RL.thkreal[j]
            print 'thk > threal in the leastsq solution'
#        step=np.diag(np.array([0.001, 0.1, 0.001, 10.])**2)
#        cost_accepted=np.array([])
#        variables1D_accepted=np.array([np.empty_like(RL.variables1D)])
#        agebot_accepted=np.array([])
#        agebot=RL.agebot[j]
#        accu_accepted=np.array([])
#        G0_accepted=np.array([])
#        melt_accepted=np.array([])
#        pprime_accepted=np.array([])
#        age_accepted=np.transpose(np.array([RL.age[:,j]]))
#        accu=RL.a[j]
#        G0=RL.G0[j]
#        melt=RL.m[j]
#        pprime=RL.pprime[j]
#        age=np.transpose(np.empty_like([RL.age[:,j]]))
        cost=RL.cost_fct(RL.variables1D,j)

        cost_accepted=np.array([cost])
        variables1D_accepted=np.array([RL.variables1D])
        agebot_accepted=np.array([RL.agebot[j]])
        accu_accepted=np.array([RL.a[j]])
        G0_accepted=np.array([RL.G0[j]])
        melt_accepted=np.array([RL.m[j]])
        pprime_accepted=np.array([RL.pprime[j]])
        age_accepted=np.transpose(np.array([RL.age[:,j]]))  #FIXME: This does not work since the depth grid varies!

        agebot=RL.agebot[j]
        accu=RL.a[j]
        G0=RL.G0[j]
        melt=RL.m[j]
        pprime=RL.pprime[j]
        age=np.transpose(np.array([RL.age[:,j]]))

        for iter in range(RL.MHnbiter):
#            print 'iteration no:',iter
            if iter==RL.MHiter_adapt1 or iter==RL.MHiter_adapt2:
                step=np.cov(np.transpose(variables1D_accepted))
                cost_accepted=np.array([cost])
                variables1D_accepted=np.array([RL.variables1D])
                agebot_accepted=np.array([RL.agebot[j]])
                accu_accepted=np.array([RL.a[j]])
                G0_accepted=np.array([RL.G0[j]])
                melt_accepted=np.array([RL.m[j]])
                pprime_accepted=np.array([RL.pprime[j]])
                age_accepted=np.transpose(np.array([RL.age[:,j]]))
            RL.variables1Dtest=np.random.multivariate_normal(RL.variables1D,step)
#            print RL.variables1Dtest[3],RL.invert_thk,RL.thkreal[j],RL.iso[-1,j]
            if (RL.invert_thk and RL.variables1Dtest[3]<=RL.thkreal[j] and RL.variables1Dtest[3]>max(RL.iso[:,j])) or not RL.invert_thk: #This is dirty, it does not work when we invert the thickness but not the GF!
                costtest=RL.cost_fct(RL.variables1Dtest,j)
#                print costtest
                if m.log(random.uniform(0,1))<=cost-costtest:
                    cost=costtest
                    RL.variables1D=RL.variables1Dtest
                    agebot=RL.agebot[j]
                    accu=RL.a[j]
                    G0=RL.G0[j]
                    melt=RL.m[j]
                    pprime=RL.pprime[j]
                    age=np.transpose(np.array([RL.age[:,j]]))
            
            variables1D_accepted=np.vstack((variables1D_accepted,RL.variables1D))
            cost_accepted=np.append(cost_accepted,cost)
            agebot_accepted=np.append(agebot_accepted,agebot)
            accu_accepted=np.append(accu_accepted,accu)
            G0_accepted=np.append(G0_accepted,G0)
            melt_accepted=np.append(melt_accepted,melt)
            pprime_accepted=np.append(pprime_accepted,pprime)
            age_accepted=np.hstack((age_accepted,age))
#            print RL.variables1D
        RL.agebotmin[j]=np.percentile(agebot_accepted,15)
        RL.sigma_a[j]=np.std(accu_accepted)
        RL.sigma_G0[j]=np.std(G0_accepted)
        RL.sigma_m[j]=np.std(melt_accepted)
        RL.sigma_pprime[j]=np.std(pprime_accepted)
        RL.sigma_age[:,j]=np.std(age_accepted, axis=1)
        print 'min age at 85%',RL.agebotmin[j]
        RL.variables1D=variables1D_accepted[np.argmin(cost_accepted),:]
#        print variables1D_accepted
#        print cost_accepted
#        RL.variables1D,RL.hess1D,infodict,mesg,ier=leastsq(RL.residuals1D, RL.variables1D, args=(j), full_output=1)
        RL.residuals1D(RL.variables1D,j)


else:
    print RL.opt_method,': Optimization method not recognized.'
    quit()
print 'calculating per layer accumulations'
#RL.model()

RL.accu_layers()





print 'Model display'
RL.model_display()
print 'parameters display'
RL.parameters_display()
RL.parameters_save()
if RL.is_EDC:
    RL.EDC()
#RL.max_age()
RL.bot_age_save()
RL.hor_age_save()
RL.iso_age_save()
RL.twtt_save()
print 'Program execution time: ', time.time() - start_time, 'seconds' 
