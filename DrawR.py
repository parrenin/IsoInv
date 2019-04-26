import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

RLlabel=sys.argv[1]
if RLlabel[-1]!='/':
    RLlabel=RLlabel+'/'
print('Radar line is: ',RLlabel)

#Reading the AICC2012 dataset, calculation of steady age and interpolation
readarray=np.loadtxt(RLlabel+'/AICC2012.txt')
AICC2012_depth=readarray[:,0]
AICC2012_iedepth=readarray[:,1]
AICC2012_accu=readarray[:,2]
AICC2012_age=readarray[:,3]
AICC2012_sigma=readarray[:,4]

AICC2012_averageaccu=np.sum((AICC2012_age[1:]-AICC2012_age[:-1])*AICC2012_accu[:-1])/(AICC2012_age[-1]-AICC2012_age[0])

R=AICC2012_accu/AICC2012_averageaccu

plt.step(np.concatenate((AICC2012_age/1000,np.array([1000]))),np.concatenate((R,np.array([1]))))
x1,x2,y1,y2 = plt.axis()
plt.axis((-0.050,1000,y1,y2))
plt.xlabel('age (kyr)')
plt.ylabel('R (no unit)')

pp=PdfPages(RLlabel+'R-vs-t.pdf')
pp.savefig()
pp.close()

