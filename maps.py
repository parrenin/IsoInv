from mpl_toolkits.basemap import Basemap
#from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LogNorm
from matplotlib.colors import Normalize
#from PIL import Image
from osgeo import gdalconst
#from scipy.interpolate import interp1d
import numpy as np
import matplotlib.pyplot as plt
import gdal
import sys
import yaml

#write_data=True

#lat1=-75.5
#lon1=124.
#lat2=-75.
#lon2=121.
lon1_bm2=-135
lat1_bm2=-48.458667
lon2_bm2=45
lat2_bm2=-48.458667
##

#Setting RadarLines directory
RLDir=sys.argv[1]
if RLDir[-1]!='/':
    RLDir=RLDir+'/'

#exec(open(RLDir+'parameters_maps.py').read())
run_model = False
output_format = "png"
lat1 = -75.35
lon1 = 126.7
lat2 = -75
lon2 = 120.5
dotsize = 4.
list_RL = []
list_RL_extra = []
list_RL_highlight = []
nbiso = 0
nbhor = 0
is_drill = True
lon_drill = 123.+21./60.
lat_drill = -75.1
lat_step = 0.25
lon_step = 1.
scale = False
labels = []

data = yaml.load(open(RLDir+'parameters_maps.yml').read(), Loader=yaml.FullLoader)
globals().update(data)

#Reading isochrones' ages
readarray=np.loadtxt(RLDir+'ages.txt')
iso_age=np.concatenate((np.array([0]),readarray[:,0]))

#Running model for each radar line
if run_model:
    for i,RLlabel in enumerate(list_RL+list_RL_extra):
        directory=RLDir+RLlabel
        sys.argv=['age_model2.py', directory]
        exec(open('age_model2.py').read())
        plt.close("all")

#Reading data for each radar line
for i,RLlabel in enumerate(list_RL):
    directory=RLDir+RLlabel
    accu_array1=np.loadtxt(directory+'/a.txt')
    botage_array1=np.loadtxt(directory+'/agebottom.txt')
    m_array1=np.loadtxt(directory+'/m.txt')
    G0_array1=np.loadtxt(directory+'/G0.txt')
    pprime_array1=np.loadtxt(directory+'/p.txt')
    hor_array1=np.loadtxt(directory+'/agehorizons.txt')
    if i==0:
        accu_array=accu_array1
        botage_array=botage_array1
        m_array=m_array1
        G0_array=G0_array1
        pprime_array=pprime_array1
        hor_array=hor_array1
    else:
        accu_array=np.concatenate((accu_array,accu_array1))
        botage_array=np.concatenate((botage_array,botage_array1))
        m_array=np.concatenate((m_array,m_array1))
        G0_array=np.concatenate((G0_array,G0_array1))
        pprime_array=np.concatenate((pprime_array,pprime_array1))
        hor_array=np.concatenate((hor_array,hor_array1))

for i,RLlabel in enumerate(list_RL_highlight):
    directory=RLDir+RLlabel
    botage_array1=np.loadtxt(directory+'/agebottom.txt')
    if i==0:
        botage_array_highlight=botage_array1
    else:
        botage_array_highlight=np.concatenate((botage_array,botage_array1))

#Reading data for extra radar lines
for  i,RLlabel in enumerate(list_RL_extra):
    directory=RLDir+RLlabel
    botage_array1=np.loadtxt(directory+'/agebottom.txt')
    m_array1=np.loadtxt(directory+'/m.txt')
    G0_array1=np.loadtxt(directory+'/G0.txt')
    pprime_array1=np.loadtxt(directory+'/pprime.txt')

    botage_array=np.concatenate((botage_array,botage_array1))
    m_array=np.concatenate((m_array,m_array1))
    G0_array=np.concatenate((G0_array,G0_array1))
    pprime_array=np.concatenate((pprime_array,pprime_array1))

#Reading AICC2012 file
#readarray=np.loadtxt(RLDir+'/AICC2012.txt')
#AICC2012_depth=readarray[:,0]
#AICC2012_iedepth=readarray[:,1]
#AICC2012_accu=readarray[:,2]
#AICC2012_age=readarray[:,3]
#AICC2012_sigma=readarray[:,4]

#AICC2012_averageaccu=np.sum((AICC2012_age[1:]-AICC2012_age[:-1])*AICC2012_accu[:-1])/(AICC2012_age[-1]-AICC2012_age[0])
#AICC2012_steadyage=np.cumsum(np.concatenate((np.array([AICC2012_age[0]]),(AICC2012_age[1:]-AICC2012_age[:-1])*AICC2012_accu[:-1]/AICC2012_averageaccu)))


#
list_maps=['resolution-1.5Myr','radar-lines','bottom-age','accu-steady','pprime', 'geothermal-heat-flux','melting','melting-sigma','Height-Above-Bed-0.8Myr','Height-Above-Bed-1Myr','Height-Above-Bed-1.2Myr','Height-Above-Bed-1.5Myr','min-bottom-age','age-100m','age-150m','age-200m','age-250m', 'resolution-1Myr','resolution-1.2Myr','geothermal-heat-flux-sigma','pprime-sigma','accu-sigma']
list_length=len(list_maps)
for i in range(nbiso):
    list_maps.append('accu-layer'+ "%02i"%(i+1) +'_'+str(int(iso_age[i]/1000.))+'-'+str(int(iso_age[i+1]/1000.))+'kyr' )
for i in range(nbhor):
    list_maps.append('age-hor'+"%02i"%(i+1))

for i,MapLabel in enumerate(list_maps):

    print(MapLabel+' map')

    fig=plt.figure(MapLabel,figsize=(21/2.54,21/2.54)) 
    plt.title(MapLabel, y=1.05)
    map0 = Basemap(projection='stere', lat_ts=-71, lat_0=-90, lon_0=180, llcrnrlon=lon1_bm2,llcrnrlat=lat1_bm2, urcrnrlon=lon2_bm2,urcrnrlat=lat2_bm2, rsphere=(6378137.00,6356752.3142))
    map1 = Basemap(projection='stere', lat_ts=-71, lat_0=-90, lon_0=180, llcrnrlat=lat1, llcrnrlon=lon1, urcrnrlat=lat2, urcrnrlon=lon2, rsphere=(6378137.00,6356752.3142))

    map1.drawparallels(np.arange(-90.,81.,lat_step), labels=[True, False, False, True], dashes=[1, 5], color='0.5')
    map1.drawmeridians(np.arange(-180.,180.,lon_step), latmax=85., labels=[False, True, True, False], dashes=[1, 5], color='0.5')
    if len(scale) == 3:
        map1.drawmapscale(scale[0], scale[1], scale[0], scale[1], scale[2], yoffset=10., barstyle='simple')


    ##Draw bed topography
    #raster = gdal.Open('bedmap2/bedmap2_bed.txt')
    #band = raster.GetRasterBand(1)
    #array = band.ReadAsArray()
    #array=np.where(array==-9999,np.nan,array)

    #map1.imshow(array[::-1,:])
    #map1.colorbar()


    ##Draw surface contours
    raster2 = gdal.Open(RLDir+'../bedmap2/bedmap2_surface.txt')
    band2 = raster2.GetRasterBand(1)
    array2 = band2.ReadAsArray()
    array2=np.where(array2==-9999,np.nan,array2)
    zz=array2[::-1,:]

    x = np.linspace(0, map0.urcrnrx, array2.shape[1])
    y = np.linspace(0, map0.urcrnry, array2.shape[0])
    x1,y1=map0(lon1,lat1)
    x2,y2=map0(lon2,lat2)
    x=x-x1
    y=y-y1
    xx, yy = np.meshgrid(x, y)


    if MapLabel[:4] != 'accu':
        levels=np.concatenate(( np.arange(3150, 3260, 10),np.arange(3260,3270, 2) ))
    else:
        levels=np.concatenate(( np.arange(3150, 3260, 2),np.arange(3260,3270, 2) ))
    cs=map1.contour(xx,yy, zz, colors='0.5', levels=levels, alpha=0.4)
    plt.clabel(cs, inline=1, fontsize=10,fmt='%1.0f')

    ##Draw bedrock contours

    raster2 = gdal.Open(RLDir+'../bedmap2/bedmap2_bed.txt')
    band2 = raster2.GetRasterBand(1)
    array2 = band2.ReadAsArray()
    array2=np.where(array2==-9999,np.nan,array2)
    zz=array2[::-1,:]

    x = np.linspace(0, map0.urcrnrx, array2.shape[1])
    y = np.linspace(0, map0.urcrnry, array2.shape[0])
    x1,y1=map0(lon1,lat1)
    x2,y2=map0(lon2,lat2)
    x=x-x1
    y=y-y1
    xx, yy = np.meshgrid(x, y)


    levels=np.arange(-1000., 900., 100.)

    if MapLabel=='radar-lines':
        cmap='terrain'
        alpha=0.4
    else:
        cmap='gray'
        alpha=0.6
    plt.imshow(zz[::-1,:], extent=[max(x),min(x),max(y),min(y)], cmap=cmap, norm=Normalize(vmin=-700, vmax=600), alpha=alpha)

    ##Draw OIA refined bedrock
#    img = Image.open(RLDir+'bedmap2/Bed_BlobA.tiff')
#    arr = np.asarray(img)
#    arr=np.where(arr==-9999,np.nan,arr)
    
    def readRasterBandAsArray(filename, bandnum):
        raster = gdal.Open(filename, gdalconst.GA_ReadOnly)
        rasterBand = raster.GetRasterBand(bandnum)
        rasterBandArray = rasterBand.ReadAsArray(0, 0, raster.RasterXSize, raster.RasterYSize).astype(np.float)
         
        rasterBandNoDataValue = rasterBand.GetNoDataValue()
        if rasterBandNoDataValue is not None:
            rasterBandArray[rasterBandArray == rasterBandNoDataValue] = np.nan
                             
        return rasterBandArray
 
    rasterBandArray=readRasterBandAsArray(RLDir+'../bedmap2/Bed_BlobA_Geoid4.tif',1)

    latmax=-75.1164861
    latmin=-75.5905194
    lonmax=121.1456639
    lonmin=124.3964778
    hmin,vmin=map1(lonmin,latmin)
    hmax,vmax=map1(lonmax,latmax)
    extent=(hmin, hmax, vmin, vmax)
    xborders=np.array([hmin,hmax,hmax,hmin,hmin])
    yborders=np.array([vmin,vmin,vmax,vmax,vmin])

    ##Draw color bar
    if MapLabel=='radar-lines':
        cb0=plt.colorbar(orientation='horizontal', shrink=0.7, pad=0.05)
    else:
        cb0=plt.colorbar(orientation='horizontal', shrink=0.7, pad=0)
    cb0.set_label('Bedrock elevation (m)')

    #Draw continent's contour
    #raster3 = gdal.Open('bedmap2/bedmap2_icemask_grounded_and_shelves.txt')
    #band3 = raster3.GetRasterBand(1)
    #array3 = band3.ReadAsArray()

    #x = np.linspace(0, map1.urcrnrx, array3.shape[1])
    #y = np.linspace(0, map1.urcrnry, array3.shape[0])
    #xx, yy = np.meshgrid(x, y)
    #map1.contour(xx,yy, array3[::-1,:], colors='k')

    levels='auto'


    if MapLabel=='radar-lines':

        LON=botage_array[:,0]
        LAT=botage_array[:,1]
        x,y=map1(LON,LAT)

        map1.scatter(x,y, c='b', marker='o', lw=0., edgecolor='', s=dotsize)

        LON=botage_array_highlight[:,0]
        LAT=botage_array_highlight[:,1]
        x,y=map1(LON,LAT)

        map1.scatter(x,y, c='r', marker='o', lw=0., edgecolor='', s=dotsize)


#        ax2 = plt.axes()
#        lon=124.4520
#        lat=-74.8901
#        x,y=map1(lon,lat)
#        plt.text(x,y,'A', color='red', fontweight='bold')
#
#        lon=121.7368
#        lat=-75.4481
#        x,y=map1(lon,lat)
#        plt.text(x,y,"A'", color='red', fontweight='bold')
#        
#        lon=124.7
#        lat=-75.1
#        x,y=map1(lon,lat)
#        plt.text(x,y,'Concordia Ridge',horizontalalignment='center',verticalalignment='center',rotation=-32)
#        lon=124.2
#        lat=-75.2
#        x,y=map1(lon,lat)
#        plt.text(x,y,'Concordia Subglacial Trench',horizontalalignment='center',verticalalignment='center',rotation=-32)
#
#        candidates=readRasterBandAsArray(RLDir+'../bedmap2/candidates_Brice_clipped.tif',1)
#        candidates1=np.where(candidates != np.nan,candidates*0+10,np.nan)
#        latmax=-75.263709
#        latmin=-73.750157
#        lonmax=132.986033
#        lonmin=115.231155
#        hmin,vmin=map1(lonmin,latmin)
#        hmax,vmax=map1(lonmax,latmax)
#        extent=(hmax, hmin, vmin, vmax)
#        plt.imshow(candidates1, origin='lower', cmap='Oranges_r', extent=extent, alpha=0.2, interpolation='none')
#
#        lon=122
#        lat=-75.1
#        x,y=map1(lon,lat)
#        plt.text(x,y,'A',horizontalalignment='center',verticalalignment='center',color='orange',bbox=dict(facecolor='white', edgecolor='#767876'))
#        lon=124.9
#        lat=-75.2
#        x,y=map1(lon,lat)
#        plt.text(x,y,'B',horizontalalignment='center',verticalalignment='center',color='orange',bbox=dict(facecolor='white', edgecolor='#767876'))
#        lon=124.95
#        lat=-75
#        x,y=map1(lon,lat)
#        plt.text(x,y,'C',horizontalalignment='center',verticalalignment='center',color='orange',bbox=dict(facecolor='white', edgecolor='#767876'))
#        lon=124.95
#        lat=-74.85
#        x,y=map1(lon,lat)
#        plt.text(x,y,'D',horizontalalignment='center',verticalalignment='center',color='orange',bbox=dict(facecolor='white', edgecolor='#767876'))
#
#        lon=121.
#        lat=-75.05
#        x,y=map1(lon,lat)
#        plt.text(x,y,'E',horizontalalignment='center',verticalalignment='center',color='orange',bbox=dict(facecolor='white', edgecolor='#767876'))


    if MapLabel=='bottom-age':

        LON=botage_array[:,0]
        LAT=botage_array[:,1]
        botage=botage_array[:,4]
        x,y=map1(LON,LAT)

        norm = LogNorm(vmin=0.7,vmax=5.)
        map1.scatter(x,y, c=botage/1e6, marker='o', lw=0., edgecolor='', norm = norm, s=dotsize)
        cblabel='Bottom age (Myr)'
        levels=np.array([0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 2, 2.5, 3, 3.5, 4, 5])

        ax2 = plt.axes()
        lon=122.4
        lat=-75.32
        x,y=map1(lon,lat)
        plt.text(x,y,'LDCP',horizontalalignment='center',verticalalignment='center')
        lon=124.02
        lat=-74.985
        x,y=map1(lon,lat)
        plt.text(x,y,'NP',horizontalalignment='center',verticalalignment='center')


#        output=np.transpose(np.vstack((LON,LAT,botage)))
#        with open(RLDir+'agebottom.txt','w') as f:
#            f.write('#LON\tLAT\tbottom age (yr)\n')
#            np.savetxt(f,output, delimiter="\t")

    if MapLabel=='min-bottom-age':

        LON=botage_array[:,0]
        LAT=botage_array[:,1]
        minbotage=botage_array[:,5]
        x,y=map1(LON,LAT)

        norm = LogNorm(vmin=0.7,vmax=5.)
        map1.scatter(x,y, c=minbotage/1e6, marker='o', lw=0., edgecolor='', norm = norm, s=dotsize)
        cblabel='Minimum bottom age (Myr)'
        levels=np.array([0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 2, 2.5, 3, 3.5, 4, 5])

        ax2 = plt.axes()
        lon=122.4
        lat=-75.32
        x,y=map1(lon,lat)
        plt.text(x,y,'LDCP',horizontalalignment='center',verticalalignment='center')
        lon=124.02
        lat=-74.985
        x,y=map1(lon,lat)
        plt.text(x,y,'NP',horizontalalignment='center',verticalalignment='center')


#        output=np.transpose(np.vstack((LON,LAT,minbotage)))
#        with open(RLDir+'minagebottom.txt','w') as f:
#            f.write('#LON\tLAT\tmin bottom age (yr)\n')
#            np.savetxt(f,output, delimiter="\t")

    if MapLabel=='age-100m':

        LON=botage_array[:,0]
        LAT=botage_array[:,1]
        botage=botage_array[:,6]
        x,y=map1(LON,LAT)

        norm = LogNorm(vmin=0.7,vmax=5.)
        map1.scatter(x,y, c=botage/1e6, marker='o', lw=0., edgecolor='', norm = norm, s=dotsize)
        cblabel='Age (Myr)'
        levels=np.array([0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 2, 2.5, 3, 3.5, 4, 5])

    if MapLabel=='age-150m':

        LON=botage_array[:,0]
        LAT=botage_array[:,1]
        botage=botage_array[:,7]
        x,y=map1(LON,LAT)

        norm = LogNorm(vmin=0.7,vmax=5.)
        map1.scatter(x,y, c=botage/1e6, marker='o', lw=0., edgecolor='', norm = norm, s=dotsize)
        cblabel='Age (Myr)'
        levels=np.array([0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 2, 2.5, 3, 3.5, 4, 5])

    if MapLabel=='age-200m':

        LON=botage_array[:,0]
        LAT=botage_array[:,1]
        botage=botage_array[:,8]
        x,y=map1(LON,LAT)

        norm = LogNorm(vmin=0.7,vmax=5.)
        map1.scatter(x,y, c=botage/1e6, marker='o', lw=0., edgecolor='', norm = norm, s=dotsize)
        cblabel='Age (Myr)'
        levels=np.array([0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 2, 2.5, 3, 3.5, 4, 5])

    if MapLabel=='age-250m':

        LON=botage_array[:,0]
        LAT=botage_array[:,1]
        botage=botage_array[:,9]
        x,y=map1(LON,LAT)

        norm = LogNorm(vmin=0.7,vmax=5.)
        map1.scatter(x,y, c=botage/1e6, marker='o', lw=0., edgecolor='', norm = norm, s=dotsize)
        cblabel='Age (Myr)'
        levels=np.array([0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 2, 2.5, 3, 3.5, 4, 5])

    if MapLabel=='resolution-1Myr':

        LON=botage_array[:,0]
        LAT=botage_array[:,1]
        resolution=botage_array[:,10]
        x,y=map1(LON,LAT)

        norm = Normalize()
        map1.scatter(x,y, c=resolution/1e3, marker='o', lw=0., edgecolor='', norm = norm, s=dotsize)
        cblabel='Resolution at 1Myr (kyr m$^{-1}$)'
#        levels=np.array([1., 2., 4., 6., 8., 10., 20., 40.])

#        output=np.transpose(np.vstack((LON,LAT,resolution/1e3)))
#        with open(RLDir+'resolution1Myr.txt','w') as f:
#            f.write('#LON\tLAT\tresolution (kyr/m)\n')
#            np.savetxt(f,output, delimiter="\t")

    if MapLabel=='resolution-1.2Myr':

        LON=botage_array[:,0]
        LAT=botage_array[:,1]
        resolution=botage_array[:,11]
        x,y=map1(LON,LAT)

        norm = Normalize()
        map1.scatter(x,y, c=resolution/1e3, marker='o', lw=0., edgecolor='', norm = norm, s=dotsize)
        cblabel='Resolution at 1.2Myr (kyr m$^{-1}$)'
#        levels=np.array([1., 2., 4., 6., 8., 10., 20., 40.])

#        output=np.transpose(np.vstack((LON,LAT,resolution/1e3)))
#        with open(RLDir+'resolution1.2Myr.txt','w') as f:
#            f.write('#LON\tLAT\tresolution (kyr/m)\n')
#            np.savetxt(f,output, delimiter="\t")

    if MapLabel=='resolution-1.5Myr':

        LON=botage_array[:,0]
        LAT=botage_array[:,1]
        resolution=botage_array[:,12]
        x,y=map1(LON,LAT)

#        norm = Normalize(vmax=20.)
        norm = Normalize()
        map1.scatter(x,y, c=resolution/1e3, marker='o', lw=0., edgecolor='', norm = norm, s=dotsize)
        cblabel='Resolution at 1.5Myr (kyr m$^{-1}$)'
#        levels=np.array([1., 2., 4., 6., 8., 10., 20., 40.])

        ax2 = plt.axes()
#        lon=122.4
#        lat=-75.32
#        x,y=map1(lon,lat)
#        plt.text(x,y,'LDCP',horizontalalignment='center',verticalalignment='center')
#        lon=124.02
#        lat=-74.985
#        x,y=map1(lon,lat)
#        plt.text(x,y,'NP',horizontalalignment='center',verticalalignment='center')


#        output=np.transpose(np.vstack((LON,LAT,resolution/1e3)))
#        with open(RLDir+'resolution1.5Myr.txt','w') as f:
#            f.write('#LON\tLAT\tresolution (kyr/m)\n')
#            np.savetxt(f,output, delimiter="\t")

    if MapLabel=='Height-Above-Bed-0.8Myr':

        LON=botage_array[:,0]
        LAT=botage_array[:,1]
        height=botage_array[:,14]
        x,y=map1(LON,LAT)

        res=map1.scatter(x,y, c=height, marker='o', lw=0., edgecolor='', s=dotsize)
        cblabel='Height above bed (m)'

    if MapLabel=='Height-Above-Bed-1Myr':

        LON=botage_array[:,0]
        LAT=botage_array[:,1]
        height=botage_array[:,15]
        x,y=map1(LON,LAT)

        res=map1.scatter(x,y, c=height, marker='o', lw=0., edgecolor='', s=dotsize)
        cblabel='Height above bed (m)'
        

    if MapLabel=='Height-Above-Bed-1.2Myr':

        LON=botage_array[:,0]
        LAT=botage_array[:,1]
        height=botage_array[:,16]
        x,y=map1(LON,LAT)

        res=map1.scatter(x,y, c=height, marker='o', lw=0., edgecolor='', s=dotsize)
        cblabel='Height above bed (m)'
        

    if MapLabel=='Height-Above-Bed-1.5Myr':

        LON=botage_array[:,0]
        LAT=botage_array[:,1]
        height=botage_array[:,17]
        x,y=map1(LON,LAT)

        res=map1.scatter(x,y, c=height, marker='o', lw=0., edgecolor='', s=dotsize)
        cblabel='Height above bed (m)'
        
        ax2 = plt.axes()
        lon=122.4
        lat=-75.32
        x,y=map1(lon,lat)
        plt.text(x,y,'LDCP',horizontalalignment='center',verticalalignment='center')
        lon=124.02
        lat=-74.985
        x,y=map1(lon,lat)
        plt.text(x,y,'NP',horizontalalignment='center',verticalalignment='center')


#        levels=np.array([1., 2., 4., 6., 8., 10., 20., 40.])
#        cb.set_ticks(levels)
#        cb.set_ticklabels(levels)



    elif MapLabel=='melting':

        LON=m_array[:,0]
        LAT=m_array[:,1]
        melting=m_array[:,3]
        x,y=map1(LON,LAT)

        norm = Normalize(vmin=0.,vmax=3.)

        map1.scatter(x,y, c=melting*1e3, marker='o', lw=0., edgecolor='', s=dotsize, norm=norm)
        cblabel='$\overline{m}$ (mm yr$^{-1}$)'
        

#        output=np.transpose(np.vstack((LON,LAT,melting*1e3)))
#        with open(RLDir+'m.txt','w') as f:
#            f.write('#LON\tLAT\tmelting (mm/yr)\n')
#            np.savetxt(f,output, delimiter="\t")

    elif MapLabel=='melting-sigma':

        LON=m_array[:,0]
        LAT=m_array[:,1]
        sigma_melting=m_array[:,4]
        x,y=map1(LON,LAT)

        norm = Normalize(vmin=0.,vmax=1.)


        map1.scatter(x,y, c=sigma_melting*1e3, marker='o', lw=0., edgecolor='', s=dotsize, norm=norm)
        cblabel='$\sigma$ Melting (mm yr$^{-1}$)'
        

    elif MapLabel=='geothermal-heat-flux':

        LON=G0_array[:,0]
        LAT=G0_array[:,1]
        G0=G0_array[:,3]
        x,y=map1(LON,LAT)
        norm = Normalize(vmin=40.,vmax=100.)

        map1.scatter(x,y, c=G0*1e3, marker='o', lw=0., edgecolor='', s=dotsize, norm=norm)
        cblabel='$G_0$ (mW m$^{-2}$)'
        

    elif MapLabel=='geothermal-heat-flux-sigma':

        LON=G0_array[:,0]
        LAT=G0_array[:,1]
        sigma_G0=G0_array[:,4]
        x,y=map1(LON,LAT)

        map1.scatter(x,y, c=sigma_G0*1e3, marker='o', lw=0., edgecolor='', s=dotsize)
        cblabel='$\sigma_{G0}$ (mW m$^{-2}$)'
        

    elif MapLabel=='pprime':

        LON=pprime_array[:,0]
        LAT=pprime_array[:,1]
        pprime=pprime_array[:,3]
        x,y=map1(LON,LAT)

#        levels=np.arange(-1,3.1, 0.1)
        norm = Normalize(vmin=-1.,vmax=3.)
        map1.scatter(x,y, c=pprime, marker='o', lw=0., edgecolor='', s=dotsize, norm=norm)
        cblabel='$p$\''
        

#        cb.set_ticks(levels)

#        output=np.transpose(np.vstack((LON,LAT,pprime)))
#        with open(RLDir+'p.txt','w') as f:
#            f.write('#LON\tLAT\tpprime\n')
#            np.savetxt(f,output, delimiter="\t")


    elif MapLabel=='pprime-sigma':

        LON=pprime_array[:,0]
        LAT=pprime_array[:,1]
        sigma_pprime=pprime_array[:,4]
        x,y=map1(LON,LAT)

        norm = Normalize(vmin=0.,vmax=1.)
        map1.scatter(x,y, c=sigma_pprime, marker='o', lw=0., edgecolor='', s=dotsize, norm=norm)
        cblabel='$\sigma$ pprime'
        

    elif MapLabel[0:4]=='accu':

        LON=accu_array[:,0]
        LAT=accu_array[:,1]
        x,y=map1(LON,LAT)
        if MapLabel=='accu-steady':
            norm = Normalize(vmin=10.,vmax=25.)
            accu=accu_array[:,3]
        elif MapLabel=='accu-sigma':
            accu=accu_array[:,4]
            norm = Normalize(vmin=0.,vmax=1.)
#            output=np.transpose(np.vstack((LON,LAT,accu*100)))
#            with open(RLDir+'a.txt','w') as f:
#                f.write('#LON\tLAT\taccu (cm/yr)\n')
#                np.savetxt(f,output, delimiter="\t")

        else:
            norm=Normalize()
            i=int(MapLabel[10:12])
            accu=accu_array[:,i+4]
#            f=interp1d(AICC2012_age,AICC2012_steadyage)
#            accu=accu*(f(iso_age[i])-f(iso_age[i-1]))/(iso_age[i]-iso_age[i-1])

        map1.scatter(x,y, c=accu*1000*0.917, marker='o', lw=0., edgecolor='', s=dotsize, norm=norm)
        cblabel='$\overline{a}$ (mm-we yr$^{-1}$)'
        
    elif MapLabel[0:7]=='age-hor':
        LON=hor_array[:,0]
        LAT=hor_array[:,1]
        x,y=map1(LON,LAT)

        age=hor_array[:,int(MapLabel[7:9])+2]

        if np.all(np.isnan(age)):
            norm=Normalize(vmin=0.,vmax=1.)
        else:
            norm=Normalize(vmin=np.nanmin(age/1000.),vmax=np.nanmax(age/1000.))

        map1.scatter(x,y, c=age/1000., marker='o', lw=0., edgecolor='', s=dotsize, norm=norm)
        cblabel='age (kyr B1950)'


    if MapLabel != 'radar-lines':
        cb=plt.colorbar(orientation='horizontal', shrink=0.7, pad=0.05)
        cb.set_label(cblabel)
        if levels != 'auto':
            cb.set_ticks(levels)
            cb.set_ticklabels(levels)


    if is_drill:
        xdrill,ydrill=map1(lon_drill,lat_drill)
        map1.scatter(xdrill,ydrill, marker='*', c='r', edgecolor='r', s=64)

    ax2 = plt.axes()
    for i in range(len(labels)):
        x,y=map1(labels[i][0],labels[i][1])
        plt.text(x,y+10,labels[i][2],horizontalalignment='center',verticalalignment='bottom',
                 color='grey')
        map1.scatter(x,y, marker='o', s=dotsize, color='k')


    plt.tight_layout()

#    pp=PdfPages(RLDir+MapLabel+'.pdf')
#    pp.savefig(plt.figure(MapLabel))
#    pp.close()
    plt.savefig(RLDir+MapLabel+'.'+output_format, format=output_format, bbox_inches='tight')
    plt.close(fig)

plt.show()
