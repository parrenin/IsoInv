from mpl_toolkits.basemap import Basemap
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LogNorm
import numpy as np
import matplotlib.pyplot as plt
import gdal
import sys

run_model=False
lat1=-75.5
lon1=128.
lat2=-74.8
lon2=118.
lonEDC=123.+21./60.
latEDC=-75.1
dotsize=4.


list_RL=['icp7.jkb2n.edmc02a','icp7.jkb2n.f16t04a','icp7.jkb2n.ridge1a','mcm.jkb1a.edmc01a','oia.jkb2n.x39a','oia.jkb2n.x45a',\
'oia.jkb2n.x48a','oia.jkb2n.x54a','oia.jkb2n.x57a','oia.jkb2n.x60a','oia.jkb2n.x63a','oia.jkb2n.x66a','oia.jkb2n.x69a','oia.jkb2n.x72a',\
'oia.jkb2n.y15a','oia.jkb2n.y52a','oia.jkb2n.y60a','oia.jkb2n.y64a','oia.jkb2n.y68a','oia.jkb2n.y72a','oia.jkb2n.y74a','oia.jkb2n.y75a',\
'oia.jkb2n.y76a','oia.jkb2n.y77a','oia.jkb2n.y78a','oia.jkb2n.y79a','oia.jkb2n.y81a','oia.jkb2n.y82a','oia.jkb2n.y84a','oia.jkb2n.y86a','oia.jkb2n.y88a','oia.jkb2n.y90a','vcd.jkb2g.dvd01a']
for i,RLlabel in enumerate(list_RL):
    directory='explore.layers.bycolumn.'+RLlabel+'.llxydzn'
    if run_model:
        sys.argv=['AgeModel.py',directory]
        execfile('AgeModel.py')
        plt.close("all")
    accu_array1=np.loadtxt(directory+'/a.txt')
    botage_array1=np.loadtxt(directory+'/agebottom.txt')
    m_array1=np.loadtxt(directory+'/m.txt')
    G0_array1=np.loadtxt(directory+'/G0.txt')
    pprime_array1=np.loadtxt(directory+'/pprime.txt')
    if i==0:
        accu_array=accu_array1
        botage_array=botage_array1
        m_array=m_array1
        G0_array=G0_array1
        pprime_array=pprime_array1
    else:
        accu_array=np.concatenate((accu_array,accu_array1))
        botage_array=np.concatenate((botage_array,botage_array1))
        m_array=np.concatenate((m_array,m_array1))
        G0_array=np.concatenate((G0_array,G0_array1))
        pprime_array=np.concatenate((pprime_array,pprime_array1))



list_maps=['bottom-age','sigma-bottom-age','min-bottom-age','melting', 'geothermal-heat-flux','pprime','steady-accu']
for i in range(17):
    list_maps.append('accu-layer'+str(i+1))


for i,MapLabel in enumerate(list_maps):

    print MapLabel+' map'

    fig=plt.figure(MapLabel)
    plt.title(MapLabel, y=1.05)
    map0 = Basemap(projection='spstere', lat_ts=-71, boundinglat=-60, lon_0=180)
    map1 = Basemap(projection='stere', lat_ts=-71, lat_0=-90, lon_0=180, llcrnrlat=lat1, llcrnrlon=lon1, urcrnrlat=lat2, urcrnrlon=lon2)
    #map1 = Basemap(projection='spstere', boundinglat=-60, lon_0=180, llcrnrx=-4.5e6, llcrnry=-2.3e6, urcrnrx=-5e6, urcrnry=-2.8e6)
    #m = Basemap(projection='stere', lat_0=-75, lon_0=123., width=1e6, height=1e6)
    #m.drawcoastlines()
    #m.fillcontinents(color='white',lake_color='aqua')
    #m.drawmapboundary(fill_color='aqua')

    map1.drawparallels(np.arange(-90.,81.,1.), labels=[True, False, False, True], dashes=[1, 5], color='0.5')
    map1.drawmeridians(np.arange(-180.,180.,2.), latmax=85., labels=[False, True, True, False], dashes=[1, 5], color='0.5')
    map1.drawmapscale(lon1-1.2, lat1+0.2, lon1, lat1, 50)



    ##Draw bed topography
    #raster = gdal.Open('bedmap2/bedmap2_bed.txt')
    #band = raster.GetRasterBand(1)
    #array = band.ReadAsArray()
    #array=np.where(array==-9999,np.nan,array)

    #map1.imshow(array[::-1,:])
    #map1.colorbar()


    ##Draw surface contours
    raster2 = gdal.Open('bedmap2/bedmap2_surface.txt')
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
    #i1=np.min(np.where(x>=0))
    #i2=np.max(np.where(x<=x2-x1))
    #j1=np.min(np.where(y>=0))
    #j2=np.max(np.where(y<=y2-y1))
    #xxp=xx[i1:i2+1,j1:j2+1]
    #yyp=yy[i1:i2+1,j1:j2+1]
    #zzp=zz[i1:i2+1,j1:j2+1]


    levels=np.concatenate(( np.arange(3150, 3260, 10),np.arange(3260,3270, 2) ))
    cs=map1.contour(xx,yy, zz, colors='0.75', levels=levels, alpha=0.5)
    plt.clabel(cs, inline=1, fontsize=10,fmt='%1.0f')

    #Draw continent's contour
    #raster3 = gdal.Open('bedmap2/bedmap2_icemask_grounded_and_shelves.txt')
    #band3 = raster3.GetRasterBand(1)
    #array3 = band3.ReadAsArray()

    #x = np.linspace(0, map1.urcrnrx, array3.shape[1])
    #y = np.linspace(0, map1.urcrnry, array3.shape[0])
    #xx, yy = np.meshgrid(x, y)
    #map1.contour(xx,yy, array3[::-1,:], colors='k')

    if MapLabel=='bottom-age':

        LON=botage_array[:,0]
        LAT=botage_array[:,1]
        botage=botage_array[:,3]
        x,y=map1(LON,LAT)

        map1.scatter(x,y, c=botage/1e6, marker='o', lw=0., edgecolor='', norm = LogNorm(), s=dotsize)
        cb=map1.colorbar(pad='12%')
        cb.set_label('Bottom age (Myr)')
        levels=np.array([0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 2, 2.5, 3, 3.5, 4, 5])
        cb.set_ticks(levels)
        cb.set_ticklabels(levels)

    if MapLabel=='sigma-bottom-age':

        LON=botage_array[:,0]
        LAT=botage_array[:,1]
        sigmabotage=botage_array[:,4]
        x,y=map1(LON,LAT)

        map1.scatter(x,y, c=sigmabotage/1e6, marker='o', lw=0., edgecolor='', norm = LogNorm(), s=dotsize)
        cb=map1.colorbar(pad='12%')
        cb.set_label('Sigma bottom age (Myr)')
#        levels=np.array([0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 2, 2.5, 3, 3.5, 4, 5])
#        cb.set_ticks(levels)
#        cb.set_ticklabels(levels)

    if MapLabel=='min-bottom-age':

        LON=botage_array[:,0]
        LAT=botage_array[:,1]
        minbotage=botage_array[:,5]
        x,y=map1(LON,LAT)

        map1.scatter(x,y, c=minbotage/1e6, marker='o', lw=0., edgecolor='', s=dotsize)
        cb=map1.colorbar(pad='12%')
        cb.set_label('Minimum bottom age (Myr)')
#        levels=np.array([0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 2, 2.5, 3, 3.5, 4, 5])
#        cb.set_ticks(levels)
#        cb.set_ticklabels(levels)

    elif MapLabel=='melting':

        LON=m_array[:,0]
        LAT=m_array[:,1]
        melting=m_array[:,3]
        x,y=map1(LON,LAT)

        map1.scatter(x,y, c=melting*1e3, marker='o', lw=0., edgecolor='', s=dotsize)
        cb=map1.colorbar(pad='12%')
        cb.set_label('Melting (mm/yr)')

    elif MapLabel=='geothermal-heat-flux':

        LON=G0_array[:,0]
        LAT=G0_array[:,1]
        G0=G0_array[:,3]
        x,y=map1(LON,LAT)

        map1.scatter(x,y, c=G0*1e3, marker='o', lw=0., edgecolor='')
        cb=map1.colorbar(pad='12%')
        cb.set_label('G0 (mW/m$^2$)')

    elif MapLabel=='pprime':

        LON=pprime_array[:,0]
        LAT=pprime_array[:,1]
        pprime=pprime_array[:,3]
        x,y=map1(LON,LAT)

        map1.scatter(x,y, c=pprime, marker='o', lw=0., edgecolor='', s=dotsize)
        cb=map1.colorbar(pad='12%')
        cb.set_label('pprime')

    elif i>=6:

        LON=accu_array[:,0]
        LAT=accu_array[:,1]
        if i==6:
            accu=accu_array[:,3]
        else:
            accu=accu_array[:,i-2]

        x,y=map1(LON,LAT)

        map1.scatter(x,y, c=accu*100, marker='o', lw=0., edgecolor='', s=dotsize)
        cb=map1.colorbar(pad='12%')
        cb.set_label('accu (cm/yr)')



    xEDC,yEDC=map1(lonEDC,latEDC)
    map1.scatter(xEDC,yEDC, marker='o', c='k', s=5.)

    pp=PdfPages(MapLabel+'.pdf')
    pp.savefig(plt.figure(MapLabel))
    pp.close()
    plt.close(fig)

plt.show()
