from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
import gdal
import sys

RLDir=sys.argv[1]
if RLDir[-1]!='/':
    RLDir=RLDir+'/'


m = Basemap(projection='stere', lat_ts=-71, lat_0=-90, lon_0=180, llcrnrlon=-135,llcrnrlat=-48.458667, urcrnrlon=45,urcrnrlat=-48.458667, rsphere=(6378137.00,6356752.3142))
#m = Basemap(projection='spstere', boundinglat=-60, lon_0=180)
#m.drawcoastlines()
#m.fillcontinents(color='white',lake_color='aqua')
#m.drawmapboundary(fill_color='aqua')

m.drawparallels(np.arange(-90.,81.,5.))
m.drawmeridians(np.arange(-180.,181.,10.), latmax=85.)

#Draw bed topography
raster = gdal.Open(RLDir+'bedmap2/bedmap2_bed.txt')
band = raster.GetRasterBand(1)
array = band.ReadAsArray()
array=np.where(array==-9999,np.nan,array)

m.imshow(array[::-1,:])
m.colorbar()


#Draw surface contours
raster2 = gdal.Open(RLDir+'bedmap2/bedmap2_surface.txt')
band2 = raster2.GetRasterBand(1)
array2 = band2.ReadAsArray()
array2=np.where(array2==-9999,np.nan,array2)

x = np.linspace(0, m.urcrnrx, array2.shape[1])
y = np.linspace(0, m.urcrnry, array2.shape[0])

xx, yy = np.meshgrid(x, y)

m.contour(xx,yy, array2[::-1,:], colors='k')

#Draw continent's contour
raster3 = gdal.Open(RLDir+'bedmap2/bedmap2_icemask_grounded_and_shelves.txt')
band3 = raster3.GetRasterBand(1)
array3 = band3.ReadAsArray()
#array3=np.where(array3==-9999,np.nan,array3)

x = np.linspace(0, m.urcrnrx, array3.shape[1])
y = np.linspace(0, m.urcrnry, array3.shape[0])
xx, yy = np.meshgrid(x, y)
m.contour(xx,yy, array3[::-1,:], colors='k')

m.scatter(-75., 123.)


plt.show()
