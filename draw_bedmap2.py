from mpl_toolkits.basemap import Basemap,cm
from matplotlib.colors import Normalize
from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
import gdal
import sys

#lat1=-75.5
#lon1=128.
#lat2=-74.8
#lon2=118.1
lat1=-75.35
lon1=126.7
lat2=-75
lon2=120.5
lonEDC=123.+21./60.
latEDC=-75.1

output_format="svg"

RLDir=sys.argv[1]
if RLDir[-1]!='/':
    RLDir=RLDir+'/'

#m = Basemap(projection='stere', lat_ts=-71, lat_0=-90, lon_0=0.,llcrnrlon=-135,llcrnrlat=-48.458667, urcrnrlon=45,urcrnrlat=-48.458667, rsphere=(6378137.00,6356752.3142))
m = Basemap(projection='stere', lat_ts=-71, lat_0=-90, lon_0=180, llcrnrlon=-135,llcrnrlat=-48.458667, urcrnrlon=45,urcrnrlat=-48.458667, rsphere=(6378137.00,6356752.3142))
#m = Basemap(projection='spstere', boundinglat=-60, lon_0=180)
m.drawcoastlines()
m.fillcontinents(color='white',lake_color='aqua')
#m.drawmapboundary(fill_color='aqua')

m.drawparallels(np.arange(-90.,81.,5.), dashes=[1, 2], color='0.5', linewidths=0.5)
m.drawmeridians(np.arange(-180.,181.,10.), latmax=85., dashes=[1, 2], color='0.5', linewidths=0.5)

#Draw bed shade
#I = plt.imread(RLDir+'bedmap2/bedmap2_bed_shade.tif')
#I=I[::-1,:]
##I=np.where(I==-9999,np.nan,I)
#m.imshow(I)

#img = Image.open(RLDir+'bedmap2/bedmap2_bed_shade.tif')
#arr = np.asarray(img)
#arr=np.where(arr==-9999,np.nan,arr)
#m.imshow(arr, cmap='Greys', origin='upper', vmin=180., vmax=220.)

#Draw bed topography
#raster = gdal.Open(RLDir+'bedmap2/bedmap2_bed.txt')
#band = raster.GetRasterBand(1)
#array = band.ReadAsArray()
#array=np.where(array==-9999,np.nan,array)
#norm = Normalize(vmin=-3000.,vmax=3000.)
#m.imshow(array, cmap='terrain', norm=norm, alpha=0.5, origin='upper')
#m.colorbar()



#Draw surface contours
raster2 = gdal.Open(RLDir+'/bedmap2_surface.txt')
band2 = raster2.GetRasterBand(1)
array2 = band2.ReadAsArray()
array2=np.where(array2==-9999,np.nan,array2)

x = np.linspace(0, m.urcrnrx, array2.shape[1])
y = np.linspace(0, m.urcrnry, array2.shape[0])

xx, yy = np.meshgrid(x, y)

m.contour(xx,yy, array2[::-1,:], colors='0.5', linewidths=0.5)

#Draw continent's contour
raster3 = gdal.Open(RLDir+'/bedmap2_icemask_grounded_and_shelves.txt')
band3 = raster3.GetRasterBand(1)
array3 = band3.ReadAsArray()
#array3=np.where(array3==-9999,np.nan,array3)

#Draw bed on Blob A
#img = Image.open(RLDir+'bedmap2/Bed_BlobA.tiff')
#arr = np.asarray(img)
##arr=np.where(arr==-9999,np.nan,arr)
#hmin=1298450.
#hmax=1391550.
#vmin=-840950.
#vmax=-888950.
#norm = Normalize(vmin=-3000.,vmax=3000.)
#m.imshow(arr, origin='upper', cmap='terrain', extent=[hmax, hmin, vmax, vmin], norm=norm)
##plt.xlim=(x1,x2)
##plt.ylim=(y1,y2)

x = np.linspace(0, m.urcrnrx, array3.shape[1])
y = np.linspace(0, m.urcrnry, array3.shape[0])
xx, yy = np.meshgrid(x, y)
m.contour(xx,yy, array3[::-1,:], colors='0.5', linewidths=0.5)

xEDC,yEDC=m(lonEDC,latEDC)
#m.plot(xEDC,yEDC, marker='*', c='r')
x1,y1=m(lon1,lat1)
x2,y2=m(lon2,lat2)
xsquare=np.array([x1,x1,x2,x2,x1])
ysquare=np.array([y1,y2,y2,y1,y1])
m.plot(xsquare,ysquare, linestyle='solid', color='r')


plt.savefig(RLDir+'Bedmap2.'+output_format, format=output_format, bbox_inches='tight')
