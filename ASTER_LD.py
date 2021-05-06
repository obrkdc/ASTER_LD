# Add in modules
import os
import warnings
import glob
from typing import Union
import fiona
import cv2
from matplotlib import pyplot as plt
import numpy.ma as ma
import xarray as xr
import rioxarray as rxr
from geopandas import GeoDataFrame
from pandas import DataFrame
from shapely.geometry import mapping, box
import geopandas as gpd
import earthpy as et
import earthpy.plot as ep
import sys
from tkinter import Tcl
import gdal
from gdalconst import GA_ReadOnly
import rasterio as rio
from rasterio.plot import plotting_extent
import earthpy.spatial as es
import os
import numpy as np
import rasterio as rio
from rasterio.plot import show
from rasterio.mask import mask
from shapely.geometry import mapping
import matplotlib.pyplot as plt
import geopandas as gpd
import earthpy as et
import earthpy.plot as ep
import earthpy.spatial as es
from WBT.whitebox_tools import WhiteboxTools
wbt = WhiteboxTools()
wbt.work_dir = os.path.join("/SRS_Processing_Data", 'Output')

scene_path = os.path.join("/SRS_Processing_Data")
output_scene_path = os.path.join("/SRS_Processing_Data", 'Output')

# This callback function suppresses WBT printing progress updates,
# which always use the '%' character. The callback function
# approach is flexible and allows for any level of complex
# interaction with tool outputs.
def my_callback(value):
   if not "%" in value:
       print(value)

# can be used once the code is written and working to stop warnings popping up to the user
# warnings.simplefilter('ignore')

# check if working directory exists
if os.path.exists(scene_path):
    print('Data directory exists, all OK')
else:
    print('Create directory (C:\SRS_Processing_Data) and add multispectral scene for processing into this folder (.hdf file)')

# Open hdf file, convert at sensor radiance to top of atmosphere reflectance and create .tif files for each band
os.system("python ASTERL1T_DN2REF.py C:/SRS_Processing_Data/")

#Clean up files in folder and delete those that are not reflectance
wanted_file_list = glob.glob(os.path.join(output_scene_path, '*_reflectance.tif'))
file_list = glob.glob(os.path.join(output_scene_path, '*.tif'))
for filename in file_list:
   if filename not in wanted_file_list:
       os.remove(filename)

#Prepare Image Files into groups for merging
im_group_b = [0]*10
for i in range(1,10):
    im_group_b[i] = glob.glob(os.path.join(output_scene_path, '*ImageData' + str(i) + '*_reflectance.tif'))
    Tcl().call('lsort', '-dict', im_group_b[i]) #Sort Images in List
    print(im_group_b[i])

def mosaic_bands(i):
    x = 0
    if (x+1) < len(im_group_b[i]):
        wbt.mosaic_with_feathering(
            im_group_b[i][x],
            im_group_b[i][x+1],
            os.path.join(output_scene_path, 'im_group_b' + str(i) + 'image1-' + str(x+2) + '.tif'),
            method="cc",
            weight=4.0,
            callback=my_callback)
        x = x + 1
        while (x+1) < len(im_group_b[i]):
            wbt.mosaic_with_feathering(
                im_group_b[i][x + 1],
                os.path.join(output_scene_path, 'im_group_b' + str(i) + 'image1-' + str(x + 1) + '.tif'),
                os.path.join(output_scene_path, 'im_group_b' + str(i) + 'image1-' + str(x + 2) + '.tif'),
                method="cc",
                weight=4.0,
                callback=my_callback)
            x = x + 1

for i in range(1,10):
    mosaic_bands(i)

scene_count = len(glob.glob1(scene_path,"*.hdf"))
print('--------------------------------')
print('--------------------------------')
print('Number of scenes in input folder =' + str(scene_count))

#Clean up mosaic files in folder, delete those not needed
full_mosaic_file_list = glob.glob(os.path.join(output_scene_path, 'im_group_b*image1-' + str(scene_count) + '.tif'))
file_list1 = glob.glob(os.path.join(output_scene_path, 'im_group*' + '*.tif'))
for filename in file_list1:
   if filename not in full_mosaic_file_list:
       os.remove(filename)

#Histogram matching between images before merging - this was supposed to make the image contrast comparable between
#imaged before merging them, but in scenes with varying amount of water (ocean) the histogram match changes lancover
#pixels too much, forming images that are further from each other than previously.
#May want to try this over areas with no ocean?
"""x = 0
while x < len(b1_images_list):
    wbt.histogram_matching_two_images(
        b1_images_list[x],
        b1_images_list[0],
        (os.path.join(output_scene_path, 'Band1_' + str(x + 1) + 'HistMatched.tif')),
        callback=my_callback)
    x = x + 1"""

#Show a plot of the hisogram of the pixel values
#img = cv2.imread((os.path.join(output_scene_path,'AST_L1T_00309222001083017_20150419191238_32915_ImageData1_reflectance.tif')), cv2.IMREAD_ANYDEPTH)
#plt.hist(img.ravel(),256,[0,1]); plt.show()






#Make a list of all the reflectance.tif file names and sort by band number
#ASTER_band_list = [i for i in os.listdir(output_scene_path) if os.path.isfile(os.path.join(output_scene_path,i)) and \
                  # 'reflectance' in i]
#ASTER_band_list.sort()
#print(ASTER_band_list)

# Resample the SWIR bands to 15m resolution (instead of using cell_size=15 to resample, use base="reference image" as this then shifts the pixels to ensure they overlap exactly (rather than half a pixel out) this took me days to try and match the different band image sizes and origins.
for i in range(4,10):
    SWIR_band_list = glob.glob(os.path.join(output_scene_path,"im_group_b" + str(i) + '*.tif'))
    #in_file_name = ("*" + str(i) + '_reflectance.tif') #Dont know why the wildcard '*' wont work here? Its called a shell wildcard
    for f in SWIR_band_list:
        in_file_name = f
        out_file_name = 'b' + str(i) + '_15m_WBT.tif'
        wbt.resample(in_file_name,out_file_name,cell_size=None,base=(os.path.join(output_scene_path, 'im_group_b1image1-' + str(scene_count) + '.tif')),method="nn")
# NEED TO SET PROJECTION

os.rename((os.path.join(output_scene_path, 'im_group_b1image1-' + str(scene_count) + '.tif')), (os.path.join(output_scene_path,'b1.tif')))
os.rename((os.path.join(output_scene_path, 'im_group_b2image1-' + str(scene_count) + '.tif')), (os.path.join(output_scene_path,'b2.tif')))
os.rename((os.path.join(output_scene_path, 'im_group_b3image1-' + str(scene_count) + '.tif')), (os.path.join(output_scene_path,'b3.tif')))


resampled_band_list_path = glob.glob(os.path.join(output_scene_path, '*15m_WBT.tif'))
resampled_band_list_path.sort()
print(resampled_band_list_path)

# Simplify names for bands
b1 = (os.path.join(output_scene_path,'b1.tif'))
b2 = (os.path.join(output_scene_path,'b2.tif'))
b3 = (os.path.join(output_scene_path,'b3.tif'))
b4 = resampled_band_list_path [0]
b5 = resampled_band_list_path [1]
b6 = resampled_band_list_path [2]
b7 = resampled_band_list_path [3]
b8 = resampled_band_list_path [4]
b9 = resampled_band_list_path [5]



"""# Open bands in GDAL as datasets 'db1'
db1 = gdal.Open(b1)
db2 = gdal.Open(b2)
db3 = gdal.Open(b3)
db4 = gdal.Open(b4)
db5 = gdal.Open(b5)
db6 = gdal.Open(b6)
db7 = gdal.Open(str(output_scene_path) + '\output.tif')
db8 = gdal.Open(b8)
db9 = gdal.Open(b9)

#After hours of trying to get this to work realised this was not for Python
#gdalbuildvrt separate, stack.vrt, b7, b2, b1
#plt.imshow stack.vrt
#gdal_translate (stack.vrt, 7-2-1.tif)

print('Number of Raster Bands in image = ' + str(db1.RasterCount))

# False Colour Composite 7-2-1
# Fetch bands
band7 = db7.GetRasterBand(1)
band2 = db2.GetRasterBand(1)
band1 = db1.GetRasterBand(1)

# Read bands as Numpy arrays
a7 = band7.ReadAsArray()
a2 = band2.ReadAsArray()
a1 = band1.ReadAsArray()

# Plot the arrays
img721 = np.dstack((a7, a2, a1))
f = plt.figure()
plt.imshow(img721)
plt.savefig('721_fig.tif') # Saves low resolution figure of the plot
plt.show()
"""
#open all the bands in rasterio
band1 = rio.open(b1)
band2 = rio.open(b2)
band3 = rio.open(b3)
band4 = rio.open(b4)
band5 = rio.open(b5)
band6 = rio.open(b6)
band7 = rio.open(b7)
band8 = rio.open(b8)
band9 = rio.open(b9)
print('Band 1 metadata ', band1.meta)
print('Band 7 Metadata ', band7.meta)

#468 composite for hydrated mineral discrimination (hydrothermal alteration)
_468=rio.open(os.path.join(output_scene_path,'ASTER_468.tiff'), 'w', driver='Gtiff',
                        width=band1.width, height=band1.height,
                        count=3,
                        crs=band1.crs,
                        transform=band1.transform,
                        dtype='float32',
                        nodata= 0.0)
_468.write(band4.read(1),1)
_468.write(band6.read(1),2)
_468.write(band8.read(1),3)
_468.close()

#721 composite for lithological discrimination in arid to semi-arid environments
_721=rio.open(os.path.join(output_scene_path,'ASTER_721.tiff'), 'w', driver='Gtiff',
                        width=band1.width, height=band1.height,
                        count=3,
                        crs=band1.crs,
                        transform=band1.transform,
                        dtype='float32',
                        nodata= 0.0)
_721.write(band7.read(1),1)
_721.write(band2.read(1),2)
_721.write(band1.read(1),3)
_721.close()

wbt.create_colour_composite(
    'b7_15m_WBT.tif',
    b2,
    b1,
    'ASTER_721_WBT.tiff',
    opacity=None,
    enhance=True,
    zeros=True,
    callback=my_callback
)

#631 composite for lithological discrimination
_631=rio.open(os.path.join(output_scene_path,'ASTER_631.tiff'), 'w', driver='Gtiff',
                        width=band1.width, height=band1.height,
                        count=3,
                        crs=band1.crs,
                        transform=band1.transform,
                        dtype='float32',
                        nodata= 0.0)
_631.write(band6.read(1),1)
_631.write(band3.read(1),2)
_631.write(band1.read(1),3)
_631.close()

#321 composite for vegetation mapping, 'standard' NIR composite
_321=rio.open(os.path.join(output_scene_path,'ASTER_321.tiff'), 'w', driver='Gtiff',
                        width=band1.width, height=band1.height,
                        count=3,
                        crs=band1.crs,
                        transform=band1.transform,
                        dtype='float32',
                        nodata= 0.0)
_321.write(band3.read(1),1)
_321.write(band2.read(1),2)
_321.write(band1.read(1),3)
_321.close()

"""# Resample bands 4-9 from 30m to 15m resolution, to match bands 1-3, does not shift the pixel location, but new image is half a pixel larger compared to original band!?!?
gdal.Warp((os.path.join(output_scene_path, 'b4_15m.tif')), b4, xRes=15, yRes=15) # Output name, input name
gdal.Warp((str(output_scene_path) + '/b5_15m.tif'), b5, xRes=15, yRes=15)
gdal.Warp((str(output_scene_path) + '/b6_15m.tif'), b6, xRes=15, yRes=15)
gdal.Warp((str(output_scene_path) + '/b7_15m.tif'), b7, xRes=15, yRes=15)
gdal.Warp((str(output_scene_path) + '/b8_15m.tif'), b8, xRes=15, yRes=15)
gdal.Warp((str(output_scene_path) + '/b9_15m.tif'), b9, xRes=15, yRes=15)
"""







"""# Crop/clip images to b1 size
# First produce a shp file for clipping
#get crs of the raster data to set shp file crs
with rio.open(b1) as raster_crs:
    clip_raster_profile = raster_crs.profile

# Added in the CRS and adjusted from code to produce outlines for multiple rasters from:
# https://gis.stackexchange.com/questions/281986/is-there-a-python-equivalent-of-gdaltindex)
df = gpd.GeoDataFrame(columns=['location','geometry'])
df.crs = clip_raster_profile["crs"]
bounds = rio.open(b1).bounds
fname = b1
df = df.append({'location':fname, 'geometry': box(bounds[0], bounds[1], bounds[2], bounds[3])},ignore_index=True)
df.to_file(os.path.join(output_scene_path, 'b1_outline.shp'))


# Open clip extent shp file in geopandas
clip_extent = gpd.read_file(os.path.join(output_scene_path, 'b1_outline.shp'))


# read imagery file
with rio.open(os.path.join(output_scene_path, 'b4_15m.tif')) as src:
    out_image, out_meta = es.crop_image(src, clip_extent)

# Update the metadata to have the new shape (x and y and affine information)
out_meta.update({"driver": "GTiff",
                 "height": out_image.shape[1],
                 "width": out_image.shape[2],
                 "transform": out_meta['transform']})

with rio.open(os.path.join(output_scene_path, 'b4_15m_clip.tif'), "w", **out_meta) as dest:
    dest.write(out_image)
"""

# Clips output raster to the same dimensions as the input raster but does not take into account the origin
"""data = gdal.Open(b1, GA_ReadOnly)
geoTransform = data.GetGeoTransform()
minx = geoTransform[0]
maxy = geoTransform[3]
maxx = minx + geoTransform[1] * data.RasterXSize
miny = maxy + geoTransform[5] * data.RasterYSize
os.system('gdal_translate -projwin ' + ' '.join([str(x) for x in [minx, maxy, maxx, miny]]) + ' -of GTiff ' + (os.path.join(output_scene_path, 'b4_15m.tif')) + ' ' + (os.path.join(output_scene_path, 'b4_15m_clip.tif')))"""

""" #Another version, output variable throws up error?
clipArea = gdal.Open(b1, GA_ReadOnly)
projection = clipArea.GetProjectionRef()
geoTransform = clipArea.GetGeoTransform()
minx = geoTransform[0]
maxy = geoTransform[3]
maxx = minx + geoTransform[1] * clipArea.RasterXSize
miny = maxy + geoTransform[5] * clipArea.RasterYSize

input = gdal.Open(os.path.join(output_scene_path, 'b4_15m.tif', GA_ReadOnly) #Your data the one you want to clip ##need to add in here the path to the file
output = os.path.join(output_scene_path, 'b4_15m_clip.tif') #output file
gdal.Translate(output, input, format='GTiff', projWin=[minx, maxy, maxx, miny], outputSRS=projection)"""


print("Performing PCA...")
wbt.verbose = False  # PCA report to be automatically displayed
wbt.principal_component_analysis(
   inputs='b1.tif;b2.tif;b3.tif;b4_15m_WBT.tif;b5_15m_WBT.tif;b6_15m_WBT.tif;b7_15m_WBT.tif;b8_15m_WBT.tif;b9_15m_WBT.tif',
   output="pca_report.html",
   num_comp=9,
   standardized=False
)
# Open PCA images in rasterio
PCA1 = rio.open(os.path.join(output_scene_path, 'PCA_component1.tif'))
PCA2 = rio.open(os.path.join(output_scene_path, 'PCA_component2.tif'))
PCA3 = rio.open(os.path.join(output_scene_path, 'PCA_component3.tif'))
PCA4 = rio.open(os.path.join(output_scene_path, 'PCA_component4.tif'))
PCA5 = rio.open(os.path.join(output_scene_path, 'PCA_component5.tif'))

#PCA 321 colour composite
PCA_321=rio.open(os.path.join(output_scene_path, 'ASTER_PCA_321.tiff'), 'w', driver='Gtiff',
                        width=band1.width, height=band1.height,
                        count=3,
                        crs=band1.crs,
                        transform=band1.transform,
                        dtype='float32',
                        nodata= 0.0)
PCA_321.write(PCA3.read(1),1)
PCA_321.write(PCA2.read(1),2)
PCA_321.write(PCA1.read(1),3)
PCA_321.close()

#PCA 123 colour composite
PCA_123=rio.open(os.path.join(output_scene_path, 'ASTER_PCA_123.tiff'), 'w', driver='Gtiff',
                        width=band1.width, height=band1.height,
                        count=3,
                        crs=band1.crs,
                        transform=band1.transform,
                        dtype='float32',
                        nodata= 0.0)
PCA_123.write(PCA1.read(1),1)
PCA_123.write(PCA2.read(1),2)
PCA_123.write(PCA3.read(1),3)
PCA_123.close()

#PCA 432 colour composite
PCA_123=rio.open(os.path.join(output_scene_path, 'ASTER_PCA_432.tiff'), 'w', driver='Gtiff',
                        width=band1.width, height=band1.height,
                        count=3,
                        crs=band1.crs,
                        transform=band1.transform,
                        dtype='float32',
                        nodata= 0.0)
PCA_123.write(PCA4.read(1),1)
PCA_123.write(PCA3.read(1),2)
PCA_123.write(PCA2.read(1),3)
PCA_123.close()

#PCA 542 colour composite
PCA_123=rio.open(os.path.join(output_scene_path, 'ASTER_PCA_542.tiff'), 'w', driver='Gtiff',
                        width=band1.width, height=band1.height,
                        count=3,
                        crs=band1.crs,
                        transform=band1.transform,
                        dtype='float32',
                        nodata= 0.0)
PCA_123.write(PCA5.read(1),1)
PCA_123.write(PCA4.read(1),2)
PCA_123.write(PCA2.read(1),3)
PCA_123.close()

"""os.mkdir(os.path.join("/SRS_Processing_Data", 'Output', 'PCA_SWIR'))

wbt.work_dir = os.path.join("/SRS_Processing_Data", 'Output', 'PCA_SWIR')

#NEED to copy across the PCA input files into the new folder for this to work or rename the PCA output files before doing the second PCA so they are not overwritten.

print("Performing PCA SWIR...")
wbt.verbose = False  # PCA report to be automatically displayed
wbt.principal_component_analysis(
   inputs='b4_15m_WBT.tif;b5_15m_WBT.tif;b6_15m_WBT.tif;b7_15m_WBT.tif;b8_15m_WBT.tif;b9_15m_WBT.tif',
   output="pca_SWIR_report.html",
   num_comp=6,
   standardized=False
)

wbt.work_dir = os.path.join("/SRS_Processing_Data", 'Output')"""

#Balance Contrast Enhancement to reduce colour bias in a colour composite image based on Liu (1991)


wbt.direct_decorrelation_stretch(
    'ASTER_721_WBT.tiff',
    'ASTER_721_WBT_DDS.tiff',
    k=0.5,
    clip=1.0,
    callback=my_callback
)
