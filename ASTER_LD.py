# Add in modules
import glob
import sys
import rasterio as rio
from tkinter import Tcl
import os
from WBT.whitebox_tools import WhiteboxTools
wbt = WhiteboxTools()

# Set up and paths ------------------------------------------------------------------------------------------------------
wbt.work_dir = os.path.join("/SRS_Processing_Data", 'Output')
scene_path = os.path.join("/SRS_Processing_Data")
output_scene_path = os.path.join("/SRS_Processing_Data", 'Output')

# This callback function suppresses WBT printing progress updates which always use the '%' character. The callback
# function approach is flexible and allows for any level of complex interaction with tool outputs.
def my_callback(value):
   if not "%" in value:
       print(value)


# check if working directory exists
if os.path.exists(scene_path):
    print('Data directory exists, all OK')
else:
    print('Create directory (C:\SRS_Processing_Data) and add scene for processing into this folder (.hdf file)')

# ------------------------------------------------------------------------------------------------------------------
# Open hdf file, convert at sensor radiance to top of atmosphere reflectance and create .tif files for each band
os.system("python ASTERL1T_DN2REF.py C:/SRS_Processing_Data/")

#Clean up files in folder and delete those that are not reflectance
wanted_file_list = glob.glob(os.path.join(output_scene_path, '*_reflectance.tif'))
file_list = glob.glob(os.path.join(output_scene_path, '*.tif'))
for filename in file_list:
   if filename not in wanted_file_list:
       os.remove(filename)

# ------------------------------------------------------------------------------------------------------------------
scene_count = len(glob.glob1(scene_path,"*.hdf"))
print('--------------------------------')
print('--------------------------------')
print('Number of scenes in input folder =' + str(scene_count))

# Merge band groups function
def mosaic_bands(i):
    x = 0
    if (x + 1) < len(im_group_b[i]):
        wbt.mosaic_with_feathering(
            im_group_b[i][x],
            im_group_b[i][x + 1],
            os.path.join(output_scene_path, 'im_group_b' + str(i) + 'image1-' + str(x + 2) + '.tif'),
            method="cc",
            weight=4.0,
            callback=my_callback)
        x = x + 1
        while (x + 1) < len(im_group_b[i]):
            wbt.mosaic_with_feathering(
                im_group_b[i][x + 1],
                os.path.join(output_scene_path, 'im_group_b' + str(i) + 'image1-' + str(x + 1) + '.tif'),
                os.path.join(output_scene_path, 'im_group_b' + str(i) + 'image1-' + str(x + 2) + '.tif'),
                method="cc",
                weight=4.0,
                callback=my_callback)
            x = x + 1

#Prepare Image Files into groups and merge band images for more than 1 hdf file input, for 1 hdf file input rename files
if scene_count > 1:
    im_group_b = [0]*10
    for i in range(1,10):
        im_group_b[i] = glob.glob(os.path.join(output_scene_path, '*ImageData' + str(i) + '*_reflectance.tif'))
        Tcl().call('lsort', '-dict', im_group_b[i]) #Sort Images in List
        print(im_group_b[i])
    for i in range(1,10):
        mosaic_bands(i)
elif scene_count == 1:
    band_list = glob.glob(os.path.join(output_scene_path, '*reflectance.tif'))
    band_list.sort()
    for index, file in enumerate(band_list):
        os.rename(file, (os.path.join(output_scene_path, 'im_group_b' + str(index + 1) + 'image1-1.tif')))



#Clean up mosaic files in folder, delete those not needed
full_mosaic_file_list = glob.glob(os.path.join(output_scene_path, 'im_group_b*image1-' + str(scene_count) + '.tif'))
file_list1 = glob.glob(os.path.join(output_scene_path, 'im_group*' + '*.tif'))
for filename in file_list1:
   if filename not in full_mosaic_file_list:
       os.remove(filename)

# ---------------------------------------------------------------------------------------------------------------------
# Resample the SWIR bands to 15m resolution (instead of using cell_size=15 to resample, use base="reference image" as
# this then shifts the pixels to ensure they overlap exactly (rather than half a pixel out)
for i in range(4,10):
    SWIR_band_list = glob.glob(os.path.join(output_scene_path,"im_group_b" + str(i) + '*.tif'))
    for f in SWIR_band_list:
        in_file_name = f
        out_file_name = 'b' + str(i) + '_15m_WBT.tif'
        wbt.resample(in_file_name,out_file_name,cell_size=None,base=(os.path.join(output_scene_path, 'im_group_b1image1-' + str(scene_count) + '.tif')),method="nn")

# Rename bands 1-3 (already 15 m resolution) to simplify filenames
os.rename((os.path.join(output_scene_path, 'im_group_b1image1-' + str(scene_count) + '.tif')), (os.path.join(output_scene_path,'b1.tif')))
os.rename((os.path.join(output_scene_path, 'im_group_b2image1-' + str(scene_count) + '.tif')), (os.path.join(output_scene_path,'b2.tif')))
os.rename((os.path.join(output_scene_path, 'im_group_b3image1-' + str(scene_count) + '.tif')), (os.path.join(output_scene_path,'b3.tif')))

#List resampled images and sort list
resampled_band_list_path = glob.glob(os.path.join(output_scene_path, '*15m_WBT.tif'))
resampled_band_list_path.sort()
print(resampled_band_list_path)

# Create false colour composite images with band combinations below: -------------------------------------------------
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

#Alternative way of making composite image, image in correct format to do DDS enhancement
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

band1.close()
band2.close()
band3.close()
band4.close()
band5.close()
band6.close()
band7.close()
band8.close()
band9.close()

# Perform principal component analysis (PCA) on all 9 bands ----------------------------------------------------------
print("Performing 9 band PCA...")
wbt.verbose = False  # PCA report to be automatically displayed
wbt.principal_component_analysis(
   inputs='b1.tif;b2.tif;b3.tif;b4_15m_WBT.tif;b5_15m_WBT.tif;b6_15m_WBT.tif;b7_15m_WBT.tif;b8_15m_WBT.tif;b9_15m_WBT.tif',
   output="pca_report.html",
   num_comp=9,
   standardized=False
)

# Create PCA false colour composites ----------------------------------------------------------------------------------
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

PCA1.close()
PCA2.close()
PCA3.close()
PCA4.close()
PCA5.close()

# End Code ---------------------------------------------------------------------------------------------------------
sys.exit()

# Next developments:------------------------------------------------------------------------------------------------

# PCA on SWIR bands only
os.mkdir(os.path.join("/SRS_Processing_Data", 'Output', 'PCA_SWIR'))

def PC_Rename (outputFilename):
    for i in range(10):
        PCA_component_list = glob.glob(os.path.join(output_scene_path,"PCA_component" + str(i) + '.tif'))
        for f in SWIR_band_list:
            in_file_name = f
            os.rename(in_file_name, {outputFilename})

PC_Rename('PC' + str(i) + '_VNIR_SWIR.tif')



#NEED to copy across the PCA input files into the new folder for this to work or rename the PCA output files before doing the second PCA so they are not overwritten.

print("Performing PCA on SWIR bands...")
wbt.verbose = False  # PCA report to be automatically displayed
wbt.principal_component_analysis(
   inputs='b4_15m_WBT.tif;b5_15m_WBT.tif;b6_15m_WBT.tif;b7_15m_WBT.tif;b8_15m_WBT.tif;b9_15m_WBT.tif',
   output="pca_SWIR_report.html",
   num_comp=6,
   standardized=False
)

wbt.work_dir = os.path.join("/SRS_Processing_Data", 'Output')

#Balance Contrast Enhancement to reduce colour bias in a colour composite image based on Liu (1991)
# Enhancement of colour composite images, direct decorrelation stretch
wbt.direct_decorrelation_stretch(
    'ASTER_721_WBT.tiff',
    'ASTER_721_WBT_DDS.tiff',
    k=0.5,
    clip=1.0,
    callback=my_callback
)
