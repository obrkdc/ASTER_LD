# ASTER_LD
 
**1	ASTER LD OVERVIEW**
ASTER Lithological Discrimination (ASTER LD) allows users to quickly produce a range of  composite false colour images highlighting differences between rock types for use in geological mapping.  
ASTER LD reads Advanced Spaceborne Thermal Emission and Reflection Radiometer (ASTER) multispectral image Hierarchical Data Format (.hdf) files and processes the visible near-infrared (VNIR) and short-wave infrared (SWIR) bands to output false colour composites including principal component analysis (PCA) images suitable for lithological discrimination.
ASTER LD can process and merge data from multiple overlapping ASTER scenes, the main steps are summarised as:
1.	Read ASTER L1T HDF-EOS VNIR/SWIR datasets and convert at sensor radiance to top of atmosphere (TOA) reflectance, utilising the ASTERL1T_DN2REF.py tool from the Land Processes Distributed Active Archive Center (LPDAA) at the United States Geological Survey (USGS), (Krehbiel, 2017).
2.	Merge individual band rasters if there are multiple ASTER scenes, to output 9 merged greyscale rasters, one for each band.
3.	Resample the SWIR bands to match VNIR band extent and spatial resolution (15 m).
4.	Produce a set of false colour red-green-blue (RGB) composite images combining bands 468, 721, and 631 for lithological discrimination and 321 for vegetation mapping.
5.	Perform principal component analysis (PCA) to highlight uncorrelated components of the 9 spectral bands and use the principal components to produce false colour (RGB) composite images, PCA123, PCA321, PCA432 and PCA542.
All output images are in georeferenced tagged image file format (GeoTiff, .tif) suitable for viewing, enhancement and interpretation in geographic information system (GIS) software such as QGIS.
In addition to the TOA reflectance conversion tool, ASTER LD utilises a number of functions from the WhiteboxTools module and acknowledges their author Dr John Lindsay, (Lindsay, 2020).

**2	INSTALLATION AND SETUP**
ASTER LD was tested in Python 3.8 on Windows 10 operating system, in an Anaconda environment with the following Python packages/modules:
•	glob, os, datetime, sys, getopt, re
•	rasterio 1.2.3
•	gdal 3.3 
•	numpy 1.20.2
•	argparse 1.4.0
•	whiteboxtools 1.4.0

**2.1	Setup Requirements**
•	Create a folder with the following or operating system (OS) comparable directory path “C:\SRS_Processing_Data”. 
•	Copy ASTER Level 1 Precision Terrain Corrected Registered At-Sensor Radiance (AST_L1T)  scene files in HDF4_EOS (.hdf) format for processing into this folder. 
o	Multiple scenes should overlap by at least 100 pixels (3 km).
o	Ideally multiple scenes should be collected at similar times of the year/season.
o	Scenes can be downloaded from USGS Earth Explorer (https://earthexplorer.usgs.gov/). Test scenes used included AST_L1T_00303052001084132_20150501094701_93840.hdf and AST_L1T_00309222001083017_20150419191238_32915.hdf in Eastern Egypt. These scenes can be located in Earth Explorer by selecting the ASTER Level 1T data set and searching by the filename (without .hdf) in the ‘Entity ID’ box in the ‘Additional Criteria’ tab. 
•	Run the ASTER_LD.py script:
o	Clone/copy ASTER_LD.py and ASTERL1T_DN2REF.py from the following GitHub repository https://github.com/obrkdc/ASTER_LD into a local directory.
o	Install a Python 3 interpreter and an optional distribution such as Anaconda (which has a graphical user interface, Anaconda Navigator) to assist in loading modules and packages and setting up the environment.
o	Load the Python modules/packages listed above into a working environment or install directly in OS. WhiteBoxTools package can be downloaded from https://jblindsay.github.io/ghrg/WhiteboxTools/download.html. 
o	Start python and run the ASTER_LD.py file within the working environment/ command prompt.
 
**3	OUTPUTS**
Output files produced by ASTER LD are listed in Table 3 1. Note: output files total size are approximately ~6 Gb per input scene, depending on scene overlap.
Table 3 1:	Output File Nomenclature
File Name Format	Details
AST_L1T…    …ImageData1_reflectance.tif	TOA reflectance. One set of 9 images per input scene. Number in filename after ImageData = band number.
im_group_b5image1-3.tif	Merged TOA reflectance data for each band b5 = band number; image1-3 = merged scenes i.e. scenes 1 to 3 merged.
b4_15m_WBT.tif	Merged band resampled to 15 m spatial resolution. Bands 4-9 only.
b2.tif	Merged band renamed. Bands 1-3 only (15 m native resolution)
pca_report.html	Principal component analysis report, detailing variance and factor loadings.
PCA_component6.tif	Merged Principal Component (PC) image, set of 9 PC’s, numbered 1-9
ASTER_721.tiff	Merged false colour composite, 721 = band 7 in red channel, band 2 in green and band 1 in blue (RGB) 
ASTER_PCA_321.tiff	Merged false colour composite image using PC’s, 321 = PC3, PC2, PC1 in RGB channels.

