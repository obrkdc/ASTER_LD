# -*- coding: utf-8 -*-
"""
-------------------------------------------------------------------------------
 How to Convert ASTER L1T HDF-EOS VNIR/SWIR datasets from Radiance Stored as
 Digital Numbers to GeoTIFF files of Top of Atmosphere Reflectance

 This tool imports ASTER L1T HDF-EOS files, converts data values from DN to
 radiance and TOA reflectance, georeferences, and  exports as GeoTIFF files.
-------------------------------------------------------------------------------
 Author: Cole Krehbiel
 Contact: LPDAAC@usgs.gov
 Organization: Land Processes Distributed Active Archive Center
 Date last modified: 03-06-2017
-------------------------------------------------------------------------------
 DESCRIPTION:
 This script demonstrates how to convert ASTER L1T data from Digital Number(DN)
 to radiance (w/m2/sr/µm), and from radiance into Top of Atmosphere (TOA)
 reflectance.

 The script uses an ASTER L1T HDF-EOS file (.hdf) as the input and outputs
 georeferenced tagged image file format (GeoTIFF) files for each of the VNIR
 and SWIR science datasets contained in the original ASTER L1T file.

 Results from these tutorials are output in UTM with WGS84 datum as GeoTIFFs.
 The output GeoTIFFs for each band include:
 1. Original_AST_L1T_Filename_ImageDataband#.tif
      Data is At-Sensor radiance stored as DN
 2. Original_AST_L1T_Filename_ImageDataband#_radiance.tif
      Data is At-Sensor radiance that has been converted out of DN
 3. Original_AST_L1T_Filename_ImageDataband#_reflectance.tif
      Data is TOA reflectance
 This tool will batch process ASTER L1T files if more than 1 is located in the
 working directory. Output directory will be 'inputfiledirectory'+'/output/'

 This tool was specifically developed for ASTER L1T HDF-EOS files and should
 only be used for those data products.
-------------------------------------------------------------------------------
 PREREQUISITES:
 Discaimer: This script was tested in the following environments:
   -  Python: Version 3.4.5 and Version 2.7.12 (Anaconda 4.1.1) on Windows OS
   -  Geospatial Data Abstraction Library (GDAL): Version 2.0.0 and 2.1.0

 Python Packages/Modules:
 osgeo with gdal and osr – 2.0.0, 2.1.0
 numpy – 1.11.1, 1.12.0
 argparse - 1.1
 re - 2.2.1
 os, glob, datetime, sys, getopt
-------------------------------------------------------------------------------
 ADDITIONAL INFORMATION:
 LP DAAC ASTER L1T Product Page:
 https://lpdaac.usgs.gov/dataset_discovery/aster/aster_products_table/ast_l1t
 LP DAAC ASTER L1T User Guide:
 https://lpdaac.usgs.gov/sites/default/files/public/product_documentation/
 aster_l1t_users_guide.pdf

 Search for other tools at https://lpdaac.usgs.gov/
-------------------------------------------------------------------------------
 PROCEDURES:
 1.	Copy/clone ASTERL1T_DN2REF.py from LP DAAC Recipes & Tutorials Repository
 2.	Download ASTER L1T data from the LP DAAC to a local directory
 3.	Open a Command Prompt window and navigate to the directory where you
      downloaded the ASTERL1T_DN2REF.py script
 4.	Activate python in the Command Prompt window
     1.  > activate [python environment name]
 5.	Once activated, run the script with the following in your Command Prompt:
     1.  > python ASTERL1T_DN2REF.py [insert input dir with AST_L1T files here]
             1.	Example of input directory: C:/users/johndoe/ASTERL1T/
-------------------------------------------------------------------------------
"""
# Load necessary packages into Python
from osgeo import gdal, osr
import numpy as np
from datetime import datetime
import os, glob, sys, getopt, argparse, re


# ------------------------------------------------------------------------------
# Define Script and handle errors
def main(argv):
    parser = argparse.ArgumentParser()
    try:
        opts, args = getopt.getopt(argv, "hi:", ["input_directory"])
        if len(sys.argv[1:]) == 0:
            class MyParser(argparse.ArgumentParser):
                def error(self, message):
                    sys.stderr.write('error: %s\n' % message)
                    self.print_help()
                    sys.exit(2)

            parser = MyParser()
            parser.add_argument('input_directory', nargs='+')
            args = parser.parse_args()
        elif "'" in sys.argv[1] or '"' in sys.argv[1]:
            parser.error('error: Do not include quotes in input directory argument')
        elif len(sys.argv) > 2:
            parser.error('error: Only 1 Argument is allowed (input_directory)')
        elif sys.argv[1][-1] != '/' and sys.argv[1][-1] != '\\':
            parser.error('error: Please end your directory location with / or \\')
    except getopt.GetoptError:
        print('error: Invalid option passed as argument')
        print('ASTERL1T_DN2REF.py <input_directory>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('ASTERL1T_DN2REF.py <input_directory>')
            sys.exit()
    try:
        os.chdir(sys.argv[1])
    except FileNotFoundError:
        print('error: input_directory provided does not exist or was not found')
        sys.exit(2)
    # ------------------------------------------------------------------------------
    # Set input/current working directory from user defined argument
    in_dir = sys.argv[1]
    os.chdir(in_dir)

    # Create and set output directory
    out_dir = os.path.normpath((os.path.split(in_dir)[0] + os.sep + 'Output'))
    if not os.path.exists(out_dir): os.makedirs(out_dir)

    # Create a list of ASTER L1T HDF files in the directory
    file_list = glob.glob('AST_L1T_**.hdf')

    if len(file_list) == 0:
        print('Error: no ASTER L1T hdf files were found in this directory')
        sys.exit(2)
    ucc = np.matrix(([[0.676, 1.688, 2.25, 0.0], \
                      [0.708, 1.415, 1.89, 0.0], \
                      [0.423, 0.862, 1.15, 0.0], \
                      [0.1087, 0.2174, 0.2900, 0.2900], \
                      [0.0348, 0.0696, 0.0925, 0.4090], \
                      [0.0313, 0.0625, 0.0830, 0.3900], \
                      [0.0299, 0.0597, 0.0795, 0.3320], \
                      [0.0209, 0.0417, 0.0556, 0.2450], \
                      [0.0159, 0.0318, 0.0424, 0.2650]]))

    # Thome et al. is used, which uses spectral irradiance values from MODTRAN
    # Ordered b1, b2, b3N, b4, b5...b9
    irradiance = [1848, 1549, 1114, 225.4, 86.63, 81.85, 74.85, 66.49, 59.85]

    def dn2rad(x):
        rad = (x - 1.) * ucc1
        return rad

    def rad2ref(rad):
        ref = (np.pi * rad * (esd * esd)) / (irradiance1 * np.sin(np.pi * sza /
                                                                  180))
        return ref

    # ------------------------------------------------------------------------------
    # Loop through all ASTER L1T hdf files in the directory
    for k in range(len(file_list)):

        # Maintains original filename convention
        file_name = file_list[k]

        print('Processing File: ' + file_name + ' (' + str(k + 1) + ' out of '
              + str(len(file_list)) + ')')

        # Read in the file and metadata
        aster = gdal.Open(file_name)
        aster_sds = aster.GetSubDatasets()
        meta = aster.GetMetadata()
        date = meta['CALENDARDATE']
        dated = datetime.strptime(date, '%Y%m%d')
        day = dated.timetuple()
        doy = day.tm_yday

        # Calculate Earth-Sun Distance
        esd = 1.0 - 0.01672 * np.cos(np.radians(0.9856 * (doy - 4)))

        del date, dated, day, doy

        # Need SZA--calculate by grabbing solar elevation info
        sza = [np.float(x) for x in meta['SOLARDIRECTION'].split(', ')][1]

        # Query gain data for each  band, needed for UCC
        gain_list = [g for g in meta.keys() if 'GAIN' in g]  ###### AARON HERE
        gain_info = []
        for f in range(len(gain_list)):
            gain_info1 = meta[gain_list[f]].split(', ')  # [0] ###### AARON HERE
            gain_info.append(gain_info1)
        gain_dict = dict(gain_info)

        # Define UL, LR, UTM zone
        ul = [np.float(x) for x in meta['UPPERLEFTM'].split(', ')]
        lr = [np.float(x) for x in meta['LOWERRIGHTM'].split(', ')]
        utm = np.int(meta['UTMZONENUMBER'])
        n_s = np.float(meta['NORTHBOUNDINGCOORDINATE'])

        # Create UTM zone code numbers
        utm_n = [i + 32600 for i in range(60)]
        utm_s = [i + 32700 for i in range(60)]

        # Define UTM zone based on North or South
        if n_s < 0:
            utm_zone = utm_s[utm]
        else:
            utm_zone = utm_n[utm]

        del utm_n, utm_s, utm, meta
        # ------------------------------------------------------------------------------
        # Loop through all ASTER L1T SDS (bands)
        for e in range(len(aster_sds)):
            gname = str(aster_sds[e])
            # Maintain original dataset name

            aster_sd = gname.split(',')[0]

            vnir = re.search("(VNIR.*)", aster_sd)
            swir = re.search("(SWIR.*)", aster_sd)
            if vnir or swir:
                # Generate output name for tif
                aster_sd2 = aster_sd.split('(')[1]
                aster_sd3 = aster_sd2[1:-1]
                band = aster_sd3.split(':')[-1]
                out_filename = '{}/{}_{}.tif'.format(out_dir, file_name.split('.hdf')[0], band)
                out_filename_rad = '{}_radiance.tif'.format(out_filename.split('.tif')[0])
                out_filename_ref = '{}_reflectance.tif'.format(out_filename.split('.tif')[0])
                # out_filename = out_dir + file_name.split('.hdf')[0] + '_' + band + '.tif'
                # out_filename_rad = out_filename.split('.tif')[0] + '_radiance.tif'
                # out_filename_ref = out_filename.split('.tif')[0] + '_reflectance.tif'
                # Open SDS and create array
                band_ds = gdal.Open(aster_sd3, gdal.GA_ReadOnly)
                sds = band_ds.ReadAsArray().astype(np.uint16)

                del aster_sd, aster_sd2, aster_sd3

                # Define extent and provide offset for UTM South zones
                if n_s < 0:
                    ul_y = ul[0] + 10000000
                    ul_x = ul[1]

                    lr_y = lr[0] + 10000000
                    lr_x = lr[1]

                # Define extent for UTM North zones
                else:
                    ul_y = ul[0]
                    ul_x = ul[1]

                    lr_y = lr[0]
                    lr_x = lr[1]

                # Query raster dimensions and calculate raster x & y resolution
                ncol, nrow = sds.shape
                y_res = -1 * round((max(ul_y, lr_y) - min(ul_y, lr_y)) / ncol)
                x_res = round((max(ul_x, lr_x) - min(ul_x, lr_x)) / nrow)

                # Define UL x and y coordinates based on spatial resolution
                ul_yy = ul_y - (y_res / 2)
                ul_xx = ul_x - (x_res / 2)
                # ------------------------------------------------------------------------------
                # Start conversions by band (1-9)
                if band == 'ImageData1':
                    bn = -1 + 1
                    # Query for gain specified in file metadata (by band)
                    if gain_dict['01'] == 'HGH':
                        ucc1 = ucc[bn, 0]
                    elif gain_dict['01'] == 'NOR':
                        ucc1 = ucc[bn, 1]
                    else:
                        ucc1 = ucc[bn, 2]

                if band == 'ImageData2':
                    bn = -1 + 2
                    # Query for gain specified in file metadata (by band)
                    if gain_dict['02'] == 'HGH':
                        ucc1 = ucc[bn, 0]
                    elif gain_dict['02'] == 'NOR':
                        ucc1 = ucc[bn, 1]
                    else:
                        ucc1 = ucc[bn, 2]

                if band == 'ImageData3N':
                    bn = -1 + 3
                    # Query for gain specified in file metadata (by band)
                    if gain_dict['3N'] == 'HGH':
                        ucc1 = ucc[bn, 0]
                    elif gain_dict['3N'] == 'NOR':
                        ucc1 = ucc[bn, 1]
                    else:
                        ucc1 = ucc[bn, 2]

                if band == 'ImageData4':
                    bn = -1 + 4
                    # Query for gain specified in file metadata (by band)
                    if gain_dict['04'] == 'HGH':
                        ucc1 = ucc[bn, 0]
                    elif gain_dict['04'] == 'NOR':
                        ucc1 = ucc[bn, 1]
                    elif gain_dict['04'] == 'LO1':
                        ucc1 = ucc[bn, 2]
                    else:
                        ucc1 = ucc[bn, 3]

                if band == 'ImageData5':
                    bn = -1 + 5
                    # Query for gain specified in file metadata (by band)
                    if gain_dict['05'] == 'HGH':
                        ucc1 = ucc[bn, 0]
                    elif gain_dict['05'] == 'NOR':
                        ucc1 = ucc[bn, 1]
                    elif gain_dict['05'] == 'LO1':
                        ucc1 = ucc[bn, 2]
                    else:
                        ucc1 = ucc[bn, 3]

                if band == 'ImageData6':
                    bn = -1 + 6
                    # Query for gain specified in file metadata (by band)
                    if gain_dict['06'] == 'HGH':
                        ucc1 = ucc[bn, 0]
                    elif gain_dict['06'] == 'NOR':
                        ucc1 = ucc[bn, 1]
                    elif gain_dict['06'] == 'LO1':
                        ucc1 = ucc[bn, 2]
                    else:
                        ucc1 = ucc[bn, 3]

                if band == 'ImageData7':
                    bn = -1 + 7
                    # Query for gain specified in file metadata (by band)
                    if gain_dict['07'] == 'HGH':
                        ucc1 = ucc[bn, 0]
                    elif gain_dict['07'] == 'NOR':
                        ucc1 = ucc[bn, 1]
                    elif gain_dict['07'] == 'LO1':
                        ucc1 = ucc[bn, 2]
                    else:
                        ucc1 = ucc[bn, 3]

                if band == 'ImageData8':
                    bn = -1 + 8
                    # Query for gain specified in file metadata (by band)
                    if gain_dict['08'] == 'HGH':
                        ucc1 = ucc[bn, 0]
                    elif gain_dict['08'] == 'NOR':
                        ucc1 = ucc[bn, 1]
                    elif gain_dict['08'] == 'LO1':
                        ucc1 = ucc[bn, 2]
                    else:
                        ucc1 = ucc[bn, 3]

                if band == 'ImageData9':
                    bn = -1 + 9
                    # Query for gain specified in file metadata (by band)
                    if gain_dict['09'] == 'HGH':
                        ucc1 = ucc[bn, 0]
                    elif gain_dict['09'] == 'NOR':
                        ucc1 = ucc[bn, 1]
                    elif gain_dict['09'] == 'LO1':
                        ucc1 = ucc[bn, 2]
                    else:
                        ucc1 = ucc[bn, 3]
                        # ------------------------------------------------------------------------------
                # Set irradiance value for specific band
                irradiance1 = irradiance[bn]

                # Generate output Geotiff files
                # First, Radiance (stored as DN)
                driver = gdal.GetDriverByName('GTiff')
                dn = driver.Create(out_filename, nrow, ncol, 1, gdal.GDT_UInt16)

                # Define output GeoTiff CRS and extent properties
                srs = osr.SpatialReference()
                srs.ImportFromEPSG(utm_zone)
                dn.SetProjection(srs.ExportToWkt())
                dn.SetGeoTransform((ul_xx, x_res, 0., ul_yy, 0., y_res))

                # Write SDS array to output GeoTiff
                #outband = dn.GetRasterBand(1)
                #outband.SetNoDataValue(0)
                #outband.WriteArray(sds)
                dn = None

                # Convert from DN to Radiance
                rad = dn2rad(sds)
                rad[rad == dn2rad(0)] = 0
                del sds

                # Next, Radiance (w/m2/sr/µm)
                out_rad = driver.Create(out_filename_rad, nrow, ncol, 1, \
                                        gdal.GDT_Float32)

                # Define output GeoTiff CRS and extent properties
                out_rad.SetProjection(srs.ExportToWkt())
                out_rad.SetGeoTransform((ul_xx, x_res, 0., ul_yy, 0., y_res))

                # Write SDS array to output GeoTiff
                #outband = out_rad.GetRasterBand(1)
                #outband.SetNoDataValue(0)
                #outband.WriteArray(rad)
                out_rad = None

                # Convert from Radiance to TOA Reflectance
                ref = rad2ref(rad)
                del rad

                # Lastly, Reflectance (w/m2/sr/µm)
                out_ref = driver.Create(out_filename_ref, nrow, ncol, 1, \
                                        gdal.GDT_Float32)

                # Define output GeoTiff CRS and extent properties
                out_ref.SetProjection(srs.ExportToWkt())
                out_ref.SetGeoTransform((ul_xx, x_res, 0., ul_yy, 0., y_res))

                # Write SDS array to output GeoTiff
                outband = out_ref.GetRasterBand(1)
                outband.SetNoDataValue(0)
                outband.WriteArray(ref)
                out_ref = None

                del band, bn, gname, out_filename, \
                    out_filename_rad, out_filename_ref, ref


if __name__ == "__main__":
    main(sys.argv[1:])