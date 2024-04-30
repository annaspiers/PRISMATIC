# from Shashi Konduri
# adapted by Anna Spiers April 2024

import pandas as pd
import argparse
import os, sys
import glob
import ray
import rasterio
import numpy as np
import hytools as ht
import json
from hytools.io.envi import WriteENVI
import warnings
from hytools.io.envi import *
from hytools.topo import calc_topo_coeffs
from hytools.brdf import calc_brdf_coeffs
from hytools.glint import set_glint_parameters
from hytools.masks import mask_create
import h5py #ais for move_subset
import shutil #ais for move_subset
import time


# from Shashi Konduri
# adapted by Anna Spiers April 2024

import pandas as pd
import argparse
import os, sys
import glob
import ray
import rasterio
import numpy as np
import hytools as ht
import json
from hytools.io.envi import WriteENVI
import warnings
from hytools.io.envi import *
from hytools.topo import calc_topo_coeffs
from hytools.brdf import calc_brdf_coeffs
from hytools.glint import set_glint_parameters
from hytools.masks import mask_create
import h5py #ais for move_subset
import shutil #ais for move_subset
import time



########################## Part 1: Convert NEON HDFs to ENVI files (NO CHANGE NEEDED) #########################

def move_red_yellow_subset(site, path): #ais adapted function from https://stackoverflow.com/questions/61963546/how-to-conditionally-move-files-from-one-directory-to-another
    dest = os.path.join(path, 'cloud_condition_red_yellow')
    if not os.path.exists(dest):
        os.makedirs(dest)
    count = 0
    for file in os.listdir(path):
        if file.endswith('.h5'):
            f = h5py.File(os.path.join(path, file), 'r')
            if f[site]['Reflectance']['Reflectance_Data'].attrs['Cloud conditions'] != 'Green (<10%) cloud cover':
                target = os.path.join(dest, file)
                if os.path.exists(target):
                    print('file ' + file + ' already moved')
                shutil.move(os.path.join(path, file), os.path.join(dest, file))
                print(file, f[site]['Reflectance']['Reflectance_Data'].attrs['Cloud conditions'],
                f[site]['Reflectance']['Reflectance_Data'].attrs['Cloud type'])
                count += 1
            print('moved: %d' % count)
    print('Red/yellow file movement finished')

@ray.remote 
def convert_hdf_to_envi(path):
    '''This function exports NEON AOP HDF imaging spectroscopy data
    to an ENVI formated binary file, with the option of also exporting
    ancillary data following formatting used by NASA JPL for AVIRIS
    observables. The script utilizes ray to export images in parralel.
    '''
    print('begin convert_hdf_to_envi')
    images = glob.glob("*.h5") 
    output_dir = path
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if ray.is_initialized():
        ray.shutdown()
    ray.init(num_cpus = len(images))

    hytool = ray.remote(ht.HyTools)
    actors = [hytool.remote() for image in images]
    _ = ray.get([a.read_file.remote(image,'neon') for a,image in zip(actors,images)])

    def neon_to_envi(hy_obj):
        basemame = os.path.basename(os.path.splitext(hy_obj.file_name)[0])
        print("Exporting %s " % basemame)
        output_name = output_dir + basemame
        writer = WriteENVI(output_name,hy_obj.get_header())
        iterator = hy_obj.iterate(by = 'chunk')
        pixels_processed = 0
        while not iterator.complete:
            chunk = iterator.read_next()
            pixels_processed += chunk.shape[0]*chunk.shape[1]
            writer.write_chunk(chunk,iterator.current_line,iterator.current_column)
            if iterator.complete:
                writer.close()

    def export_anc(hy_obj):
        anc_header = hy_obj.get_header()
        anc_header['bands'] = 10
        anc_header['band_names'] = ['path length', 'to-sensor azimuth',
                                    'to-sensor zenith','to-sun azimuth',
                                      'to-sun zenith','phase', 'slope',
                                      'aspect', 'cosine i','UTC time']
        anc_header['wavelength units'] = np.nan
        anc_header['wavelength'] = np.nan
        anc_header['data type'] = 4

        output_name = output_dir + os.path.basename(os.path.splitext(hy_obj.file_name)[0])
        writer = WriteENVI(output_name + "_ancillary", anc_header)
        writer.write_band(hy_obj.get_anc("path_length"),0)
        writer.write_band(hy_obj.get_anc("sensor_az",radians = False),1)
        writer.write_band(hy_obj.get_anc("sensor_zn",radians = False),2)
        writer.write_band(hy_obj.get_anc("solar_az",radians = False),3)
        writer.write_band(hy_obj.get_anc("solar_zn",radians = False),4)
        #writer.write_band(hy_obj.get_anc("phase placeholder"),5)
        writer.write_band(hy_obj.get_anc("slope",radians = False),6)
        writer.write_band(hy_obj.get_anc("aspect",radians = False),7)
        writer.write_band(hy_obj.cosine_i(),8)
        #writer.write_band('UTC time placeholder',9)
        writer.close()
    
    _ = ray.get([a.do.remote(neon_to_envi) for a in actors])
    
    print("\nExporting ancillary data")
    _ = ray.get([a.do.remote(export_anc) for a in actors])

    print("Export convert_hdf_to_envi complete.") #ais takes 20-25min for SOAP

########################## Part 2: Create config file with all the parameters for BRDF correction (TWEAK AS NECESSARY) #################################

def create_config(path):
    '''Template script for generating image_correct configuration JSON files.
       These settiing are meant only as an example, are not appropriate for
       all situations and may need to be adjusted
    '''
    print('begin create_config')

   # Output path for configuration file
    config_file = path + "config.json"
    config_dict = {}
   
   # Only coefficients for good bands will be calculated
    config_dict['bad_bands'] = [] #[[300,400],[1337,1430],[1800,1960],[2450,2600]]   ####### Change this as needed #######
   #ais make blank so that zenith tif will be 427th layer

   # Input data settings for ENVI
    config_dict['file_type'] = 'neon'
    aviris_anc_names = ['path_length','sensor_az','sensor_zn',
                        'solar_az', 'solar_zn','phase','slope',
                        'aspect', 'cosine_i','utc_time']
    images= glob.glob("NEON*reflectance")
    images.sort()
    config_dict["input_files"] = images
    config_dict["anc_files"] = {}
    anc_files = glob.glob(path + "*ancillary")
    anc_files.sort()
    for i,image in enumerate(images):
        config_dict["anc_files"][image] = dict(zip(aviris_anc_names,
                                                               [[anc_files[i],a] for a in range(len(aviris_anc_names))]))
   
   # Export settings
    ''' Options for subset waves:
       1. List of subset wavelenths
       2. Empty list, this will output all good bands, if resampler is
       set it will also resample.
       - Currently resampler cannot be used in conjuction with option 1
    '''
   
    config_dict['export'] = {}
    config_dict['export']['coeffs']  = False
    config_dict['export']['image']  = True
    config_dict['export']['masks']  = False
    config_dict['export']['subset_waves']  = []
    config_dict['export']['output_dir'] = path
    config_dict['export']["suffix"] = 'BRDF_topo_corrected'
   
   #Corrections
    ''' Specify correction(s) to be applied, corrections will be applied
   in the order they are specified.
   Options include:
       ['topo']
       ['brdf']
       ['glint']
       ['topo','brdf']
       ['brdf','topo']
       ['brdf','topo','glint']
       [] <---Export uncorrected images
    '''
   
    config_dict["corrections"] = ['brdf','topo','glint']
   
    #Topographic Correction options
    '''
    Types supported:
        - 'cosine'
        - 'c'
        - 'scs
        - 'scs+c'
        - 'mod_minneart'
        - 'precomputed'
    Apply and calc masks are only needed for C and SCS+C corrections. They will
    be ignored in all other cases and correction will be applied to all
    non no-data pixels.
    'c_fit_type' is only applicable for the C or SCS+C correction type. Options
    include 'ols' or 'nnls'. Choosing 'nnls' can limit overcorrection.
    For precomputed topographic coefficients 'coeff_files' is a
    dictionary where each key is the full the image path and value
    is the full path to coefficients file, one per image.
    '''
    config_dict["topo"] =  {}
    config_dict["topo"]['type'] =  'scs+c'
    config_dict["topo"]['calc_mask'] = [["ndi", {'band_1': 850,'band_2': 660,
                                                'min': 0.1,'max': 1.0}]]
   
    config_dict["topo"]['apply_mask'] = [["ndi", {'band_1': 850,'band_2': 660,
                                                'min': 0.1,'max': 1.0}]] #could also add something about cloud_mask - check hytools GitHub, 
    config_dict["topo"]['c_fit_type'] = 'nnls'
   
    #BRDF Correction options
    '''
    Types supported:
       - 'universal': Simple kernel multiplicative correction.
       - 'local': Correction by class. (Future.....)
       - 'flex' : Correction by NDVI class
       - 'precomputed' : Use precomputed coefficients
    If 'bin_type' == 'user'
    'bins' should be a list of lists, each list the NDVI bounds [low,high]
    Object shapes ('h/b','b/r') only needed for Li kernels.
    For precomputed topographic coefficients 'coeff_files' is a
    dictionary where each key is the full the image path and value
    is the full path to coefficients file, one per image.
    '''
   
    config_dict["brdf"]  = {}
    
    # Options are 'line','scene', or a float for a custom solar zn
    # Custom solar zenith angle should be in radians
    config_dict["brdf"]['solar_zn_type'] ='scene'
    
    # Flex BRDF configs
    config_dict["brdf"]['type'] =  'flex'
    config_dict["brdf"]['grouped'] =  True
    config_dict["brdf"]['geometric'] = 'li_sparse'
    config_dict["brdf"]['volume'] = 'ross_thick'
    config_dict["brdf"]["b/r"] = 10
    config_dict["brdf"]["h/b"] = 2
    config_dict["brdf"]['sample_perc'] = 0.1
    config_dict["brdf"]['interp_kind'] = 'linear'
    config_dict["brdf"]['calc_mask'] = [["ndi", {'band_1': 850,'band_2': 660,
                                                    'min': 0.1,'max': 1.0}]]
    config_dict["brdf"]['apply_mask'] = [["ndi", {'band_1': 850,'band_2': 660,
                                                 'min': 0.1,'max': 1.0}]]
   # Flex dynamic NDVI params
    config_dict["brdf"]['bin_type'] = 'dynamic'
    config_dict["brdf"]['num_bins'] = 18
    config_dict["brdf"]['ndvi_bin_min'] = 0.05
    config_dict["brdf"]['ndvi_bin_max'] = 1.0
    config_dict["brdf"]['ndvi_perc_min'] = 10
    config_dict["brdf"]['ndvi_perc_max'] = 95
   
    #config_dict["resample"]  = False #ais resampling really only important when fitting to many sites at once
    config_dict["resample"] =True
    config_dict["resampler"] = {}
    config_dict["resampler"]['type'] = 'cubic'
    config_dict["resampler"]['out_waves'] = [x for x in range(400,2500,5)]  ################ change this as needed #################
    config_dict["resampler"]['out_fwhm'] = []
    config_dict['num_cpus'] = len(images)
   
    with open(config_file, 'w') as outfile:
       json.dump(config_dict,outfile,indent=3)

    print('finished create_config')

################################ Part 3: Perform BRDF image correction (No change needed) ###########################

def export_coeffs(hy_obj,export_dict):
       '''Export correction coefficients to file.
       '''
       for correction in hy_obj.corrections:
           coeff_file = export_dict['output_dir']
           coeff_file += os.path.splitext(os.path.basename(hy_obj.file_name))[0]
           coeff_file += "_%s_coeffs_%s.json" % (correction,export_dict["suffix"])
   
           with open(coeff_file, 'w') as outfile:
               if correction == 'topo':
                   corr_dict = hy_obj.topo
               elif correction == 'glint':
                   continue
               else:
                   corr_dict = hy_obj.brdf
               json.dump(corr_dict,outfile)
   
def apply_corrections(hy_obj,config_dict):
       '''Apply correction to image and export
           to file.
       '''
   
       header_dict = hy_obj.get_header()
       header_dict['data ignore value'] = hy_obj.no_data
       header_dict['data type'] = 2 
   
       output_name = config_dict['export']['output_dir']
       output_name += os.path.splitext(os.path.basename(hy_obj.file_name))[0]
       output_name +=  "_%s" % config_dict['export']["suffix"]
   
       #Export all wavelengths
       if len(config_dict['export']['subset_waves']) == 0:
   
           if config_dict["resample"] == True:
               hy_obj.resampler = config_dict['resampler']
               waves= hy_obj.resampler['out_waves']
           else:
               waves = hy_obj.wavelengths
   
           header_dict['bands'] = len(waves)
           header_dict['wavelength'] = waves
   
           writer = WriteENVI(output_name,header_dict)
           iterator = hy_obj.iterate(by='line', corrections=hy_obj.corrections,
                                     resample=config_dict['resample'])
           while not iterator.complete:
               line = iterator.read_next()
               writer.write_line(line,iterator.current_line)
           writer.close()
   
       #Export subset of wavelengths
       else:
           waves = config_dict['export']['subset_waves']
           bands = [hy_obj.wave_to_band(x) for x in waves]
           waves = [round(hy_obj.wavelengths[x],2) for x in bands]
           header_dict['bands'] = len(bands)
           header_dict['wavelength'] = waves
   
           writer = WriteENVI(output_name,header_dict)
           for b,band_num in enumerate(bands):
               band = hy_obj.get_band(band_num,
                                      corrections=hy_obj.corrections)
               writer.write_band(band, b)
           writer.close()
   
       #Export masks
       if (config_dict['export']['masks']) and (len(config_dict["corrections"]) > 0):
           masks = []
           mask_names = []
   
           for correction in config_dict["corrections"]:
               for mask_type in config_dict[correction]['apply_mask']:
                   mask_names.append(correction + '_' + mask_type[0])
                   masks.append(mask_create(hy_obj, [mask_type]))
   
           header_dict['data type'] = 1
           header_dict['bands'] = len(masks)
           header_dict['band names'] = mask_names
           header_dict['samples'] = hy_obj.columns
           header_dict['lines'] = hy_obj.lines
           header_dict['wavelength'] = []
           header_dict['fwhm'] = []
           header_dict['wavelength units'] = ''
           header_dict['data ignore value'] = 255
   
   
           output_name = config_dict['export']['output_dir']
           output_name += os.path.splitext(os.path.basename(hy_obj.file_name))[0]
           output_name +=  "_%s_mask" % config_dict['export']["suffix"]
   
           writer = WriteENVI(output_name,header_dict)
   
           for band_num,mask in enumerate(masks):
               mask = mask.astype(int)
               mask[~hy_obj.mask['no_data']] = 255
               writer.write_band(mask,band_num)
   
           del masks

def implement_brdf_correction(path):
   
   print('begin implement_brdf_correction')
   start=time.time()
   warnings.filterwarnings("ignore")
   np.seterr(divide='ignore', invalid='ignore')
    
   config_file = path + "config.json"

   with open(config_file, 'r') as outfile:
       config_dict = json.load(outfile)

   images = config_dict["input_files"]

   if ray.is_initialized():
       ray.shutdown()
   print("Using %s CPUs." % config_dict['num_cpus'])
   ray.init(num_cpus = config_dict['num_cpus'])

   HyTools = ray.remote(ht.HyTools)
   actors = [HyTools.remote() for image in images]

   if config_dict['file_type'] == 'envi':
       anc_files = config_dict["anc_files"]
       _ = ray.get([a.read_file.remote(image,config_dict['file_type'],
                                       anc_files[image]) for a,image in zip(actors,images)])

   elif config_dict['file_type'] == 'neon':
       _ = ray.get([a.read_file.remote(image,config_dict['file_type']) for a,image in zip(actors,images)])

   _ = ray.get([a.create_bad_bands.remote(config_dict['bad_bands']) for a in actors])

   for correction in config_dict["corrections"]:
       if correction =='topo':
           calc_topo_coeffs(actors,config_dict['topo'])
       elif correction == 'brdf':
           calc_brdf_coeffs(actors,config_dict)
       elif correction == 'glint':
           set_glint_parameters(actors,config_dict)

   if config_dict['export']['coeffs'] and len(config_dict["corrections"]) > 0:
       print("Exporting correction coefficients.")
       _ = ray.get([a.do.remote(export_coeffs,config_dict['export']) for a in actors])

   if config_dict['export']['image']:
       print("Exporting corrected images.")
       _ = ray.get([a.do.remote(apply_corrections,config_dict) for a in actors])

   ray.shutdown()
   end = time.time()
   total_time = end - start
   print("\n" + str(total_time)) 



def convert_envi_to_tif(site, year, data_int_path, envi_path):
    # from Shashi Konduri
    # takes about 1min to convert each flightline from envi to tif
    print("converting envi to tif")
    tif_path = os.path.join(data_int_path, site, year, 'hs_tif_flightline/')
    
    pattern = "*_BRDF_topo_corrected"

    for envi_file in glob.glob(pattern):

        # check if that tif file has been generated
        if os.path.exists(tif_path + envi_file + '.tif'):

            print(envi_file + ' tiff already exists in ' + envi_path)

        else:
            
            print("Working on %s" %envi_file)
            
            ## read the envi file using rasterio
            src = rasterio.open(envi_file, driver="ENVI")
            
            ## Set the name of the output tif file. This will save the file in the "location" specified above
            output_name = envi_path + envi_file + ".tif"
            
            ## This will copy all raster info from the ENVI file like height, width, coordinate reference system etc. into the new tif file 
            output = rasterio.open(output_name, 'w', driver="GTiff", height=src.height, width=src.width, count=src.count+1, #+1 for zen layer
                                dtype='float32', crs=src.crs, transform=src.transform, nodata=src.nodata)
            
            ## Save each of the 426 bands in the new tif file.
            ## src.read(i) reads the ith band from ENVI file
            for i in range(1,src.count+1):
                output.write(src.read(i),i)
                print(i)

            # Add another layer with to-sensor zenith angle, or should I just save it as a separate tiff file?
            # could probably do this with envi files, but not sure how... 
            # ais reflectance data stored in corrected envi file, and sensor_zn stored in anc envi file
            h5_path = envi_path + envi_file.split(pattern[1:])[0] + '.h5' 
            f = h5py.File(h5_path, 'r') # load in h5 from the directory
            zen = f['SOAP']['Reflectance']['Metadata']['to-sensor_Zenith_Angle'][:] # extract zenith angle
            output.write(zen, src.count+1) # save as an extra layer in tiff

            output.close()



# ########################## Part 1: Convert NEON HDFs to ENVI files (NO CHANGE NEEDED) #########################

# def move_red_yellow_subset(site, path): #ais adapted function from https://stackoverflow.com/questions/61963546/how-to-conditionally-move-files-from-one-directory-to-another
#     dest = os.path.join(path, 'cloud_condition_red_yellow')
#     if not os.path.exists(dest):
#         os.makedirs(dest)
#     count = 0
#     for file in os.listdir(path):
#         if file.endswith('.h5'):
#             f = h5py.File(os.path.join(path, file), 'r')
#             if f[site]['Reflectance']['Reflectance_Data'].attrs['Cloud conditions'] != 'Green (<10%) cloud cover':
#                 target = os.path.join(dest, file)
#                 if os.path.exists(target):
#                     print('file ' + file + ' already moved')
#                 shutil.move(os.path.join(path, file), os.path.join(dest, file))
#                 print(file, f[site]['Reflectance']['Reflectance_Data'].attrs['Cloud conditions'],
#                 f[site]['Reflectance']['Reflectance_Data'].attrs['Cloud type'])
#                 count += 1
#             print('moved: %d' % count)
#     print('Red/yellow file movement finished')

# def convert_hdf_to_envi(site, year, flightline_h5_path, data_int_path):
#     '''This function exports NEON AOP HDF imaging spectroscopy data
#     to an ENVI formated binary file, with the option of also exporting
#     ancillary data following formatting used by NASA JPL for AVIRIS
#     observables. The script utilizes ray to export images in parralel.
#     '''
#     print('begin convert_hdf_to_envi')
#     images = glob.glob(flightline_h5_path + "*.h5") 
#     output_dir = os.path.join(data_int_path, site, year, 'hs_envi_flightline/')
#     if not os.path.exists(output_dir):
#         os.makedirs(output_dir)

#     if ray.is_initialized():
#         ray.shutdown()
#     ray.init(num_cpus = len(images))

#     hytool = ray.remote(ht.HyTools)
#     actors = [hytool.remote() for image in images]
#     _ = ray.get([a.read_file.remote(image,'neon') for a,image in zip(actors,images)])

#     def neon_to_envi(hy_obj):
#         basemame = os.path.basename(os.path.splitext(hy_obj.file_name)[0])
#         print("Exporting %s " % basemame)
#         output_name = output_dir + basemame
#         writer = WriteENVI(output_name,hy_obj.get_header())
#         iterator = hy_obj.iterate(by = 'chunk')
#         pixels_processed = 0
#         while not iterator.complete:
#             chunk = iterator.read_next()
#             pixels_processed += chunk.shape[0]*chunk.shape[1]
#             writer.write_chunk(chunk,iterator.current_line,iterator.current_column)
#             if iterator.complete:
#                 writer.close()

#     def export_anc(hy_obj):
#         anc_header = hy_obj.get_header()
#         anc_header['bands'] = 10
#         anc_header['band_names'] = ['path length', 'to-sensor azimuth',
#                                     'to-sensor zenith','to-sun azimuth',
#                                       'to-sun zenith','phase', 'slope',
#                                       'aspect', 'cosine i','UTC time']
#         anc_header['wavelength units'] = np.nan
#         anc_header['wavelength'] = np.nan
#         anc_header['data type'] = 4

#         output_name = output_dir + os.path.basename(os.path.splitext(hy_obj.file_name)[0])
#         writer = WriteENVI(output_name + "_ancillary", anc_header)
#         writer.write_band(hy_obj.get_anc("path_length"),0)
#         writer.write_band(hy_obj.get_anc("sensor_az",radians = False),1)
#         writer.write_band(hy_obj.get_anc("sensor_zn",radians = False),2)
#         writer.write_band(hy_obj.get_anc("solar_az",radians = False),3)
#         writer.write_band(hy_obj.get_anc("solar_zn",radians = False),4)
#         #writer.write_band(hy_obj.get_anc("phase placeholder"),5)
#         writer.write_band(hy_obj.get_anc("slope",radians = False),6)
#         writer.write_band(hy_obj.get_anc("aspect",radians = False),7)
#         writer.write_band(hy_obj.cosine_i(),8)
#         #writer.write_band('UTC time placeholder',9)
#         writer.close()
    
#     _ = ray.get([a.do.remote(neon_to_envi) for a in actors])
    
#     print("\nExporting ancillary data")
#     _ = ray.get([a.do.remote(export_anc) for a in actors])

#     return(output_dir)
#     print("Export convert_hdf_to_envi complete.") #ais takes 20-25min for SOAP

# ########################## Part 2: Create config file with all the parameters for BRDF correction (TWEAK AS NECESSARY) #################################

# def create_config(input_path, output_path):
#     '''Template script for generating image_correct configuration JSON files.
#        These settiing are meant only as an example, are not appropriate for
#        all situations and may need to be adjusted
#     '''
#     print('begin create_config')

#    # Output path for configuration file
#     config_file = output_path + "config.json"
#     config_dict = {}
   
#    # Only coefficients for good bands will be calculated
#     config_dict['bad_bands'] = [] #[[300,400],[1337,1430],[1800,1960],[2450,2600]]   ####### Change this as needed #######
#    #ais make blank so that zenith tif will be 427th layer

#    # Input data settings for ENVI
#     config_dict['file_type'] = 'envi'
#     aviris_anc_names = ['path_length','sensor_az','sensor_zn',
#                         'solar_az', 'solar_zn','phase','slope',
#                         'aspect', 'cosine_i','utc_time']
#     images= glob.glob(input_path + "NEON*reflectance")
#     images.sort()
#     config_dict["input_files"] = images
#     config_dict["anc_files"] = {}
#     anc_files = glob.glob(output_path + "*ancillary")
#     anc_files.sort()
#     for i,image in enumerate(images):
#         config_dict["anc_files"][image] = dict(zip(aviris_anc_names,
#                                                                [[anc_files[i],a] for a in range(len(aviris_anc_names))]))
   
#    # Export settings
#     ''' Options for subset waves:
#        1. List of subset wavelenths
#        2. Empty list, this will output all good bands, if resampler is
#        set it will also resample.
#        - Currently resampler cannot be used in conjuction with option 1
#     '''
   
#     config_dict['export'] = {}
#     config_dict['export']['coeffs']  = False
#     config_dict['export']['image']  = True
#     config_dict['export']['masks']  = False
#     config_dict['export']['subset_waves']  = []
#     config_dict['export']['output_dir'] = output_path
#     config_dict['export']["suffix"] = 'BRDF_topo_corrected'
   
#    #Corrections
#     ''' Specify correction(s) to be applied, corrections will be applied
#    in the order they are specified.
#    Options include:
#        ['topo']
#        ['brdf']
#        ['glint']
#        ['topo','brdf']
#        ['brdf','topo']
#        ['brdf','topo','glint']
#        [] <---Export uncorrected images
#     '''
   
#     config_dict["corrections"] = ['brdf','topo','glint']
   
#     #Topographic Correction options
#     '''
#     Types supported:
#         - 'cosine'
#         - 'c'
#         - 'scs
#         - 'scs+c'
#         - 'mod_minneart'
#         - 'precomputed'
#     Apply and calc masks are only needed for C and SCS+C corrections. They will
#     be ignored in all other cases and correction will be applied to all
#     non no-data pixels.
#     'c_fit_type' is only applicable for the C or SCS+C correction type. Options
#     include 'ols' or 'nnls'. Choosing 'nnls' can limit overcorrection.
#     For precomputed topographic coefficients 'coeff_files' is a
#     dictionary where each key is the full the image path and value
#     is the full path to coefficients file, one per image.
#     '''
#     config_dict["topo"] =  {}
#     config_dict["topo"]['type'] =  'scs+c'
#     config_dict["topo"]['calc_mask'] = [["ndi", {'band_1': 850,'band_2': 660,
#                                                 'min': 0.1,'max': 1.0}]]
   
#     config_dict["topo"]['apply_mask'] = [["ndi", {'band_1': 850,'band_2': 660,
#                                                 'min': 0.1,'max': 1.0}]] #could also add something about cloud_mask - check hytools GitHub, 
#     config_dict["topo"]['c_fit_type'] = 'nnls'
   
#     #BRDF Correction options
#     '''
#     Types supported:
#        - 'universal': Simple kernel multiplicative correction.
#        - 'local': Correction by class. (Future.....)
#        - 'flex' : Correction by NDVI class
#        - 'precomputed' : Use precomputed coefficients
#     If 'bin_type' == 'user'
#     'bins' should be a list of lists, each list the NDVI bounds [low,high]
#     Object shapes ('h/b','b/r') only needed for Li kernels.
#     For precomputed topographic coefficients 'coeff_files' is a
#     dictionary where each key is the full the image path and value
#     is the full path to coefficients file, one per image.
#     '''
   
#     config_dict["brdf"]  = {}
    
#     # Options are 'line','scene', or a float for a custom solar zn
#     # Custom solar zenith angle should be in radians
#     config_dict["brdf"]['solar_zn_type'] ='scene'
    
#     # Flex BRDF configs
#     config_dict["brdf"]['type'] =  'flex'
#     config_dict["brdf"]['grouped'] =  True
#     config_dict["brdf"]['geometric'] = 'li_sparse'
#     config_dict["brdf"]['volume'] = 'ross_thick'
#     config_dict["brdf"]["b/r"] = 10
#     config_dict["brdf"]["h/b"] = 2
#     config_dict["brdf"]['sample_perc'] = 0.1
#     config_dict["brdf"]['interp_kind'] = 'linear'
#     config_dict["brdf"]['calc_mask'] = [["ndi", {'band_1': 850,'band_2': 660,
#                                                     'min': 0.1,'max': 1.0}]]
#     config_dict["brdf"]['apply_mask'] = [["ndi", {'band_1': 850,'band_2': 660,
#                                                  'min': 0.1,'max': 1.0}]]
#    # Flex dynamic NDVI params
#     config_dict["brdf"]['bin_type'] = 'dynamic'
#     config_dict["brdf"]['num_bins'] = 18
#     config_dict["brdf"]['ndvi_bin_min'] = 0.05
#     config_dict["brdf"]['ndvi_bin_max'] = 1.0
#     config_dict["brdf"]['ndvi_perc_min'] = 10
#     config_dict["brdf"]['ndvi_perc_max'] = 95
   
#     #config_dict["resample"]  = False #ais resampling really only important when fitting to many sites at once
#     config_dict["resample"] =True
#     config_dict["resampler"] = {}
#     config_dict["resampler"]['type'] = 'cubic'
#     config_dict["resampler"]['out_waves'] = [x for x in range(400,2500,5)]                  ################ change this as needed #################
#     config_dict["resampler"]['out_fwhm'] = []
#     config_dict['num_cpus'] = len(images)
   
#     with open(config_file, 'w') as outfile:
#        json.dump(config_dict,outfile,indent=3)

#     print('finished create_config')

# ################################ Part 3: Perform BRDF image correction (No change needed) ###########################

# def export_coeffs(hy_obj,export_dict):
#        '''Export correction coefficients to file.
#        '''
#        for correction in hy_obj.corrections:
#            coeff_file = export_dict['output_dir']
#            coeff_file += os.path.splitext(os.path.basename(hy_obj.file_name))[0]
#            coeff_file += "_%s_coeffs_%s.json" % (correction,export_dict["suffix"])
   
#            with open(coeff_file, 'w') as outfile:
#                if correction == 'topo':
#                    corr_dict = hy_obj.topo
#                elif correction == 'glint':
#                    continue
#                else:
#                    corr_dict = hy_obj.brdf
#                json.dump(corr_dict,outfile)
   
# def apply_corrections(hy_obj,config_dict):
#        '''Apply correction to image and export
#            to file.
#        '''
   
#        header_dict = hy_obj.get_header()
#        header_dict['data ignore value'] = hy_obj.no_data
#        header_dict['data type'] = 2 
   
#        output_name = config_dict['export']['output_dir']
#        output_name += os.path.splitext(os.path.basename(hy_obj.file_name))[0]
#        output_name +=  "_%s" % config_dict['export']["suffix"]
   
#        #Export all wavelengths
#        if len(config_dict['export']['subset_waves']) == 0:
   
#            if config_dict["resample"] == True:
#                hy_obj.resampler = config_dict['resampler']
#                waves= hy_obj.resampler['out_waves']
#            else:
#                waves = hy_obj.wavelengths
   
#            header_dict['bands'] = len(waves)
#            header_dict['wavelength'] = waves
   
#            writer = WriteENVI(output_name,header_dict)
#            iterator = hy_obj.iterate(by='line', corrections=hy_obj.corrections,
#                                      resample=config_dict['resample'])
#            while not iterator.complete:
#                line = iterator.read_next()
#                writer.write_line(line,iterator.current_line)
#            writer.close()
   
#        #Export subset of wavelengths
#        else:
#            waves = config_dict['export']['subset_waves']
#            bands = [hy_obj.wave_to_band(x) for x in waves]
#            waves = [round(hy_obj.wavelengths[x],2) for x in bands]
#            header_dict['bands'] = len(bands)
#            header_dict['wavelength'] = waves
   
#            writer = WriteENVI(output_name,header_dict)
#            for b,band_num in enumerate(bands):
#                band = hy_obj.get_band(band_num,
#                                       corrections=hy_obj.corrections)
#                writer.write_band(band, b)
#            writer.close()
   
#        #Export masks
#        if (config_dict['export']['masks']) and (len(config_dict["corrections"]) > 0):
#            masks = []
#            mask_names = []
   
#            for correction in config_dict["corrections"]:
#                for mask_type in config_dict[correction]['apply_mask']:
#                    mask_names.append(correction + '_' + mask_type[0])
#                    masks.append(mask_create(hy_obj, [mask_type]))
   
#            header_dict['data type'] = 1
#            header_dict['bands'] = len(masks)
#            header_dict['band names'] = mask_names
#            header_dict['samples'] = hy_obj.columns
#            header_dict['lines'] = hy_obj.lines
#            header_dict['wavelength'] = []
#            header_dict['fwhm'] = []
#            header_dict['wavelength units'] = ''
#            header_dict['data ignore value'] = 255
   
   
#            output_name = config_dict['export']['output_dir']
#            output_name += os.path.splitext(os.path.basename(hy_obj.file_name))[0]
#            output_name +=  "_%s_mask" % config_dict['export']["suffix"]
   
#            writer = WriteENVI(output_name,header_dict)
   
#            for band_num,mask in enumerate(masks):
#                mask = mask.astype(int)
#                mask[~hy_obj.mask['no_data']] = 255
#                writer.write_band(mask,band_num)
   
#            del masks

# def implement_brdf_correction(path):
   
#    print('begin implement_brdf_correction')
#    start=time.time()
#    warnings.filterwarnings("ignore")
#    np.seterr(divide='ignore', invalid='ignore')
    
#    config_file = path + "config.json"

#    with open(config_file, 'r') as outfile:
#        config_dict = json.load(outfile)

#    images = config_dict["input_files"]

#    if ray.is_initialized():
#        ray.shutdown()
#    print("Using %s CPUs." % config_dict['num_cpus'])
#    ray.init(num_cpus = config_dict['num_cpus'])

#    HyTools = ray.remote(ht.HyTools)
#    actors = [HyTools.remote() for image in images]

#    if config_dict['file_type'] == 'envi':
#        anc_files = config_dict["anc_files"]
#        _ = ray.get([a.read_file.remote(image,config_dict['file_type'],
#                                        anc_files[image]) for a,image in zip(actors,images)])

#    elif config_dict['file_type'] == 'neon':
#        _ = ray.get([a.read_file.remote(image,config_dict['file_type']) for a,image in zip(actors,images)])

#    _ = ray.get([a.create_bad_bands.remote(config_dict['bad_bands']) for a in actors])

#    for correction in config_dict["corrections"]:
#        if correction =='topo':
#            calc_topo_coeffs(actors,config_dict['topo'])
#        elif correction == 'brdf':
#            calc_brdf_coeffs(actors,config_dict)
#        elif correction == 'glint':
#            set_glint_parameters(actors,config_dict)

#    if config_dict['export']['coeffs'] and len(config_dict["corrections"]) > 0:
#        print("Exporting correction coefficients.")
#        _ = ray.get([a.do.remote(export_coeffs,config_dict['export']) for a in actors])

#    if config_dict['export']['image']:
#        print("Exporting corrected images.")
#        _ = ray.get([a.do.remote(apply_corrections,config_dict) for a in actors])

#    ray.shutdown()
#    end = time.time()
#    total_time = end - start
#    print("\n" + str(total_time)) 



# def convert_envi_to_tif(site, year, data_int_path, envi_path):
#     # from Shashi Konduri
#     # takes about 1min to convert each flightline from envi to tif

#     tif_path = os.path.join(data_int_path, site, year, 'hs_tif_flightline/')
    
#     pattern = "*_BRDF_topo_corrected"

#     for envi_file in glob.glob(pattern):

#         # check if that tif file has been generated
#         if os.path.exists(tif_path + envi_file + '.tif'):

#             print(envi_file + ' tiff already exists in ' + envi_path)

#         else:
            
#             print("Working on %s" %envi_file)
            
#             ## read the envi file using rasterio
#             src = rasterio.open(envi_file, driver="ENVI")
            
#             ## Set the name of the output tif file. This will save the file in the "location" specified above
#             output_name = envi_path + envi_file + ".tif"
            
#             ## This will copy all raster info from the ENVI file like height, width, coordinate reference system etc. into the new tif file 
#             output = rasterio.open(output_name, 'w', driver="GTiff", height=src.height, width=src.width, count=src.count+1, #+1 for zen layer
#                                 dtype='float32', crs=src.crs, transform=src.transform, nodata=src.nodata)
            
#             ## Save each of the 426 bands in the new tif file.
#             ## src.read(i) reads the ith band from ENVI file
#             for i in range(1,src.count+1):
#                 output.write(src.read(i),i)
#                 print(i)

#             # Add another layer with to-sensor zenith angle, or should I just save it as a separate tiff file?
#             # could probably do this with envi files, but not sure how... 
#             # ais reflectance data stored in corrected envi file, and sensor_zn stored in anc envi file
#             h5_path = envi_path + envi_file.split(pattern[1:])[0] + '.h5' 
#             f = h5py.File(h5_path, 'r') # load in h5 from the directory
#             zen = f['SOAP']['Reflectance']['Metadata']['to-sensor_Zenith_Angle'][:] # extract zenith angle
#             output.write(zen, src.count+1) # save as an extra layer in tiff

#             output.close()
