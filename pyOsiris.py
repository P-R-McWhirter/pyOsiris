import numpy as np
import argparse
import os
from astropy import units as u
from astropy.io import fits
from astropy.nddata import CCDData
from astropy.stats import sigma_clip
from ccdproc import Combiner
from aspired import image_reduction
from aspired import spectral_reduction
from matplotlib.pyplot import *
ion()

ap = argparse.ArgumentParser()

# Get the folder containing the spectra in queue-mode form

ap.add_argument("--input_folder", type=str, required=True, help="Folder with the queue-mode format fits files.")
ap.add_argument("--custom_std_grp", type=str, required=False, help="Group name of a custom standard star observation.")



args = vars(ap.parse_args())

# Get the current work directory

cwd = os.getcwd()

# Fetch the input folder and then combine it into a full path

input_folder = args['input_folder']

input_path = os.path.join(cwd, input_folder)

# Define functions for the reduction

def make_master_bias(bias_frames, hdunum = 2, clip_low_bias = 5, clip_high_bias = 5):

    # Creates a master bias frame using median combine

    print("Creating a master bias frame...")

    bias_CCDData = []

    for i in bias_frames:
            # Open all the bias frames
            bias = fits.open(i)[hdunum]
            bias_CCDData.append(CCDData(bias.data, unit=u.adu))

    # Put data into a Combiner
    bias_combiner = Combiner(bias_CCDData)

    # Apply sigma clipping
    bias_combiner.sigma_clipping(low_thresh=clip_low_bias, high_thresh=clip_high_bias, func=np.ma.median)

    bias_master = bias_combiner.median_combine()

    # Free memory
    del bias_CCDData
    del bias_combiner

    print("Master bias frame complete.")

    return bias_master

def make_master_flat(flat_frames, master_bias, hdunum = 2, clip_low_flat = 5, clip_high_flat = 5):

    # Creates a master flat frame using median combine

    print("Creating a master flat frame...")

    flat_CCDData = []

    for i in flat_frames:
            # Open all the flat frames
            flat = fits.open(i)[hdunum]
            flat_data = CCDData(flat.data, unit=u.adu)
  
            # Subtract the master bias frame
            flat_data = flat_data.subtract(master_bias)

            # Add to the list
            flat_CCDData.append(flat_data)

    # Put data into a Combiner
    flat_combiner = Combiner(flat_CCDData)

    # Apply sigma clipping
    flat_combiner.sigma_clipping(low_thresh=clip_low_flat, high_thresh=clip_high_flat, func=np.ma.median)

    flat_master = flat_combiner.median_combine()

    # Free memory
    del flat_CCDData
    del flat_combiner

    print("Master flat frame complete.")

    return flat_master

def get_arc_data(arc_data, master_bias, master_flat, hdunum = 2):

    arc_CCDData = []

    arc_first = fits.open(arc_data[0])
    arc_header = arc_first[0].header
    arc_header_hdunum = arc_first[hdunum].header
    arc_header.extend(arc_header_hdunum)

    for i in arc_data:

        # Open all the standard frames

        arc_file = fits.open(i)[hdunum]
        arc_data = CCDData(arc_file.data, unit=u.adu)

        # Subtract the master bias frame
        arc_data = arc_data.subtract(master_bias)

        # Divide the master flat frame
        arc_data = arc_data.divide(master_flat)

        # Add to the list
        arc_CCDData.append(arc_data)

    # Put data into a combiner
    arc_combiner = Combiner(arc_CCDData)

    # Free memory
    del arc_CCDData

    arc_master = arc_combiner.median_combine()

    # Free memory
    del arc_combiner

    return arc_master, arc_header

def get_input_data(input_data, master_bias, master_flat, hdunum = 2, clip_low_input = 5, clip_high_input = 5):

    input_CCDData = []
    input_time = []
    input_header = []

    for i in input_data:

        # Open all the standard frames

        input_file = fits.open(i)
        input_data = CCDData(input_file[hdunum].data, unit=u.adu)

        # Subtract the master bias frame
        input_data = input_data.subtract(master_bias)

        # Divide the master flat frame
        input_data = input_data.divide(master_flat)

        # Add to the list
        input_CCDData.append(input_data)

        input_header_0 = input_file[0].header
        input_header_hdunum = input_file[hdunum].header
        input_header_0.extend(input_header_hdunum)
        input_header.append(input_header_0)

        # Get the exposure time

        input_time.append(input_header_0['EXPTIME'])

    # Put data into a combiner
    input_combiner = Combiner(input_CCDData)

    # Free memory
    del input_CCDData

    # Apply sigma clipping
    input_combiner.sigma_clipping(low_thresh=clip_low_input, high_thresh=clip_high_input, func=np.ma.median)

    input_master = input_combiner.median_combine()
    input_exptime = np.median(input_time)

    # Free memory
    del input_combiner

    return input_master, input_exptime, input_header

def create_fits(input_frame, input_header, light_filename, bias_filename, flat_filename, clip_low_bias, clip_high_bias, clip_low_flat, clip_high_flat, clip_low_light, clip_high_light, exptime_light, hdunum, rotate = True):

    # Rotate the frame if needed (used for std and obj usually)

    if rotate == True:

        input_frame = np.rot90(input_frame)

    input_frame = np.array((input_frame))

    # Construct a FITS object of the reduced frame

    input_fits = fits.ImageHDU(input_frame)

    input_fits.header = input_header

    # Add the names of all the light frames to header
    if len(light_filename) > 0:
        for i in range(len(light_filename)):
            input_fits.header.set(keyword='light' + str(i + 1),
                                  value=light_filename[i],
                                  comment='Light frames')

    # Add the names of all the bias frames to header
    if len(bias_filename) > 0:
        for i in range(len(bias_filename)):
            input_fits.header.set(keyword='bias' + str(i + 1),
                                  value=bias_filename[i],
                                  comment='Bias frames')

    # Add the names of all the flat frames to header
    if len(flat_filename) > 0:
        for i in range(len(flat_filename)):
            input_fits.header.set(keyword='flat' + str(i + 1),
                                  value=flat_filename[i],
                                  comment='Flat frames')

    # Add all the other keywords
    input_fits.header.set(
        keyword='COMBTYPE',
        value="median",
        comment='Type of image combine of the light frames.')
    input_fits.header.set(
        keyword='SIGCLIP',
        value="True",
        comment='True if the light frames are sigma clipped.')
    input_fits.header.set(
        keyword='CLIPLOW',
        value=clip_low_light,
        comment='Lower threshold of sigma clipping of the light frames.')
    input_fits.header.set(
        keyword='CLIPHIG',
        value=clip_high_light,
        comment='Higher threshold of sigma clipping of the light frames.')
    input_fits.header.set(
        keyword='XPOSURE',
        value=exptime_light,
        comment='Average exposure time of the light frames.')
    input_fits.header.set(
        keyword='KEYWORD',
        value="EXPTIME",
        comment='Automatically identified exposure time keyword of the '
        'light frames.')
    input_fits.header.set(
        keyword='BCOMTYPE',
        value="median",
        comment='Type of image combine of the bias frames.')
    input_fits.header.set(
        keyword='BSIGCLIP',
        value="True",
        comment='True if the dark frames are sigma clipped.')
    input_fits.header.set(
        keyword='BCLIPLOW',
        value=clip_low_bias,
        comment='Lower threshold of sigma clipping of the bias frames.')
    input_fits.header.set(
        keyword='BCLIPHIG',
        value=clip_high_bias,
        comment='Higher threshold of sigma clipping of the bias frames.')
    input_fits.header.set(
        keyword='FCOMTYPE',
        value="median",
        comment='Type of image combine of the flat frames.')
    input_fits.header.set(
        keyword='FSIGCLIP',
        value="True",
        comment='True if the flat frames are sigma clipped.')
    input_fits.header.set(
        keyword='FCLIPLOW',
        value=clip_low_flat,
        comment='Lower threshold of sigma clipping of the flat frames.')
    input_fits.header.set(
         keyword='FCLIPHIG',
         value=clip_high_flat,
         comment='Higher threshold of sigma clipping of the flat frames.')

    return input_frame

def save_fits(input_fits, output_folder, filename, extension = 'fits', overwrite = False):

    input_fits = fits.PrimaryHDU(input_fits)

    output_path = os.path.join(output_folder, filename + '.' + extension)

    # Save file to disk
    input_fits.writeto(output_path, overwrite=overwrite)

    return output_path

def remove_overscan():

    # Might not use this.

    return None

def create_onedspec():

    return None

def do_img_red(input_path, grp_path, arc, bias, flat, stds, obj, hdunum = 2, clip_low_bias = 5, clip_high_bias = 5, clip_low_flat = 5, clip_high_flat = 5, clip_low_light = 5, clip_high_light = 5):

    # Create a master bias frame.

    bias_master = make_master_bias(bias, hdunum = hdunum, clip_low_bias = clip_low_bias, clip_high_bias = clip_high_bias)

    # Create a master flat frame.

    flat_master = make_master_flat(flat, bias_master, hdunum = hdunum, clip_low_flat = clip_low_flat, clip_high_flat = clip_high_flat)

    # Read the HDUs of the arcs, standards and objects to determine number of grating combinations.

    big_fits_list = arc + stds + obj

    grisms = []

    for j in big_fits_list:

        hdul = fits.open(j)

        grisms.append(hdul[0].header['GRISM'])

    unique_grisms = np.unique(grisms)

    # Collect the frames appropriate for each grism.

    all_master_arc_paths = []
    all_master_stds_paths = []
    all_master_obj_paths = []

    for j in unique_grisms:

        if j == "OPEN":

            all_master_arc_paths.append(None)
            all_master_stds_paths.append(None)
            all_master_obj_paths.append(None)

            continue

        arc_grism = []

        stds_grism = []

        obj_grism = []

        for k in arc:

            hdul = fits.open(k)

            grism = hdul[0].header['GRISM']

            if grism == j:

                arc_grism.append(k)

        for k in stds:

            hdul = fits.open(k)

            grism = hdul[0].header['GRISM']

            if grism == j:

                stds_grism.append(k)

        for k in obj:

            hdul = fits.open(k)

            grism = hdul[0].header['GRISM']

            if grism == j:

                obj_grism.append(k)

        # Fetch the arc data

        print("Collecting arc frames for grism " + j + "...")

        if len(arc_grism) > 0:

            arc_master, arc_header = get_arc_data(arc_grism, bias_master, flat_master, hdunum = hdunum)

        else:
        
            arc_master = None
            arc_header = None

        print("Done.")

        # Fetch the std data

        print("Collecting standard frames for grism " + j + "...")

        if len(stds_grism) > 0:

            stds_master, stds_exptime, stds_header = get_input_data(stds_grism, bias_master, flat_master, hdunum = hdunum, clip_low_input = clip_low_light, clip_high_input = clip_high_light)

        else:

            stds_master = None
            stds_exptime = None
            stds_header = None

        print("Done.")

        # Fetch the object data

        print("Collecting light frames for grism " + j + "...")

        if len(obj_grism) > 0:

            obj_master, obj_exptime, obj_header = get_input_data(obj_grism, bias_master, flat_master, hdunum = hdunum, clip_low_input = clip_low_light, clip_high_input = clip_high_light)

        else:

            obj_master = None
            obj_exptime = None
            obj_header = None

        print("Done.")

        # Create fits objects and save them for the arcs, stds and obj frames.

        # First for the arc

        arc_fits = create_fits(arc_master, arc_header, obj_grism, bias, flat, clip_low_bias, clip_high_bias, clip_low_flat, clip_high_flat, clip_low_light, clip_high_light, obj_exptime, hdunum, rotate = True)

        # Then the std frame

        stds_fits = create_fits(stds_master, stds_header[0], stds_grism, bias, flat, clip_low_bias, clip_high_bias, clip_low_flat, clip_high_flat, clip_low_light, clip_high_light, stds_exptime, hdunum, rotate = True)

        # Then the obj frame

        obj_fits = create_fits(obj_master, obj_header[0], obj_grism, bias, flat, clip_low_bias, clip_high_bias, clip_low_flat, clip_high_flat, clip_low_light, clip_high_light, obj_exptime, hdunum, rotate = True)

        # Save the master arc, master std and master obj frame to the base group directory.

        master_arc_path = save_fits(arc_fits, os.path.join(input_path, grp_path), "arc_" + j + "_master", extension = 'fits', overwrite = True)

        master_stds_path = save_fits(stds_fits, os.path.join(input_path, grp_path), "stds_" + j + "_master", extension = 'fits', overwrite = True)

        master_obj_path = save_fits(obj_fits, os.path.join(input_path, grp_path), "obj_" + j + "_master", extension = 'fits', overwrite = True)

        all_master_arc_paths.append(master_arc_path)
        all_master_stds_paths.append(master_stds_path)
        all_master_obj_paths.append(master_obj_path)

    return unique_grisms, all_master_arc_paths, all_master_stds_paths, all_master_obj_paths

def do_spec_red(input_path, grp_path, unique_grisms, all_master_arc_paths, all_master_stds_paths, all_master_obj_paths, custom_std_grp):

    for n, j in enumerate(unique_grisms):

        if j == "OPEN":

            continue

        spec_science = fits.open(all_master_obj_paths[n])
        spec_standard = fits.open(all_master_stds_paths[n])
        spec_arc = fits.open(all_master_arc_paths[n])

        spec_science[0].data[np.isnan(spec_science[0].data)] = 0.
        spec_standard[0].data[np.isnan(spec_standard[0].data)] = 0.
        spec_arc[0].data[np.isnan(spec_arc[0].data)] = 0.

        spatial_mask = np.arange(100, 1000)

        science_twodspec = spectral_reduction.TwoDSpec(
            spec_science[0].data,
            spec_science[0].header,
            spatial_mask=spatial_mask,
            saxis=1)

        standard_twodspec = spectral_reduction.TwoDSpec(
            spec_standard[0].data,
            spec_standard[0].header,
            spatial_mask=spatial_mask,
            saxis=1)

        science_twodspec.ap_trace(nspec=1)
        science_twodspec.ap_extract(display=True, apwidth = 15)

        standard_twodspec.ap_trace(nspec=1)
        standard_twodspec.ap_extract(display=True, apwidth = 15)

        onedspec = spectral_reduction.OneDSpec()
        onedspec.from_twodspec(science_twodspec, stype='science')
        onedspec.from_twodspec(standard_twodspec, stype='standard')
        onedspec.add_arc(spec_arc[0].data[spatial_mask])
        onedspec.extract_arc_spec()
        onedspec.find_arc_lines(background=1e-3, prominence=1e-2)

        try:

            atlas = [
                3650.153, 4046.563, 4077.831, 4358.328, 5460.735, 5769.598, 5790.663,
                6682.960, 6752.834, 6871.289, 6965.431, 7030.251, 7067.218, 7147.042,
                7272.936, 7383.981, 7503.869, 7514.652, 7635.106, 7723.98
            ]
            element = ['HgAr'] * len(atlas)

            onedspec.initialise_calibrator(  stype='science+standard')
            onedspec.set_hough_properties(
                num_slopes=10000,
                xbins=200,
                ybins=200,
                min_wavelength=3500.,
                max_wavelength=8000.,
                range_tolerance=500.,
                linearity_tolerance=50,
                stype='science+standard')

            onedspec.load_user_atlas(
                elements=element,
                wavelengths=atlas,
                stype='science+standard')

            onedspec.add_atlas(
                elements=['Ne'],
                min_intensity=5.,
                min_atlas_wavelength=3500.,
                max_atlas_wavelength=8000.)

            onedspec.set_ransac_properties(
                sample_size=5,
                top_n_candidate=5,
                linear=True,
                filter_close=True,
                ransac_tolerance=5,
                candidate_weighted=True,
                hough_weight=1.0,
                stype='science+standard')
            onedspec.do_hough_transform()
            onedspec.fit(max_tries=50, stype='science+standard')
            onedspec.apply_wavelength_calibration(stype='science+standard')

            onedspec.load_standard(target='Feige110', library='esoxshooter')

            onedspec.inspect_standard(save_iframe=True, filename='test/test_inspect_standard_' + j)

            onedspec.compute_sensitivity(kind='cubic', mask_fit_size = 1)
            onedspec.inspect_sensitivity(save_iframe=True, filename='test/test_sensitivity_' + j)

            onedspec.apply_flux_calibration(stype='science+standard')

            onedspec.inspect_reduced_spectrum(stype='science', save_iframe=True, filename='test/test_science_spectrum_' + j)

            onedspec.inspect_reduced_spectrum(stype='standard', save_iframe=True, filename='test/test_standard_spectrum_' + j)

            onedspec.inspect_reduced_spectrum(display=True)

        except:

            print("Ransac failed (probably R2500U) PLACEHOLDER ERROR")

            continue

    return None

#############################################



# Count the number of observation groups

input_folders = np.sort([i for i in os.listdir(input_path) if os.path.isdir(os.path.join(input_path, i))])

custom_std_grp = args['custom_std_grp']

input_grps = list(np.array(input_folders)[np.where(np.array(input_folders) != str(custom_std_grp))])

if custom_std_grp != None:
    print("Custom standard group defined: " + custom_std_grp)

    custom_std_present = os.path.isdir(os.path.join(input_path, custom_std_grp))

    if custom_std_present == False:
        raise NameError("Custom standard group not found in input folder. Exiting...")
    else:
        print("Using custom standard group...")

for i in input_folders:
    print("Beginning group: " + i)

    arc = []
    bias = []
    flat = []
    stds = []
    obj = []
    
    for file in os.listdir(os.path.join(input_path, i, "arc")):
        if file.endswith(".fits"):
            arc.append(os.path.join(input_path, i, "arc", file))

    for file in os.listdir(os.path.join(input_path, i, "bias")):
        if file.endswith(".fits"):
            bias.append(os.path.join(input_path, i, "bias", file))

    for file in os.listdir(os.path.join(input_path, i, "flat")):
        if file.endswith(".fits"):
            flat.append(os.path.join(input_path, i, "flat", file))

    for file in os.listdir(os.path.join(input_path, i, "stds")):
        if file.endswith(".fits"):
            stds.append(os.path.join(input_path, i, "stds", file))

    for file in os.listdir(os.path.join(input_path, i, "object")):
        if file.endswith(".fits"):
            obj.append(os.path.join(input_path, i, "object", file))

    # Run image reduction routine

    unique_grisms, all_master_arc_paths, all_master_stds_paths, all_master_obj_paths = do_img_red(input_path, i, arc, bias, flat, stds, obj, hdunum = 2, clip_low_bias = 5, clip_high_bias = 5, clip_low_flat = 5, clip_high_flat = 5, clip_low_light = 5, clip_high_light = 5)

    # Execute the Aspired 2d spec reduction

    do_spec_red(input_path, i, unique_grisms, all_master_arc_paths, all_master_stds_paths, all_master_obj_paths, custom_std_grp=None)



