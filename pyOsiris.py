import numpy as np
import argparse
import os
from astropy import units as u
from astropy.io import fits
from astropy.nddata import CCDData
from astropy.stats import sigma_clip
from ccdproc import Combiner
from aspired import spectral_reduction

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

def make_master_flat(flat_frames, master_bias, hdunum = 2, clip_low_bias = 5, clip_high_bias = 5):

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
    flat_combiner.sigma_clipping(low_thresh=clip_low_bias, high_thresh=clip_high_bias, func=np.ma.median)

    flat_master = flat_combiner.median_combine()

    # Free memory
    del flat_CCDData
    del flat_combiner

    print("Master flat frame complete.")

    return flat_master

def get_arc_data(arc_data, master_bias, master_flat, hdunum = 2):

    arc_CCDData = []

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

    return arc_master

def get_input_data(input_data, master_bias, master_flat, hdunum = 2, clip_low_bias = 5, clip_high_bias = 5):

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
    input_combiner.sigma_clipping(low_thresh=clip_low_bias, high_thresh=clip_high_bias, func=np.ma.median)

    input_master = input_combiner.median_combine()
    input_exptime = np.median(input_time)

    # Free memory
    del input_combiner

    return input_master, input_exptime

def remove_overscan():

    # Might not use this.

    return None

def create_onedspec():

    return None

def do_img_red(input_path, grp_path, arc, bias, flat, stds, obj):

    # Create a master bias frame.

    bias_master = make_master_bias(bias, hdunum = 2, clip_low_bias = 5, clip_high_bias = 5)

    # Create a master flat frame.

    flat_master = make_master_flat(flat, bias_master, hdunum = 2, clip_low_bias = 5, clip_high_bias = 5)

    # Read the HDUs of the arcs, standards and objects to determine number of grating combinations.

    big_fits_list = arc + stds + obj

    grisms = []

    for j in big_fits_list:

        hdul = fits.open(j)

        grisms.append(hdul[0].header['GRISM'])

    unique_grisms = np.unique(grisms)

    # Collect the frames appropriate for each grism.

    for j in unique_grisms:

        if j == "OPEN":
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

            arc_master = get_arc_data(arc_grism, bias_master, flat_master, hdunum = 2)

        else:
        
            arc_master = None

        print("Done.")

        # Fetch the std data

        print("Collecting standard frames for grism " + j + "...")

        if len(stds_grism) > 0:

            stds_master, stds_exptime = get_input_data(stds_grism, bias_master, flat_master, hdunum = 2, clip_low_bias = 5, clip_high_bias = 5)

        else:

            stds_master = None
            stds_exptime = None

        print("Done.")

        # Fetch the object data

        print("Collecting light frames for grism " + j + "...")

        if len(obj_grism) > 0:

            obj_master, obj_exptime = get_input_data(obj_grism, bias_master, flat_master, hdunum = 2, clip_low_bias = 5, clip_high_bias = 5)

        else:

            obj_master = None
            obj_exptime = None

        print("Done.")

    return None





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

    do_img_red(input_path, i, arc, bias, flat, stds, obj)

