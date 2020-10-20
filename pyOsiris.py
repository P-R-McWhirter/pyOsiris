import numpy as np
import argparse
import os
from astropy.io import fits

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

# Count the number of observation groups

input_folders = np.sort([i for i in os.listdir(input_path) if os.path.isdir(os.path.join(input_path, i))])

custom_std_grp = args['custom_std_grp']

input_grps = list(np.array(input_folders)[np.where(np.array(input_folders) != str(custom_std_grp))])

if custom_std_grp != None:
    print("Custom standard group defined: " + custom_std_grp)

for i in input_grps:
    print("Beginning group: " + i)

    arc = []
    bias = []
    flat = []
    obj = []
    stds = []

    for file in os.listdir(os.path.join(input_path, i, "arc")):
        if file.endswith(".fits"):
            arc.append(os.path.join(input_path, i, "arc", file))

    for file in os.listdir(os.path.join(input_path, i, "bias")):
        if file.endswith(".fits"):
            bias.append(os.path.join(input_path, i, "bias", file))

    for file in os.listdir(os.path.join(input_path, i, "flat")):
        if file.endswith(".fits"):
            flat.append(os.path.join(input_path, i, "flat", file))

    for file in os.listdir(os.path.join(input_path, i, "object")):
        if file.endswith(".fits"):
            obj.append(os.path.join(input_path, i, "object", file))

    for file in os.listdir(os.path.join(input_path, i, "stds")):
        if file.endswith(".fits"):
            stds.append(os.path.join(input_path, i, "stds", file))

    # Placeholder: list all the fits files in the expected folders:

    print(bias)

