# Evan Hackstadt
# Whittall Lab Dudleya Genomics
# December 2025

# TODO: organize into functions and __main__ CLI

# NOTE: assumes this specific filename format:
#       popcode_LP_xxx_Du-xxx_Sxxx_Lxxx_Rx_00x.fastq.gz

"""
A streamlined way to update the sample samples.yaml file needed for Snakemake pipeline runs
Simply point the script at a folder containing your raw data of interest, and to the old config file
Useful for folks who don't want to use VIM or a remote VSCode session to edit the file
Script currenlty cannot edit ref_genome or anc_genome paths; those must be changed manually
"""

import sys
import os
import argparse
from pathlib import Path
import random
from datetime import datetime
import re

# --- CLI ARGS ---
parser = argparse.ArgumentParser()
parser.add_argument('data_dir', help='Path to the directory containing input data (raw reads) of interest. By default, adds all samples to a new config file.', type=str)
parser.add_argument('config_dir', help='Path to the directory where the samples.yaml file should be saved.', type=str)
parser.add_argument('-f', '--filename', help='(optional) string to be the name of the new config file (do not include .yaml extension)', type=str)
parser.add_argument('-n', '--n_samples', help='(optional) integer value --> add n random samples to config file', type=int)
parser.add_argument('-c', '--custom_samples', help='(optional) allows you to enter custom samples for the config file', action='store_true')
parser.add_argument('-x', '--exclude_samples', help='(optional) allows you to specify samples to exclude from the config file', action='store_true')
parser.add_argument('-a', '--append_to_file', help='(optional) provide path to an existing .yaml file to add selected samples to', type=str)
parser.add_argument('-q', '--quiet', help='(optional) script will not print modified config file contents after writing', action='store_true')
args = parser.parse_args()

# process args
data_dir = os.path.abspath(args.data_dir)   # ensure we have absolute paths
data_parent_dir = Path(data_dir).parent     # add to config so Snakemake can store intermediates next to raw reads
config_dir = os.path.abspath(args.config_dir)
filename = args.filename if args.filename else "samples"
n_samples = args.n_samples if args.n_samples else None

data_files = [f for f in os.listdir(data_dir) 
              if os.path.isfile(os.path.join(data_dir, f))]
data_files.sort()

# improper args handling
if not os.path.isdir(data_dir):
    raise ValueError(f"{data_dir} is not a directory.")
if not os.path.isdir(config_dir):
    raise ValueError(f"{config_dir} is not a directory.")
if args.n_samples is not None:
    if args.custom_samples is True:
        raise ValueError("Cannot use both --n_samples and --custom_samples flags. Choose one or neither.")
    if args.exclude_samples is True:
        raise ValueError("Cannot use both --n_samples and --exclude_samples flags. Choose one or neither.")
elif args.custom_samples is True and args.exclude_samples is True:
    raise ValueError("Cannot use both --custom_samples and --exclude_samples flags. Choose one or neither.")

if len(data_files) == 0:
    raise ValueError(f"{data_dir} does not contain any files.")
if len(data_files) % 2 != 0:
    raise ValueError(f"{data_dir} has an odd number of files. Expected even number (R1 and R2 for each sample).")
if not os.path.isfile(os.path.join(config_dir, "samples.yaml")):
    raise ValueError(f"No samples.yaml file found in {config_dir}. An initial file is required to update.")


# --- HELPER FUNCTION ---
def custom_selection(samples):
    selected_samples = []
    print(data_samples)
    print("Enter the exact name of a sample above, or simply its unique Du-### or LP_##")
    print("Type 'save' to finish selection. Type 'quit' to cancel and exit.")
    print(">>>>>>>>")
    user_input = input()
    print("<<<<<<<<")
    while user_input != 'quit' and user_input != 'save':
        if user_input in data_samples and not in selected_samples:  # if they entered the exact sample
            selected_samples.append(user_input)
            print(f"Added {user_input}")
            print("All selected samples: ", selected_samples)
        else:
            # attempt to search
            matches = [x for x in samples if user_input in x]
            if len(matches) == 1:
                if user_input not in selected_samples:
                    selected_samples.append(matches[0])
                    print(f"Found match {matches[0]} and added it.")
                    print("All selected samples: ", selected_samples)
                else:
                    print(f"Error: You've already selected '{user_input}'")
            elif len(matches) != 1:
                print(f"Error: found {len(matches)} matches for {user_input} in possible samples.")
                print(f"Valid samples are: {[x for x in data_samples if x not in selected_samples]}") if not args.quiet else print()
            else:
                print(f"Error: '{user_input}' is not a valid sample.")
                print(f"Valid samples are: {[x for x in data_samples if x not in selected_samples]}") if not args.quiet else print()
        print("Enter another sample, or 'save' to finish, or 'quit' to exit.")
        print(">>>>>>>>")
        user_input = input()
        print("<<<<<<<<")
    
    if user_input == 'quit':
        print("Cancelling script. No config file written.")
        sys.exit()
    if user_input == 'save':
        print("--------")
        print("Done selecting samples. Final selection:")
        print(selected_samples)
    
    return selected_samples


# --- MAIN SCRIPT ---

# from filenames, map R1 and R2 files to their sample names
data_samples = []
data_dict = {}
for f in data_files:
    if '_R1' in f:
        r1_path = os.path.abspath(os.path.join(data_dir, f))
        # filename parsing assuming: "popcode_LP_xxx_Du-xxx_Sxxx_Lxxx_Rx_00x.fastq.gz"
        DU_label = LP_num = ""
        # split on LP to get popcode (which may contain '_')
        popcode_substrings = f.split('_LP_')
        popcode = popcode_substrings[0]
        # extract Du# and LP# from split on '_'
        substrings = f.split('_')
        for i, s in enumerate(substrings):
            if 'Du-' in s:
                DU_label = s
            if s == 'LP':
                LP_num = substrings[i+1]   # once we find 'LP', the next substring is the LP number
        sample_name = f"{popcode}_LP_{LP_num}_{DU_label}"
        data_samples.append(sample_name)
        # now find the corresponding R2 file (check DU# and LP# for redundancy)
        for f2 in data_files:
            if DU_label in f2 and f'LP_{LP_num}' in f2 and '_R2' in f2:
                r2_path = os.path.abspath(os.path.join(data_dir, f2))
                data_dict[sample_name] = [r1_path, r2_path]
    elif '_R2' in f:
        continue    # we process R2 above, after finding its R1
    else:
        raise ValueError(f"File found without R1/R2 designation: {f}")
data_samples.sort()

# make sure nothing got messed up
if len(data_samples)*2 != len(data_files):
    raise ValueError(f"Error: found {len(data_samples)} but {len(data_files)}. Expected twice as many files as samples.")
for sample, file_list in data_dict.items():
    if len(file_list) != 2:
        raise ValueError(f"Error: mapped {len(file_list)} files to {sample}. Expected 2 (R1 and R2).")

print("\n----ARGS----")

print(f"Path to data directory: {data_dir}")

print(f"Found {len(data_files)} files in directory: ")
print(data_files) if not args.quiet else print("(message supressed)")

print(f"Parsed into {len(data_samples)} samples: ")
print(data_samples) if not args.quiet else print("(message supressed)")

print("Mapped filepaths to samples as follows:")
print(data_dict) if not args.quiet else print("(message supressed)")



# --- SELECT SAMPLES ---
print("\n----SAMPLE SELECTION----")
selected_samples = []

if args.n_samples is not None:      # n_samples
    print("Specific number of samples requested.")
    print(f"Selecting {n_samples} random from the directory.")
    if n_samples <= len(data_samples):
        selected_samples = random.sample(data_samples, n_samples)
    else:
        print(f"n={n_samples} provided, which is > number of samples present. ",
              f"Selecting all {len(data_samples)} samples instead.")
        selected_samples = data_samples

elif args.custom_samples:         # custom_samples
    print("Custom samples requested.")
    print("Please enter samples one-by-one to add to the config file.")
    selected_samples = custom_selection(data_samples)
    
elif args.exclude_samples:
    print("Exclusion of samples requested.")
    print("Please enter samples one-by-one to EXCLUDE from the config file.")
    excluded_samples = custom_selection(data_samples)
    print("Will add all samples except selected exclusions...")
    selected_samples = [x for x in data_samples if x not in excluded_samples]
    print("Final custom samples:")
    print(selected_samples)

else:               # (default) all samples
    print("All samples will be added to config file.")
    selected_samples = data_samples

print(f"Selected {len(selected_samples)} samples to write: ", selected_samples)



# --- WRITE TO CONFIG FILE ---
print("\n----WRITING TO FILE----")
old_config_path = os.path.join(config_dir, "samples.yaml")
config_path = os.path.join(config_dir, f"{filename}.yaml")
print("Config file path: ", config_path)

if args.append_to_file is not None:
    mode = 'a'
    config_path = os.path.abspath(args.append_to_file)
    print(f"Append requested. Selected samples will be added to end of existing config file: {config_path}")
else:
    mode = 'w'
    overwrite = True
    print("Mode is write.")
    
    # extract ref_genome and anc_genome so we can write them in the new file
    with open(old_config_path, 'r') as f:
        line = f.readline()
        while line == '\n' or line[0] == '#':
            line = f.readline()
        ref_line = line
        line = f.readline()
        while line == '\n' or line[0] == '#':
            line = f.readline()
        anc_line = line
        while line == '\n' or line[0] == '#':
            line = f.readline()
        chroms_line = line
        # ensure file was properly formatted...
        if 'ref_genome: ' not in ref_line:
            raise ValueError(f"Unable to extract ref_genome from old config file. Instead read: {ref_line}")
        if 'anc_genome: ' not in anc_line:
            raise ValueError(f"Unable to extract anc_genome from old config file. Instead read: {anc_line}")
        if 'chroms_file: ' not in chroms_line:
            raise ValueError(f"Unable to extract chroms_file from old config file. Instead read: {chroms_line}")
    print(f"...Extracted {ref_line}...Extracted {anc_line}...Extracted {chroms_line}...Parsed data_parent_dir: {data_parent_dir}\n")
    
    # optionally name this dataset for unique results folder
    output_name = "results"
    print("Would you like to give this dataset a unique name?")
    print("The name will be used to name the results folder, making it unique and preserving previous results folders.")
    print("Format: results-your_dataset_name")
    print("The default name is simply 'results'")
    print("Enter a single string as the name, or '0' to use default, or 'quit' to cancel the script.")
    print(">>>>>>>>")
    user_input = input()
    print("<<<<<<<<")
    while " " in user_input:
        print("Invalid input. No spaces. Please enter a string, '0', or 'quit'.")
        print(">>>>>>>>")
        user_input = input()
        print("<<<<<<<<")
    if user_input == '0':
        print("No custom dataset name. Output folder will simply be called 'results'")
    elif user_input == 'quit':
        print("Cancelling script. No config file written.")
        sys.exit()
    else:
        print(f"Name is {user_input}. Output folder will be results-{user_input}/")
        output_name = f"results-{user_input}"
    
    # if a samples.yaml already exists, ask whether to preserve or overwrite
    if os.path.isfile(os.path.join(config_dir, f"{filename}.yaml")) and args.append_to_file:
        print(f"WARNING! Current {filename}.yaml file exists. Would you like to preserve this file or overwrite it?")
        print("(1)\tPreserve existing config file")
        print("(2)\tOverwrite (delete) existing config file")
        print("(quit)\tCancel script")
        print(">>>>>>>>")
        user_input = input()
        print("<<<<<<<<")
        while user_input != '1' and user_input != '2' and user_input != 'quit':
            print("Invalid input. Please enter '1', '2', or 'quit'.")
            print(">>>>>>>>")
            user_input = input()
            print("<<<<<<<<")
        if user_input == '1':
            overwrite = False
            print(f"Preserve requested. Current {filename}.yaml file will be renamed.")
        elif user_input == '2':
            overwrite = True
            print(f"Overwrite requested. Current {filename}.yaml file will be replaced.")
        elif user_input == 'quit':
            print("Cancelling script. No config file written.")
            sys.exit()
    if not overwrite:
        timestamp = str(datetime.now().strftime("%Y-%m-%d_%H-%M-%S"))
        old_filename = os.path.join(config_dir, f"{filename}.yaml")
        renamed = os.path.join(config_dir, f"old_samples_{timestamp}.yaml")
        os.rename(old_filename, renamed)
        print(f"Renamed {filename}.yaml --> old_samples_{timestamp}.yaml. New file will be {filename}.yaml.")

with open(config_path, mode) as f:
    print(f"...Opening {config_path}...")

    # if we are writing a new file, add ref_genome and anc_genome to the top
    if mode == "w":
        f.write(f"# {config_path}")
        f.write(f"\n# Note: changing the layout of this file might break the update_sample_config.py script")
        f.write(f"\n{ref_line}")
        f.write(f"{anc_line}")
        f.write(f"{chroms_line}")
        f.write(f'\ndata_parent_dir: "{data_parent_dir}"')
        f.write(f'\noutput_name: "{output_name}"')
        f.write("\n\nsamples:")

    # now write the selected samples line-by-line
    for sample in selected_samples:
        # ensure sample is not already in file (skip duplicates)
        with open(config_path, 'r') as f2:
            contents = f2.read()
        if sample in contents:
            print(f">>>>WARNING: {sample} is already in config file. Skipping...")
        else:
            f.write(f'\n  {sample}:')
            f.write(f'\n    r1: "{data_dict[sample][0]}"')
            f.write(f'\n    r2: "{data_dict[sample][1]}"')

    # write a comment logging number of samples added for convenience
    f.write(f"\n# Added {len(selected_samples)} samples above")
print("Finished writing!")

if not args.quiet:
    print(f"\nNEW CONFIG FILE LOOKS LIKE:\n")
    with open(config_path, 'r') as f:
        print(f.read(), "\n")