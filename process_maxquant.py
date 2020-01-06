#!/usr/bin/env python3

import numpy as np
import pandas as pd
import argparse
import sys
import operator
import re
import copy
from itertools import chain
import xlsxwriter
from natsort import natsorted
from collections import Counter

# shut off settingwithcopy warning
pd.options.mode.chained_assignment = None  # default='warn'

"""
Purpose: Pipeline for processing proteinGroups.txt files generated from MaxQuant.

This script is designed for processing MaxQuant output file, proteinGroups.txt, from triple SILAC mass spec experiments.
The proteinGroups.txt file can contain either one mix or two mixes with overlapping timepoints.
This script assumes that the data contains at least two replicates for each condition in the experimental and control
groups.

To use:
python process_maxquant.py -exp <path/to/file> -con <path/to/file> --prefix <string> --meta <path/to/file>
Example:
python process_maxquant.py -con Mix12_Con_txt/proteinGroups.txt -exp Mix12_Myd_txt/proteinGroups.txt --prefix Mix12_ConMyd --meta info.txt

Parameters:
Required:
-exp: proteinGroups.txt file for experimental group
-con: proteinGroups.txt file for control group
--meta: txt file with timepoint info
Optional:
--prefix: prefix for output file names

Files that are generated:
<prefix>_raw_filtered.xlsx: filtered, merged excel file with raw intensity values
<prefix>_counts.xlsx: excel file of replicate counts
<prefix>_exp_stats.xlsx: excel file with averages, standard deviations, and counts for experimental group
<prefix>_all_norm.xlsx: excel file with normalized values for each replicate for both control and experimental groups
<prefix>_all_norm_rand.xlsx: same as above but with random imputation for 0's
<prefix>_norm_avg.xlsx: excel file with normalized averages for both control and experimental groups

"""
__author__ = "Sara Jones"
__email__ = "jonessarae@gmail.com"
__doc__ = "Pipeline that processes files from MaxQuant."
__date__ = "1/3/2020"

def filter_maxquant_col(filename):
    """
    This function takes in MaxQuant file "proteinGroups.txt" and filters out
    columns not needed for analysis and returns a new dataframe.

    Arguments: filename (string)
    Returns: filtered_cols_df (dataframe)
    """
    # read in proteinGroups.txt
    df = pd.read_csv(filename, sep="\t", low_memory=False)

    # list of non-intensity columns to keep
    col_list = ["Majority protein IDs", "Protein names", "Gene names",
                "Fasta headers", "Peptides", "Unique peptides", "Sequence coverage [%]",
                "Mol. weight [kDa]","Sequence length", "Q-value", "Score", "Only identified by site",
                "Reverse", "Potential contaminant"]

    try:
        # create dataframe with non-intensity columns
        no_intensity_df = df[col_list]
    except KeyError:
        print("\nError")
        print("At least one of the files is not a proteinGroups.txt.")
        print("Please check you have the right files.")

    try:
        # create dataframe with intensity columns
        intensity_df = df[df.filter(regex="Intensity H |Intensity M |Intensity L ").columns]
    except KeyError:
        print("\nError")
        print("The proteinGroups.txt file must contain three isotopes.")
        print("Please check you have the right files.")

    # combine dataframes
    filtered_cols_df = pd.concat([no_intensity_df, intensity_df], axis=1)

    return filtered_cols_df

def filter_maxquant_row(df):
    """
    This function takes in a dataframe and filters out rows not needed for analysis.

    Arguments: df (dataframe)
    Returns: df (dataframe)
    """
    # remove rows with no unique peptides
    df = df[df["Unique peptides"] != 0]

    # remove rows with '+'
    df = df[df["Only identified by site"] != "+"]
    df = df[df["Reverse"] != "+"]
    df = df[df["Potential contaminant"] != "+"]

    # list of columns to remove
    col_list = ["Only identified by site", "Reverse", "Potential contaminant"]

    # drop columns
    df = df.drop(col_list, axis=1)

    return df

def sort_raw_dataframe(df, num):
    """
    This function takes in a dataframe and sorts the dataframe so that replicates
    for each stimulant/timepoint are grouped together.

    Arguments: df (dataframe)
               num (integer), number of unique replicates
    Returns: df (dataframe)
    """
    # get list of intensity column names
    intensity_col_list = list(df.columns[11:])

    # regular expression pattern for column name containing Intensity
    regex = r'Intensity\s([^\s])\s([^\s])[^\s]([0-9])([^\s])_([^\s])'

    # key for sorting columns containing Intensity by the following order:
    # stimulant, group (control vs experiment), mix, isotope, replicate
    key = [1,3,2,0,4]

    try:
        # sort list of intensity columns
        sorted_list = sorted(intensity_col_list, key=lambda name: tuple(re.findall(regex,name)[0][i] for i in key))
        # bug checker - will exit program if the following condition is not met
        assert len(sorted_list)%num==0,"Missing replicate for one of the conditions. Need {} replicates per condition.".format(num)
        # list of columns of protein info
        protein_list = list(df.columns[:11])
        # replace old columns with new columns
        df = df[protein_list+sorted_list]
    except IndexError as error:
        print("\nError")
        print("\nOne of the intensity column names does not have the following naming convention:")
        print("Example: Intensity H LM2C_B")
        sys.exit(0) # exit program
    return df

def create_count_df(df, num):
    """
    This function takes in a dataframe and creates a new dataframe counting the
    number of replicates with an intensity value greater than 0.

    If the count is less than or equal to one third of the number of replicates in the control group and
    equal to or greater than two thirds of the number of replicates in the experimental group,
    then the corresponding stimulant/timepoint is set to 1, else it is set to 0.

    For example, if the number of replicates is 3, a count of less than or equal to 1 in the control group
    and a count equal to or greater than 2 in the experimental group is set to 1.

    A column summing up every instance that a stimulant/timepoint is 1 is provided
    for filtering.

    Arguments: df (dataframe)
               num (integer), number of unique replicates
    Returns: count_df (dataframe)
    """
    # get list of intensity column names
    intensity_col_list = list(df.columns[11:])

    # create dataframe of intensity values
    intensity_df = df[intensity_col_list]

    # replace any NaN values with 0
    intensity_df[:].fillna(0, inplace=True)

    # create dataframe for counts with protein info
    count_df = df[df.columns[:11]]

    # loop through each stimulant/mix and count number of replicates with intensity
    # value greater than 0
    for i in range(len(intensity_df.columns))[::num]:
        # create new column name
        col_name_a = intensity_df.columns[i].split(" ")
        isotope_a = col_name_a[1] # isotope
        stim_mix_a = col_name_a[2].split("_")[0] # stimulant
        new_col_name_a = stim_mix_a + "_" + isotope_a # example: LMIC_H
        # create new column with count info
        count_df[new_col_name_a]=(intensity_df.iloc[:,i:i+num] != 0).astype(int).sum(axis=1)

    # get list of count column names
    count_col_list = list(count_df.columns[11:])
    # regular expression pattern for count column names
    regex = r'([^\s]*)([^\s])_([^\s])'
    # key for sorting count columns by stimulant/mix, isotope, control/experiment
    key = [0,2,1]
    # sort list of count columns
    new_list = sorted(count_col_list, key=lambda name: tuple(re.findall(regex,name)[0][i] for i in key))

    # bug checker - program will exit if the following condition is not met
    assert len(new_list)%2==0,"Missing replicate data for control and/or experimental group."

    # set count dataframe with sorted count columns and first 4 columns of protein info
    count_df = count_df[list(count_df.columns[:4]) + new_list]

    # get number of columns in count dataframe
    num_cols = len(count_df.columns)

    # thresholds for criteria
    low_thld = num*(1/3)
    high_thld = num*(2/3)

    # loop through each stimulant/mix and check if the number of replicates
    # meet criteria
    for i in range(4, num_cols)[::2]:
        # create new column name
        col_name_b = count_df.columns[i].split("C_")
        stim_mix_b = col_name_b[0] # stimulant
        isotope_b = col_name_b[1] # isotope
        new_col_name_b = stim_mix_b + "_" + isotope_b # example: PM1_M
        # create new column with value 1 if it meets condition
        count_df.loc[(count_df[count_df.columns[i]] <= low_thld) & (count_df[count_df.columns[i+1]] >= high_thld), new_col_name_b] = 1

    # loop through each column to replace any NaNs with 0
    for i in range(num_cols, len(count_df.columns)):
        count_df[count_df.columns[i]].fillna(0, inplace=True)

    # set numbers as integer type in condition columns
    count_df[list(count_df.columns[num_cols:len(count_df.columns)])] = count_df[list(count_df.columns[num_cols:len(count_df.columns)])].astype(int)

    # create new column to sum number of 1's in condition columns
    count_df["sum"] = (count_df.iloc[:,num_cols:len(count_df.columns)] != 0).astype(int).sum(axis=1)

    return count_df

def create_sample_dict(metafile):
    """
    This function takes in a tab-delmited file listing the timepoints to the isotopes/mix
    and returns a dictionary with the keys as isotopes/mix and values as timepoints.

    Arguments: metafile (text file)
    Returns: sample_dict (dictionary)
    """
    # create dictionary
    sample_dict = {}
    # read in file
    with open(metafile) as f:
        for line in f:
            # split line by tab delimiter
            (key, val) = line.split()
            # set key as istopes/mix and values as timepoints
            sample_dict[key] = int(val)

    return sample_dict

def get_common_isotope(sample_dict):
    """
    This function takes in a dictionary of isotopes/mix to timepoints and returns the
    the isotope that is common to both mixes.

    Arguments: sample_dict (dictionary)
    Returns: isotope (string)
    """
    # create dictionary
    rev_dict = {}

    # iterate through sample_dict to get keys and values and set to new dict
    for key, value in sample_dict.items():
        rev_dict.setdefault(value, set()).add(key)

    # look for keys with same values
    results = [values for key, values in rev_dict.items() if len(values) > 1]

    # check there are two keys with matching values
    if len(results) == 0:
        print("\nError")
        print("Meta file contains an error. Check that there is a shared timepoint between mixes.")
        sys.exit(0)

    # check there is just one result with 2 keys representing two mixes
    if len(results[0]) != 2:
        print("\nError")
        print("Meta file contains an error. Check that there is only one timepoint shared between mixes.")
        sys.exit(0) # exit program

    # check if the isotopes match
    if list(results[0])[0].split("_")[1] == list(results[0])[1].split("_")[1]:
        isotope = list(results[0])[0].split("_")[1]
    else:
        print("\nError")
        print("Isotopes don't match between Mix 1 and Mix 2.")
        sys.exit(0) # exit program

    return isotope

def get_average(df, num):
    """
    This function takes in a dataframe and outputs the averages of the replicates
    in a new dataframe. Note that if the dataframe has NaNs, they will not be
    part of the average.

    Arguments: df (dataframe)
               num (integer), number of unique replicates
    Returns: avg_df (dataframe)
    """
    df[df == 0] = np.nan # set 0's to NaNs

    # create average dataframe with protein id/name info only
    avg_df = df[df.columns[0:4]].copy()

    # intensity values start at 4
    for i in range(4,len(df.columns))[::num]:
        # create new column name
        col_name = df.columns[i][:-2] # example: LM_0
        # create new column of average
        avg_df[col_name] = df.iloc[:,i:i+num].mean(axis=1)

    # replace NaNs with 0's
    avg_df[avg_df.columns[4:]] = avg_df[avg_df.columns[4:]].replace({np.nan:0})

    return avg_df

def get_variance(df, num):
    """
    This function takes in a dataframe and outputs the standard deviations of
    the replicates in a new dataframe. Note that if the dataframe has NaNs, they will not be
    part of the final standard deviation.

    Arguments: df (dataframe)
               num (integer), number of unique replicates
    Returns: sd_df (dataframe)
    """
    df[df == 0] = np.nan # set 0's to NaNs

    # create standard deviation dataframe with protein id/name info only
    sd_df = df[df.columns[0:4]].copy()

    # intensity values start at 4
    for i in range(4,len(df.columns))[::num]:
        # create new column name
        col_name = df.columns[i][:-2] # example: LM_0
        # create new column of standard deviation
        sd_df[col_name] = df.iloc[:,i:i+num].std(axis=1)

    # replace NaNs with 0's
    sd_df[sd_df.columns[4:]] = sd_df[sd_df.columns[4:]].replace({np.nan:0})

    return sd_df

def get_factors(isotope, df):
    """
    This function takes in a dataframe and isotope(timepoint) shared between mixes
    to get the factors needed to normalize each timepoint in the two mixes.
    This returns a new dataframe containing the factors.

    Arguments: df (dataframe)
               isotope (string)
    Returns: factors_df (dataframe)
    """
    # create factors dataframe with protein id/name info only
    factors_df = df[df.columns[0:4]].copy()

    # create dataframe with just the intensity values of common isotope
    isotope_df = df.filter(regex='Intensity\s({}).*'.format(isotope))

    # list of column names
    col_list = isotope_df.columns
    # column name pattern
    regex = r'Intensity\s[^\s]\s([^\s])[^\s]([0-9])([^\s])_([^\s])'
    # group mixes of same stimulant/group/replicate together using key
    key = [0,2,3,1]
    # list of sorted column names
    sorted_col_list = sorted(col_list, key=lambda name: tuple(re.findall(regex,name)[0][i] for i in key))

    # dataframe with column names sorted
    isotope_df = isotope_df[sorted_col_list]

    pre_num = len(isotope_df.columns) # number of columns
    for i in range(0, pre_num)[::2]:
        # create new column name
        col = isotope_df.columns[i].split(" ")
        stim_mix = col[2].split("_")[0][0] # stimulant
        group = col[2].split("_")[0][3] # group
        rep = col[2].split("_")[1] # replicate
        new_col_name = stim_mix + group + "_" + rep + "_avg" # example: LM_A_avg
        # get average of two mixes
        isotope_df[new_col_name] = isotope_df.iloc[:,i:i+2].mean(axis=1)

    post_num = len(isotope_df.columns) # number of columns
    # get the factors by dividing average by each mix's intensity value
    for i in range(pre_num, post_num):
        col_name = isotope_df.columns[i] # column name
        mix1 = "Intensity " + isotope + " " + col_name[0] + "M1" + col_name[1] + "_" + col_name.split("_")[1] # intensity column name mix 1
        mix2 = "Intensity " + isotope + " " + col_name[0] + "M2" + col_name[1] + "_" + col_name.split("_")[1] # intensity column name mix 2
        mix1_col = col_name[0] + "M1" + col_name[1] + "_" + col_name[3] + "_factor" # mix 1 factor column name
        mix2_col = col_name[0] + "M2" + col_name[1] + "_" + col_name[3] + "_factor" # mix 2 factor column name
        isotope_df[mix1_col] = isotope_df.iloc[:,i]/isotope_df[mix1] # mix 1 division
        isotope_df[mix2_col] = isotope_df.iloc[:,i]/isotope_df[mix2] # mix 2 division

    # convert any NaNs to 1 for the factors
    isotope_df[isotope_df.columns[post_num:]] = isotope_df[isotope_df.columns[post_num:]].replace({np.nan:1})

    # convert remaining NaNs to 0
    isotope_df = isotope_df.replace({np.nan:0})
    isotope_df = isotope_df.reindex(sorted(isotope_df.columns), axis=1) # sort and reindex columns

    # add protein info
    factors_df = factors_df.join(isotope_df)

    return factors_df

def get_norm_values(df, f_df):
    """
    This function takes in two dataframes of the pre-normalized values and factors
    to return a new dataframe of normalized values

    Arguments: df (dataframe), intensity values
               f_df (dataframe), factors
    Returns: new_df (dataframe)
    """
    num_cols = len(df.columns) # number of columns
    factors = f_df[f_df.filter(regex='.*(_factor)').columns] # get factors only
    new_df = df.join(factors) # add factors to dataframe

    # normalize replicates by multiplying replicates by its factors
    for i in range(num_cols, len(new_df.columns)):
        col_name = new_df.columns[i] # column name
        repH = "Intensity H " + col_name.split("_")[0] + "_" + col_name.split("_")[1] # H column name
        repM = "Intensity M " + col_name.split("_")[0] + "_" + col_name.split("_")[1] # M column name
        repL = "Intensity L " + col_name.split("_")[0] + "_" + col_name.split("_")[1] # L column name
        new_df[repH] = new_df.iloc[:,i] * new_df[repH] # H normalization
        new_df[repM] = new_df.iloc[:,i] * new_df[repM] # M normalization
        new_df[repL] = new_df.iloc[:,i] * new_df[repL] # L normalization

    # drop factor columns
    new_df.drop(list(new_df.filter(regex = '_factor')), axis = 1, inplace = True)

    # replace NaNs with 0s
    new_df = new_df.replace({np.nan:0})

    return new_df

def get_max(isotope, df, num):
    """
    This function takes the max of the values for the common isotope for each
    stimulant and drops one of the columns.

    Arguments: df (dataframe)
               isotope (string)
               num (integer), number of unique replicates
    Returns: df (dataframe)
    """
    # if dataframe doesn't have "Intensity" in column names
    if not "Intensity" in list(df.columns)[4]:
        # get first sample's name after first four columns containing protein info
        sample = df.columns[4].split("_")[0]

        # get index of first sample with common isotope
        index = df.columns.get_loc("{}_{}".format(sample,isotope))

        # take the max of Mix 1 and Mix 2 of common isotope for each stimulant
        for i in range(index, len(df.columns)-2)[::6]:
            # get column names of mix 1 and mix 2
            mix1 = df.columns[i]
            mix2 = df.columns[i+3]
            # set mix 1 column with max from the two mixes
            df[mix1] = df[[mix1,mix2]].max(axis=1)

        # drop Mix 2 for each common isotope of each stimulant
        df.drop(list(df.filter(regex = '.M2_L')), axis = 1, inplace = True)
    else:
        # get first sample's name after first four columns containing protein info
        sample = df.columns[4].split(" ")[2]

        # get index of first sample with common isotope
        index = df.columns.get_loc("Intensity {} {}".format(isotope,sample))

        # distance from Mix 1 and Mix 2
        dist = num * 3 # number of unique replicates multiplied by 3 isotopes

        # take the max of Mix 1 and Mix 2 of the common isotope for each stimulant
        for i in range(index, len(df.columns)-dist)[::dist]:
            # iterate through each unique replicate
            for j in range(num):
                # get column names of mix 1 and mix 2 for each replicate
                mix1 = df.columns[i+j]
                mix2 = df.columns[i+dist+j]
                df[mix1] = df[[mix1,mix2]].max(axis=1)

        # drop Mix 2 for each common isotope of each stimulant
        df.drop(list(df.filter(regex = 'Intensity\s(L)\s.(M2).*')), axis = 1, inplace = True)

    return df

def rename_columns(sample_dict, df):
    """
    This function renames the columns with the timepoints.

    Arguments: sample_dict (dictionary)
               df (dataframe)
    Returns: df (dataframe)
    """
     # get names of intensity columns
    column_names = list(df.filter(regex="_").columns)

    # create new dictionary
    new_d = {}

    if not "Intensity" in column_names[0]:
        # go through sample dictionary and get new column names with timepoint info
        for key,val in sample_dict.items():
            for col in column_names:
                if key in col:
                    new_d[col] = col[0] + "_" + str(val)
    else:
        # go through sample dictionary and get new column names with timepoint info
        for key,val in sample_dict.items():
            for col in column_names:
                new_col = col.split(" ")[2][1:3] + "_" + col.split(" ")[1] # reformat column name
                if key in new_col:
                    new_d[col] = col.split(" ")[2][0] + col.split(" ")[2][3] + "_" + str(val) + "_" +  col.split(" ")[2][5]

    # rename columns
    df = df.rename(columns=new_d)
    # get new column names
    new_col_list = list(df.filter(regex="_").columns)
    # sort column names
    sorted_list = natsorted(new_col_list)
    # save dataframe with sorted columns
    df = df[list(df.columns)[:4] + sorted_list]

    return df

def random_impute(arr):
    """
    This function takes in an array and imputes a random number for any 0's.

    Argument: arr (array)
    Returns: arr (array)
    """
    # iterate through array to assign a random number from 500 to 1000 if 0
    for i in range(len(arr)):
        if arr[i] == 0:
            arr[i] = np.random.randint(low=500,high=1000)

    return arr

def colToExcel(col):
    """
    This function takes in an integer of the column index and returns the excel column letter.

    Argument: col (integer)
    Returns: excelCol (string)
    """
    # col is 1 based
    excelCol = str() # returns string
    div = col
    while div:
        (div, mod) = divmod(div-1, 26) # returns two numbers ~ div is quotient, mod is remainder
        excelCol = chr(mod + 65) + excelCol # returns a character whose unicode point is an integer (65 is A)

    return excelCol

def main(args):

    # get path to control and experiment proteinGroups.txt files
    con_file = args.con
    exp_file = args.exp

    # remove unwanted columns
    con_df = filter_maxquant_col(con_file)
    exp_df = filter_maxquant_col(exp_file)
    print("Number of protein hits for experiment group is {}.".format(exp_df.shape[0]))
    print("Number of protein hits for control group is {}.".format(con_df.shape[0]))

    # check that the right file was selected for control group
    if con_df.columns[14].split(" ")[2].split("_")[0][3] != "C":
        print("\nError")
        print("The wrong proteinGroups.txt file was selected for the control group.")
        print("Please check the path to the correct file.")
        sys.exit(0) # exit program

    # check that the right file was selected for experimental group
    if exp_df.columns[14].split(" ")[2].split("_")[0][3] == "C":
        print("\nError")
        print("The wrong proteinGroups.txt file was selected for the experimental group.")
        print("Please check the path to the correct file.")
        sys.exit(0) # exit program

    # count number of replicates per condition
    con_rep_counter = pd.Series(x[12] for x in con_df.columns[14:]).map(Counter).sum()
    exp_rep_counter = pd.Series(x[12] for x in exp_df.columns[14:]).map(Counter).sum()

    # get number of conditions/stimulants (i.e. L,O,P,R) per group
    con_stim = len(con_rep_counter)
    exp_stim = len(exp_rep_counter)

    # check if the count is the same for each group
    if con_stim != exp_stim:
        print("\nError")
        print("The number of conditions don't match between experimental and control groups.")
        sys.exit(0) # exit

    # check if the conditions are the same for experimental and control groups
    if list(con_rep_counter) != list(exp_rep_counter):
        print("\nError")
        print("The conditions are not the same between experimental and control groups.")
        sys.exit(0) # exit

    # get number of timepoints per condition by dividing by sum of replicates
    # by number of conditions and by number of isotopes for triple SILAC
    exp_count = (sum(exp_rep_counter.values())/len(exp_rep_counter))/3
    con_count = (sum(con_rep_counter.values())/len(con_rep_counter))/3

    # check if the count is the same for each group
    if exp_count != con_count:
        print("\nError")
        print("The number of timepoints don't match between experimental and control groups.")
        sys.exit(0) # exit

    # boolean if two mixes are present in proteinGroups.txt file
    is_two_mixes = True
    # number of replicates per condition/timepoint
    rep_num = int((sum(exp_rep_counter.values())/len(exp_rep_counter))/6)

    # if only 3 timepoints, then one mix is present
    if exp_count==3 and con_count==3:
        is_two_mixes = False
        rep_num = int((sum(exp_rep_counter.values())/len(exp_rep_counter))/3)

    # remove unwanted rows
    con_df = filter_maxquant_row(con_df)
    exp_df = filter_maxquant_row(exp_df)
    print("Number of protein hits for experiment group is {} after removing unwanted rows.".format(exp_df.shape[0]))
    print("Number of protein hits for control group is {} after removing unwanted rows.".format(con_df.shape[0]))

    # merge dataframes - left join on experiment dataframe
    merged_df = pd.merge(exp_df, con_df, how="left", on="Majority protein IDs")

    # delete unused dataframes
    del exp_df
    del con_df

    # remove duplicate and unneeded columns
    remove_col_list = ["Protein names_y", "Gene names_y", "Fasta headers_y",
                       "Peptides_y", "Unique peptides_y", "Sequence coverage [%]_y",
                       "Mol. weight [kDa]_y","Sequence length_y", "Q-value_y",
                       "Score_y"]
    merged_df = merged_df.drop(remove_col_list, axis=1)

    # remove _x from columns that came from merging dataframes
    merged_df = merged_df.rename(columns = lambda x : str(x).rstrip("_x"))

    # create new dataframe with replicates grouped together
    sorted_merged_df = sort_raw_dataframe(merged_df, rep_num)

    # generate count table
    merged_count_df = create_count_df(sorted_merged_df, rep_num)

    # save raw intensity values
    if args.prefix:
        sorted_merged_df.to_excel(args.prefix + "_raw_filtered.xlsx", index=False, na_rep="NaN")
    else:
        sorted_merged_df.to_excel("raw_filtered.xlsx", index=False, na_rep="NaN")

    # save counts to file
    if args.prefix:
        workbook = pd.ExcelWriter(args.prefix + "_counts.xlsx", engine="xlsxwriter")
    else:
        workbook = pd.ExcelWriter("counts.xlsx", engine="xlsxwriter")

    # write count table to worksheet
    merged_count_df.to_excel(workbook, sheet_name="Counts", index=False)

    # number of columns in count table to use for getting excel column letter
    numCols = int((len(merged_count_df.columns)-5)/3)

    # write info describing count table
    line1 = "Columns {} to {} give the number of replicates with nonzero intensity values.".format(colToExcel(5), colToExcel(numCols*2+4))
    line2 = "Columns {} to {} assigns 0 or 1 if it meets the following criteria.".format(colToExcel(numCols*2+5),colToExcel(numCols*3+4))
    line3 = "1 is given if the number of biological replicates >= 2/3 * total unique replicates and the number of control replicates <= 1/3 * total unique replicates, else 0."
    line4 = "For this dataset, the number of biological replicates must be >= {} and the number of control replicates must be <= {} to get 1.".format(rep_num*2/3, rep_num*1/3)
    line5 = "Column {} gives the sum of conditions/timepoints that meet the criteria.".format(colToExcel(numCols*3+5))

    # save as dataframe
    info_df = pd.DataFrame({"info":[line1,line2,line3,line4,line5]})

    # add new worksheet describing counts table
    info_df.to_excel(workbook, sheet_name="Info", index=False)

    # close and save counts excel file
    workbook.save()

    # filter out rows with count_sum == 0, no stimulant/timepoint met criteria
    filtered_df = sorted_merged_df[merged_count_df["sum"] != 0].reset_index()
    print("Number of protein hits after filtering out those that didn't meet count criteria is {}.".format(filtered_df.shape[0]))

    # create dictionary from meta file containing timepoint info
    samples = create_sample_dict(args.meta)

    if is_two_mixes:
        # get isotope of timepoint shared between mixes
        common_isotope = get_common_isotope(samples)

    # get counts for experiment group only and which meets criteria
    exp_count_df = merged_count_df[merged_count_df["sum"] != 0].reset_index().drop("index", axis = 1)
    exp_count_df = exp_count_df[exp_count_df.columns[0:4]].join(exp_count_df[exp_count_df.filter(regex='.*[^C\d]_.*').columns]) # remove control group columns
    exp_count_df.columns = exp_count_df.columns.str.replace("M_", "_") # remove extra letter
    if is_two_mixes:
        exp_count_df = get_max(common_isotope, exp_count_df, rep_num)
    exp_count_df = rename_columns(samples, exp_count_df)

    # delete unused dataframes
    del merged_count_df
    del merged_df

    # create dataframe of control and experimental group intensities
    both_df = filtered_df[filtered_df.columns[:5]].join(filtered_df[filtered_df.columns[12:]]).drop("index", axis=1)
    both_df[both_df == 0] = np.nan # set 0's to NaNs

    # create dataframe for normalizing replicates between mixes
    if is_two_mixes: # two mixes only
        both_factors_df = get_factors(common_isotope, both_df)
        norm_both_df = get_norm_values(both_df, both_factors_df)
        norm_both_df = get_max(common_isotope, norm_both_df, rep_num)
        norm_both_df = rename_columns(samples, norm_both_df)
    else: # one mix only
        norm_both_df = rename_columns(samples, both_df)

    # get averages of replicates
    norm_avg_df = get_average(norm_both_df,rep_num)

    # save as excel files
    if args.prefix:
        norm_both_df.to_excel(args.prefix + "_all_norm.xlsx", index=False, na_rep=0)
        norm_avg_df.to_excel(args.prefix + "_norm_avg.xlsx", index=False)
    else:
        norm_both_df.to_excel("all_norm.xlsx", index=False)
        norm_avg_df.to_excel("norm_avg.xlsx", index=False)

    # get averages from only experimental group
    exp_norm_avg_df = norm_avg_df[norm_avg_df.columns[0:4]].join(norm_avg_df[norm_avg_df.filter(regex='.*[^C]_.*').columns]) # remove control group columns
    exp_norm_avg_df.columns = exp_norm_avg_df.columns.str.replace("M_", "_") # remove extra letter

    # get standard deviations of normalized replicates
    norm_std_df = get_variance(norm_both_df, rep_num)

    # get standard deviations of just the experimental group
    exp_norm_std_df = norm_std_df[norm_std_df.columns[0:4]].join(norm_std_df[norm_std_df.filter(regex='.*[^C]_.*').columns]) # remove control group columns
    exp_norm_std_df.columns = exp_norm_std_df.columns.str.replace("M_", "_") # remove extra letter

    # save experimental group stats to file
    if args.prefix:
        writer = pd.ExcelWriter(args.prefix + '_exp_stats.xlsx', engine='xlsxwriter')
    else:
        writer = pd.ExcelWriter('_exp_stats.xlsx', engine='xlsxwriter')

    # write each dataframe to a different worksheet
    exp_norm_avg_df.to_excel(writer, sheet_name="Mean", index=False)
    exp_norm_std_df.to_excel(writer, sheet_name="Std Dev", index=False)
    exp_count_df.to_excel(writer, sheet_name="Count", index=False)

    # close the Pandas Excel writer and output the Excel file.
    writer.save()

    np.random.seed(0) # to get same random numbers every time when running the program

    norm_both_df = norm_both_df.replace({np.nan:0}) # replace NaNs with 0's

    # iterate through each column and row to replace 0's with random number
    for i in range(4,len(norm_both_df.columns[4:])):
        vals = norm_both_df[norm_both_df.columns[i]].values
        new_vals = random_impute(vals)
        norm_both_df[norm_both_df.columns[i]] = new_vals

    # save as excel file
    if args.prefix:
        norm_both_df.to_excel(args.prefix + "_all_norm_rand.xlsx", index=False)
    else:
        norm_both_df.to_excel("all_norm_rand.xlsx", index=False)

    print("\nProgram ran successfully.")


if __name__ == "__main__":

    # create arguments
    p = argparse.ArgumentParser(description=__doc__, prog="process_maxquant.py",
        usage="%(prog)s -exp <path/to/file> -con <path/to/file> --meta <path/to/file> [options]", add_help=True)
    p.add_argument("-exp", help="path to experiment proteinGroups.txt file", required=True)
    p.add_argument("-con", help="path to control proteinGroups.txt file", required=True)
    p.add_argument("--prefix", help="prefix to use for naming output files")
    p.add_argument("--meta", help="tab-delimited file listing timepoints to isotope/mix", required=True)

    # parse arguments
    args = p.parse_args()
    # run program with arguments
    main(args)
