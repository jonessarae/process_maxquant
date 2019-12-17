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

# shut off settingwithcopy warning
pd.options.mode.chained_assignment = None  # default='warn'

"""
Purpose: Pipeline for processing files generated from MaxQuant.

This script assumes that the data contains three replicates for each condition in the experimental and control
groups. It also assumes that the replicates are the same for both experimental and control
groups.

To use:
python process_maxquant.py -exp <path/to/file> -con <path/to/file> --prefix <string> --meta <path/to/file>
Example:
python process_maxquant.py -con Mix12_Con_txt/proteinGroups.txt -exp Mix12_Myd_txt/proteinGroups.txt --prefix Mix12_ConMyd --meta info.txt

Files that are generated:
<prefix>_raw_filtered.xlsx: filtered, merged excel file with raw intensity values
<prefix>_counts.xlsx: excel file of replicate counts
<prefix>_exp_stats.xlsx: excel file with averages, std dev, and counts for experimental group
<prefix>_all_norm.xlsx.xlsx: excel file with normalized values for each replicate for both control and experimental groups
<prefix>_all_norm_rand.xlsx.xlsx: same as above but with random imputation for 0's

"""
__author__ = "Sara Jones"
__email__ = "jonessarae@gmail.com"
__doc__ = "Pipeline that processes files from MaxQuant."

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

    # create dataframe with non-intensity columns
    no_intensity_df = df[col_list]

    # create dataframe with intensity columns
    intensity_df = df[df.filter(regex="Intensity H |Intensity M |Intensity L ").columns]

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

def sort_raw_dataframe(df):
    """
    This function takes in a dataframe and sorts the dataframe so that replicates
    for each stimulant/timepoint are grouped together.

    Arguments: df (dataframe)
    Returns: df (dataframe)
    """
    # get list of intensity column names
    intensity_col_list = list(df.columns[11:])
    # regular expression pattern for column name containing Intensity
    regex = r'Intensity\s([^\s])\s([^\s])[^\s]([0-9])([^\s])_([^\s])'
    # key for sorting columns containing Intensity by the following order:
    # stimulant, group (control vs experiment), mix, isotope, replicate
    key = [1,3,2,0,4]
    # sort list of intensity columns
    sorted_list = sorted(intensity_col_list, key=lambda name: tuple(re.findall(regex,name)[0][i] for i in key))

    # bug checker - will exit program if the following condition is not met
    assert len(sorted_list)%3==0,"Missing replicate for one of the conditions. Need 3 replicates per condition."

    # list of columns of protein info
    protein_list = list(df.columns[:11])

    # replace old columns with new columns
    df = df[protein_list+sorted_list]

    return df

def create_count_df(df):
    """
    This function takes in a dataframe and creates a new dataframe counting the
    number of replicates with an intensity value greater than 0.

    If the count is <= 1 in the control group and >= 2 in the experiment group,
    then the corresponding stimulant/timepoint is set to 1, else it is set to 0.

    A column summing up every instance that a stimulant/timepoint is 1 is provided
    for filtering.

    Arguments: df (dataframe)
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
    for i in range(len(intensity_df.columns))[::3]:
        # create new column name
        col_name_a = intensity_df.columns[i].split(" ")
        isotope_a = col_name_a[1]
        stim_mix_a = col_name_a[2].split("_")[0]
        new_col_name_a = stim_mix_a + "_" + isotope_a # example: LMIC_H
        # create new column with count info
        count_df[new_col_name_a]=(intensity_df.iloc[:,i:i+3] != 0).astype(int).sum(axis=1)

    # get list of count column names
    count_col_list = list(count_df.columns[11:])
    # regular expression pattern for count column names
    regex = r'([^\s]*)([^\s])_([^\s])'
    # key for sorting count columns by stimulant/mix, isotope, control/experiment
    key = [0,2,1]
    # sort list of count columns
    new_list = sorted(count_col_list, key=lambda name: tuple(re.findall(regex,name)[0][i] for i in key))

    # bug checker
    assert len(new_list)%2==0,"Missing control and/or experimental group in data."

    # set count dataframe with sorted count columns and first 4 columns of protein info
    count_df = count_df[list(count_df.columns[:4]) + new_list]

    # get number of columns in count dataframe
    num_cols = len(count_df.columns)

    # loop through each stimulant/mix and check if the number of replicates
    # for controls is <= 1 and for experiments is >= 2
    for i in range(4, num_cols)[::2]:
        # create new column name
        col_name_b = count_df.columns[i].split("C_")
        stim_mix_b = col_name_b[0]
        isotope_b = col_name_b[1]
        new_col_name_b = stim_mix_b + "_" + isotope_b # example: PM1_M
        # create new column with value 1 if it meets condition
        count_df.loc[(count_df[count_df.columns[i]] <= 1) & (count_df[count_df.columns[i+1]] >= 2), new_col_name_b] = 1

    # loop through each column to replace any NaNs with 0
    for i in range(num_cols, len(count_df.columns)):
        count_df[count_df.columns[i]].fillna(0, inplace=True)
    # set numbers as int type in condition columns
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
        print("Meta file contains an error. Check that there is a shared timepoint between mixes.")
        sys.exit(0)

    # check there is just one result with 2 keys representing two mixes
    if len(results[0]) !=2:
        print("Meta file contains an error. Check that there is only one timepoint shared between mixes.")
        sys.exit(0) # exit program

    # check if the isotopes match
    if list(results[0])[0].split("_")[1] == list(results[0])[1].split("_")[1]:
        isotope = list(results[0])[0].split("_")[1]
    else:
        print("Isotopes don't match between Mix 1 and Mix 2.")
        sys.exit(0) # exit program

    return isotope

def get_average(df):
    """
    This function takes in a dataframe and outputs the averages of three replicates
    in a new dataframe. Note that if the dataframe has NaNs, they will not be
    part of the average.

    Arguments: df (dataframe)
    Returns: avg_df (dataframe)
    """
    df[df == 0] = np.nan # set 0's to NaNs

    # create average dataframe with protein id/name info only
    avg_df = df[df.columns[0:4]].copy()

    # intensity values start at 4
    for i in range(4,len(df.columns))[::3]:
        # create new column name
        col_name = df.columns[i][:-2] # example: LM_0
        # create new column of average
        avg_df[col_name] = df.iloc[:,i:i+3].mean(axis=1)

    # replace NaNs with 0's
    avg_df[avg_df.columns[4:]] = avg_df[avg_df.columns[4:]].replace({np.nan:0})

    return avg_df

def get_variance(df):
    """
    This function takes in a dataframe and outputs the standard deviations of
    three replicates in a new dataframe. Note that if the dataframe has NaNs, they will not be
    part of the final standard deviation.

    Arguments: df (dataframe)
    Returns: sd_df (dataframe)
    """
    df[df == 0] = np.nan # set 0's to NaNs

    # create standard deviation dataframe with protein id/name info only
    sd_df = df[df.columns[0:4]].copy()

    # intensity values start at 4
    for i in range(4,len(df.columns))[::3]:
        # create new column name
        col_name = df.columns[i][:-2] # example: LM_0
        # create new column of standard deviation
        sd_df[col_name] = df.iloc[:,i:i+3].std(axis=1)

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
    regex = r'Intensity\s[^\s]\s([^\s])[^\s]([0-9])([^\s])_([^\s])'
    # group mixes of same stimulant/group/replicate together
    key = [0,2,3,1]
    # list of sorted column names
    sorted_col_list = sorted(col_list, key=lambda name: tuple(re.findall(regex,name)[0][i] for i in key))
    # dataframe with column names sorted
    isotope_df = isotope_df[sorted_col_list]
    pre_num = len(isotope_df.columns) # number of columns
    for i in range(0, pre_num)[::2]:
        # create new column name
        col = isotope_df.columns[i].split(" ")
        stim_mix = col[2].split("_")[0][0]
        group = col[2].split("_")[0][3]
        rep = col[2].split("_")[1]
        new_col_name = stim_mix + group + "_" + rep + "_avg" # example: LM_A_avg
        # get average of two mixes
        isotope_df[new_col_name] = isotope_df.iloc[:,i:i+2].mean(axis=1)
    post_num = len(isotope_df.columns) # number of column
    # get the factors by dividing combined average by each mix's average
    for i in range(pre_num, post_num):
        col_name = isotope_df.columns[i]
        mix1 = "Intensity " + isotope + " " + col_name[0] + "M1" + col_name[1] + "_" + col_name.split("_")[1]
        mix2 = "Intensity " + isotope + " " + col_name[0] + "M2" + col_name[1] + "_" + col_name.split("_")[1]
        mix1_col = col_name[0] + "M1" + col_name[1] + "_" + col_name[3] + "_factor"
        mix2_col = col_name[0] + "M2" + col_name[1] + "_" + col_name[3] + "_factor"
        isotope_df[mix1_col] = isotope_df.iloc[:,i]/isotope_df[mix1]
        isotope_df[mix2_col] = isotope_df.iloc[:,i]/isotope_df[mix2]

    # convert any NaNs to 1 for the factors
    isotope_df[isotope_df.columns[post_num:]] = isotope_df[isotope_df.columns[post_num:]].replace({np.nan:1})
    # convert remaining NaNs to 0
    isotope_df = isotope_df.replace({np.nan:0})
    isotope_df = isotope_df.reindex(sorted(isotope_df.columns), axis=1) # sort columns

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
        col_name = new_df.columns[i]
        repH = "Intensity H " + col_name.split("_")[0] + "_" + col_name.split("_")[1]
        repM = "Intensity M " + col_name.split("_")[0] + "_" + col_name.split("_")[1]
        repL = "Intensity L " + col_name.split("_")[0] + "_" + col_name.split("_")[1]
        new_df[repH] = new_df.iloc[:,i] * new_df[repH]
        new_df[repM] = new_df.iloc[:,i] * new_df[repM]
        new_df[repL] = new_df.iloc[:,i] * new_df[repL]

    # drop factors
    new_df.drop(list(new_df.filter(regex = '_factor')), axis = 1, inplace = True)

    # replace NaNs with 0s
    new_df = new_df.replace({np.nan:0})

    return new_df

def get_max(isotope, df):
    """
    This function takes the max of the values for the common isotope for each
    stimulant and drops one of the columns.

    Arguments: df (dataframe)
               isotope (string)
    Returns: df (dataframe)
    """
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

        # take the max of Mix 1 and Mix 2 of common isotope for each stimulant
        for i in range(index, len(df.columns)-7)[::9]:
            # get column names of mix 1 and mix 2 for each replicate
            # first replicate
            mix1_a = df.columns[i]
            mix2_a = df.columns[i+9]
            df[mix1_a] = df[[mix1_a,mix2_a]].max(axis=1)
            # second replicate
            mix1_b = df.columns[i+1]
            mix2_b = df.columns[i+10]
            df[mix1_b] = df[[mix1_b,mix2_b]].max(axis=1)
            # third replicate
            mix1_c = df.columns[i+2]
            mix2_c = df.columns[i+11]
            df[mix1_c] = df[[mix1_c,mix2_c]].max(axis=1)

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
                new_col = col.split(" ")[2][1:3] + "_" + col.split(" ")[1]
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

def main(args):

    # get path to control and experiment proteinGroups.txt files
    con_file = args.con
    exp_file = args.exp

    # check for proteinGroups.txt files
    if not exp_file.endswith("proteinGroups.txt"):
        print("Missing proteinGroups.txt file for experiment group.")
        print("Please check the path to the proteinGroups.txt file.")
        sys.exit(0) # exit program
    if not con_file.endswith("proteinGroups.txt"):
       print("Missing proteinGroups.txt file for control group.")
       print("Please check the path to the proteinGroups.txt file.")
       sys.exit(0) # exit program

    # remove unwanted columns
    con_df = filter_maxquant_col(con_file)
    exp_df = filter_maxquant_col(exp_file)
    print("Number of protein hits for experiment group is {}.".format(exp_df.shape[0]))
    print("Number of protein hits for control group is {}.".format(con_df.shape[0]))

    # check that the right files were selected
    if con_df.columns[14].split(" ")[2].split("_")[0][3] != "C":
        print("The wrong proteinGroups.txt file was selected for the control group.")
        print("Please check the path to the correct file.")
        sys.exit(0) # exit program

    if exp_df.columns[14].split(" ")[2].split("_")[0][3] == "C":
        print("The wrong proteinGroups.txt file was selected for the experimental group.")
        print("Please check the path to the correct file.")
        sys.exit(0) # exit program

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
    sorted_merged_df = sort_raw_dataframe(merged_df)

    # generate count table
    merged_count_df = create_count_df(sorted_merged_df)

    # save raw intensity values and counts to separate excel files
    if args.prefix:
        sorted_merged_df.to_excel(args.prefix + "_raw_filtered.xlsx", index=False, na_rep="NaN")
        merged_count_df.to_excel(args.prefix + "_counts.xlsx", index=False)
    else:
        sorted_merged_df.to_excel("raw_filtered.xlsx", index=False, na_rep="NaN")
        merged_count_df.to_excel("counts.xlsx", index=False)

    # filter out rows with count_sum == 0, no stimulant/timepoint met criteria
    filtered_df = sorted_merged_df[merged_count_df["sum"] != 0].reset_index()
    print("Number of protein hits after filtering out those that didn't meet count criteria is {}.".format(filtered_df.shape[0]))

    # create dictionary from meta file containing timepoint info
    samples = create_sample_dict(args.meta)

    # get isotope of timepoint shared between mixes
    common_isotope = get_common_isotope(samples)

    # get counts for experiment group only and which meets criteria
    exp_count_df = merged_count_df[merged_count_df["sum"] != 0].reset_index().drop("index", axis = 1)
    exp_count_df = exp_count_df[exp_count_df.columns[0:4]].join(exp_count_df[exp_count_df.filter(regex='.*[^C\d]_.*').columns])
    exp_count_df.columns = exp_count_df.columns.str.replace("M_", "_") # remove extra letter
    exp_count_df = get_max(common_isotope, exp_count_df)
    exp_count_df = rename_columns(samples, exp_count_df)

    # delete unused dataframes
    del merged_count_df
    del merged_df

    # create dataframe of control and experimental group intensities
    both_df = filtered_df[filtered_df.columns[:5]].join(filtered_df[filtered_df.columns[12:]]).drop("index", axis=1)
    both_df[both_df == 0] = np.nan # set 0's to NaNs

    # create dataframe for normalizing replicates
    both_factors_df = get_factors(common_isotope, both_df)
    norm_both_df = get_norm_values(both_df, both_factors_df)
    norm_both_df = get_max(common_isotope, norm_both_df)
    norm_both_df = rename_columns(samples, norm_both_df)

    # get averages of normalized replicates
    norm_avg_df = get_average(norm_both_df)

    # save as excel files
    if args.prefix:
        norm_both_df.to_excel(args.prefix + "_all_norm.xlsx", index=False, na_rep=0)
        norm_avg_df.to_excel(args.prefix + "_norm_avg.xlsx", index=False)
    else:
        norm_both_df.to_excel("all_norm.xlsx", index=False)
        norm_avg_df.to_excel("norm_avg.xlsx", index=False)

    # get averages from only experimental group
    exp_norm_avg_df = norm_avg_df[norm_avg_df.columns[0:4]].join(norm_avg_df[norm_avg_df.filter(regex='.*[^C]_.*').columns])
    exp_norm_avg_df.columns = exp_norm_avg_df.columns.str.replace("M_", "_") # remove extra letter

    # get standard deviations of normalized replicates
    norm_std_df = get_variance(norm_both_df)

    # get standard deviations of just the experimental group
    exp_norm_std_df = norm_std_df[norm_std_df.columns[0:4]].join(norm_std_df[norm_std_df.filter(regex='.*[^C]_.*').columns])
    exp_norm_std_df.columns = exp_norm_std_df.columns.str.replace("M_", "_") # remove extra letter

    # save experimental group stats to file
    if args.prefix:
        writer = pd.ExcelWriter(args.prefix + '_exp_stats.xlsx', engine='xlsxwriter')
    else:
        writer = pd.ExcelWriter('_exp_stats.xlsx', engine='xlsxwriter')

    # write each dataframe to a different worksheet.
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


if __name__ == "__main__":

    # create arguments
    p = argparse.ArgumentParser(description=__doc__, prog="process_maxquant.py",
        usage="%(prog)s -exp <path/to/file> -con <path/to/file> [options]", add_help=True)
    p.add_argument("-exp", help="path to experiment proteinGroups.txt file", required=True)
    p.add_argument("-con", help="path to control proteinGroups.txt file", required=True)
    p.add_argument("--prefix", help="prefix to use for naming output files")
    p.add_argument("--meta", help="tab-delimited file listing timepoints to isotope/mix", required=True)

    # parse arguments
    args = p.parse_args()
    # run program with arguments
    main(args)
