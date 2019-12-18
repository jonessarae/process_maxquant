# process_maxquant
Python script for processing MaxQuant files for downstream analysis.

The script is specific for triple SILAC quantitative proteomic analysis involving two mixes as shown in the diagram below.


<img src="https://github.com/jonessarae/process_maxquant/blob/master/triple_silac.PNG">

## Software

## How to run the program

<pre>
python process_maxquant.py -exp path/to/file -con path/to/file --meta path/to/file [options]
</pre> 

Example:
<pre>
python process_maxquant.py -con Mix12_Con/proteinGroups.txt -exp Mix12_Myd/proteinGroups.txt --prefix Mix12_ConMyd --meta info.txt
</pre> 

Parameters:
* __-exp__: path to experiment proteinGroups.txt file, *required*
* __-con__: path to control proteinGroups.txt file, *required*
* __--prefix__: prefix to use for naming output files
* __--meta__: tab-delimited file listing timepoints to isotope/mix, *required*

## How to prepare meta file

The meta file must be in the following format:

<pre>
M1_M    0	
M1_L    10
M2_L    10
M2_M    30
M1_H    60
M2_H    120
</pre>

The first column lists the isotope/mix and the second column the associated timepoints. 

The columns are separated by a TAB. 

There must be one timepoint that is shared between the two mixes. In this case, it is M1_L and M2_L.

## Files generated

* *raw_filtered.xlsx*: filtered, merged excel file with raw intensity values from control and experimental groups
* *counts.xlsx*: excel file of replicate counts
* *exp_stats.xlsx*: excel file with mean, standard deviation, and counts for each timepoint in experimental group
* *all_norm.xlsx.xlsx*: excel file with normalized values for each replicate for both control and experimental groups
* *all_norm_rand.xlsx*.xlsx: same as above but with random imputation (range 500-1000) applied to 0's
* *norm_avg.xlsx.xlsx*: excel file with normalized averages for both control and experimental groups

## Workflow

This workflow diagram shows the major steps in the python script. 

Shown in orange are all the files generated from the script, and in purple are the downstream software programs that will utilize those files.

<img src="https://github.com/jonessarae/process_maxquant/blob/master/diagram.png">

## Template for Prism

Included is a text file called *prism_output_template.txt* that contains excel formulas specific for the output file, *exp_stats.xlsx*. 

Everything in the text file can be copied to a new sheet in *exp_stats.xlsx*. 

It will then populate with values based on the protein ID (cell B1). 

The table generated can then be used as input for Prism. 

Example:

<img src="https://github.com/jonessarae/process_maxquant/blob/master/prism_table_example.PNG">

