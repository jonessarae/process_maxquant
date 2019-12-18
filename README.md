# process_maxquant
Python script for processing MaxQuant files for downstream analysis.

The script is specific for triple SILAC quantitative proteomic analysis involving two mixes as shown in the diagram below.


<img src="https://github.com/jonessarae/process_maxquant/blob/media/triple_silac.PNG">

## Table of contents

   * [Installation](#installation)
   * [Usage](#usage)
      * [How to prepare meta file](#how-to-prepare-meta-file)
   * [Output](#output)
   * [Prism template](#Prism-template)
   

## Installation

You will need Python 3 already on your computer and the following packages installed:

* numpy
* pandas
* natsort

These packages can be installed with the following command:

<pre>
pip install numpy pandas natsort
</pre>

For reference, the versions used to create the python script are listed in the file *versions.txt*.

Next, download the python script:

<pre>
git clone https://github.com/jonessarae/process_maxquant.git
</pre>

The script, *process_maxquant.py* is available in the folder *process_maxquant*.

## Usage

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

There must be one timepoint that is shared between the two mixes, M1 and M2. In this case, it is M1_L and M2_L.

An example, *info.txt*, is provided. 

## Output

* *raw_filtered.xlsx*: filtered, merged excel file with raw intensity values from control and experimental groups
* *counts.xlsx*: excel file of replicate counts
* *exp_stats.xlsx*: excel file with mean, standard deviation, and counts for each timepoint in experimental group
* *all_norm.xlsx.xlsx*: excel file with normalized values for each replicate for both control and experimental groups
* *all_norm_rand.xlsx*: same as above but with random imputation (range 500-1000) applied to 0's
* *norm_avg.xlsx.xlsx*: excel file with normalized averages for both control and experimental groups

## Workflow

This workflow diagram shows the major steps in the python script. 

Shown in orange are all the files generated from the script, and in purple are the downstream software programs that will utilize those files.

<img src="https://github.com/jonessarae/process_maxquant/blob/media/diagram.png">

## Prism template

Included is a text file called *prism_output_template.txt* that contains excel formulas specific for the output file, *exp_stats.xlsx*. 

Everything in the text file can be copied to a new sheet in *exp_stats.xlsx*. 

It will then populate with values based on the protein ID (cell B1). 

The table generated can then be used as input for Prism. 

Example:

<img src="https://github.com/jonessarae/process_maxquant/blob/media/prism_table_example.PNG">

