# process_maxquant
Python script for processing MaxQuant files for downstream analysis.

The script is specific for triple SILAC quantitative proteomic analysis involving two mixes with overlapping timepoints as shown in the diagram below. 


<img src="https://github.com/jonessarae/process_maxquant/blob/media/triple_silac1.PNG">

## Table of contents

   * [Installation](#installation)   
   * [Input](#Input)
   * [Naming scheme for MaxQuant](#Naming-scheme-for-MaxQuant)
   * [How to prepare meta file](#how-to-prepare-meta-file)
   * [Usage](#usage)
   * [Output](#output)
   * [Prism template](#Prism-template)
   * [Other resources](#Other-resources)
   

## Installation

You will need Python 3 already on your computer.

To get Python 3, you can download Anaconda at https://www.anaconda.com/distribution/.

Once installed, open the Anaconda Prompt and install *git* with the following command:

<pre>
conda install -c anaconda git
</pre>

Next, install the following Python packages:

* numpy
* pandas
* natsort

These packages can be installed with the following command:

<pre>
pip install numpy pandas natsort
</pre>

Download the python script:

<pre>
git clone https://github.com/jonessarae/process_maxquant.git
</pre>

The script, *process_maxquant.py*, is available in the folder *process_maxquant*.

For reference, the versions used to create the python script are listed in the file *versions.txt*.

## Input 

This script requires the MaxQuant output file, *proteinGroups.txt*, from the control and experimental groups. They must contain data collected from three isotopes and from at least two replicates. The *proteinGroups.txt* can contain intensity data from either one mix or two mixes. When two mixes are included, this script performs additional steps for normalizing the intensity data between the two mixes. 

The script also requires an additional file called the *meta file* that contains timepoint information and specifies the overlapping timepoint between two mixes for normalization.

## Naming scheme for MaxQuant files

This script depends on the following naming scheme for the intensity columns in the proteinGroups.txt containing either one mix or two mixes:

<img src="https://github.com/jonessarae/process_maxquant/blob/media/naming.PNG">

This information is used to sort the columns. The names must also be unique.

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


## Output

* *raw_filtered.xlsx*: filtered, merged excel file with raw intensity values from control and experimental groups
* *counts.xlsx*: excel file of replicate counts
* *all_norm.xlsx.xlsx*: excel file with normalized values for each replicate for both control and experimental groups
* *all_norm_rand.xlsx*: same as above but with random imputation (range 500-1000) applied to 0's
* *norm_avg.xlsx.xlsx*: excel file with averages for both control and experimental groups
* *exp_stats.xlsx*: excel file with mean, standard deviation, and counts for each timepoint in experimental group

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

### Other resources

To combine *proteinGroups.txt* files of Mix 1 and Mix 2 from the same group, go to https://github.com/jonessarae/merge_maxquant for the python script. 

For tutorials to run proteomic software for downstream analyses, go to https://github.com/jonessarae/proteomic_tutorials.

