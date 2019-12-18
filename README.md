# process_maxquant
Python script for processing MaxQuant files for downstream analysis.

The script is specific for triple SILAC quantitative proteomic analysis involving two mixes as shown in the diagram below.


<img src="https://github.com/jonessarae/process_maxquant/blob/master/triple_silac.PNG">

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

## Other files



## Workflow
<img src="https://github.com/jonessarae/process_maxquant/blob/master/diagram.png">
