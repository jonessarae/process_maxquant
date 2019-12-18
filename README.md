# process_maxquant
Python script for processing MaxQuant files for downstream analysis.

<img src="https://github.com/jonessarae/process_maxquant/blob/master/triple_silac.PNG">

## How to use

<pre>
python process_maxquant.py -exp path/to/file -con path/to/file --meta path/to/file [options]
</pre> 

Example:
<pre>
python process_maxquant.py -con Mix12_Con/proteinGroups.txt -exp Mix12_Myd/proteinGroups.txt --prefix Mix12_ConMyd --meta info.txt
</pre> 

Parameters:
* __-exp__: path to experiment proteinGroups.txt file, *required
* __-con__: path to control proteinGroups.txt file, *required
* __--prefix__: prefix to use for naming output files
* __--meta__: tab-delimited file listing timepoints to isotope/mix, *required


## Diagram
<img src="https://github.com/jonessarae/process_maxquant/blob/master/diagram.png">
