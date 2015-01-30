# README #

**mixTAPE: mix of Tools for Analysis of Phylogenetics and Evolution**

http://www.bitbucket.org/jbpease/mixtape

James B. Pease
jbpease@indiana.edu
https://mypage.iu.edu/~jbpease

This is a set of simple scripts to perform general tasks in data manipulation, bioinformatics, evolutionary analysis and phylogenetics.
Development of these tools is ongoing and there may be daily builds of some tools, so be sure to **update often**. 

### Requirements ###

* Python 2.6+ or 3.x

NOTE: Scipy and Numpy will be required for tools in the near future.

### How to use ###
Instructions for each tool are in the help notes.  
Just run a script with the '-h' option to see a full description of the parameters.
A list of all mixTAPE tools is at the bottom.

### Have a question?  Want a new tool? ###

Email me or visit my web page (see above).

### Citing mixTAPE ###
If you use my tool in published work, please include:
* the name "mixTAPE"
* the name of the specific tool (i.e. "count_sam") 
* the URL of this project (http://www.bitbucket.org/jbpease/mixtape)



### mixTAPE tracks ###
* clean_fasta.py: Clean up headers and sequences to prepare a FASTA file (for software with specific format requirements).
* count_sam.py: Count SAM file mapping flags and count each bit separately, with percentages and text descriptions of the flags.
* sam2seq: Convert SAM alignment file back to FASTA or FASTQ format.

