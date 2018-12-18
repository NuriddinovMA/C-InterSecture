## C-InterSecture
**C-InterSecture** (*C*omputional tool for *InterS*pecies analysis of genome archit*ecture*) pipline is python 2.7 based utilits to cross-species comparison of Hi-C map. C-InterSecture was designed to liftover contacts between species, compare 3-dimensional organization of defined genomic regions, such as TADs, and analyze statistically of individual contact frequencies.
 
# Prerequisites:
- numpy >= 1.9.0
- scipy >= 1.1.0
- [genome module](https://mirnylab.bitbucket.io/hiclib/_modules/mirnylib/genome.html) from mirnylib package with dependencies (it's minimal required module is included)

# Installation
There is no need for installation.

# Using C-InterSecture
C-InterSecture pipeline involves three step: data preprocessing, contact liftovering and, finally, visualization. Scripts for each step are placed in special folder. 

### Data preprocessing
The preprocessing scripts are located in "0_preprocessing" folder.
The preprocessing includes a contacts filtration, statistical analysis and equalization. This step requires genome files, raw and normalized contact matrices.
To filtrate contacts, the pipeline requires genome fasta-files. Using unmappedBins.py, it calculates a coverage by N-bases of each bin.
`python unmappedBins.py < unmapped.ini`
