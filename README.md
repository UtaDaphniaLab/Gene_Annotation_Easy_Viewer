# KEGG Annotation Easy Viewer (KAEV)
KAEV is a tool to help visualize BLAST results after using KEGG Automatic Annotation Server (KAAS) to annotate a region of DNA. The software focuses on organizing the gene annotation data obtained from KAAS in a gene-centric view. The goal of this software is to help researchers predict the function of genes belonging to non-model organisms.

## Starting Up

### Prerequisites
This code requires python 3 to run. It is recommended to use the most up-to-date version of python, which can be downloaded [here](https://www.python.org/downloads/).

The folder containing the KAEV software as well as input and output examples can found [here](https://github.com/UtaDaphniaLab/KEGG_Annotation_Easy_Viewer).

### Installation
KAEV requires no installation. If python 3 is properly installed and set up, then the user only needs to run KAEV.py to start the code. It is recommended to move KAEV.py and run it inside the same directory that the input file is in to make specifying relative paths more easy.

## Usage

### Input files
In the kegg_annotation_easy_viewer folder, there is an example input file. The input should be taken from the KAAS results webpage for K-code assignment. It should be formatted in two columns with no header. The first column should contain the gene ID while the second column should contain the k-code that was assigned to the gene. Not all genes may have a k-code assignment.
```
gene1
gene2	K12829
gene3	K14963
gene4	K20672
gene5	K12184
gene6
gene7	K04958
gene8
gene9	K09075
gene10
```
### Data files
After the first time an input file has been processed, all the information extracted from the KEGG servers for that set of genes will be stored on data files (.dat) with the same name as the input files had. These data files can be loaded by KAEV at any time to generate tables without needing to extract the information from the KEGG servers anymore. The option to load a data file instead of entering an input file can be found on the first menu displayed by KAEV after running.

### Output files
KAEV will output tables in HTML files. Users may choose between a gene-centric table, a pathway-centric table, or both. Users may also apply filters to the data, so that tables will be generated with only genes that contain a target term in its name, definition, or list of linked pathways. Each item in the associated pathway column of the table is an embedded url that will take the user to the pathway map on KEGG. Genes on the pathway map that were present in the original input file's genome assembly will be displayed in green. Meanwhile, the target gene is displayed in red to be easily distinguishable.  
An example output file can be found in the kegg_annotation_easy_viewer folder. It was produced using the example input file, choosing not to apply any filters, and choosing to display both genes and pathway tables.

## Authors
**Trung Huynh** - *Initial work on code*  
**Sen Xu** - *Conceptualization*

## License
This project is licensed under the MIT license and is available for free. 
