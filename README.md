# Gene Annotation Easy Viewer (GAEV)
GAEV is a tool to help visualize BLAST results after using KEGG Automatic Annotation Server (KAAS) to annotate a region of DNA. The software focuses on organizing the gene annotation data obtained from KAAS in a gene-centric view. The goal of this software is to help researchers predict the function of genes belonging to non-model organisms.

## Starting Up

### Prerequisites
This code requires python 3 to run. It is recommended to use the most up-to-date version of python, which can be downloaded [here](https://www.python.org/downloads/).

The folder containing the GAEV software as well as input and output examples can found [here](https://github.com/UtaDaphniaLab/KEGG_Annotation_Easy_Viewer).

### Installation
GAEV requires no installation. If python 3 is properly installed and set up, then the user only needs to run GAEV.py to start the code. It is recommended to move GAEV.py and run it inside the same directory that the input file is in to make specifying relative paths more easy.

## Usage

### Input files
In the gene_annotation_easy_viewer folder, there is an example input file. The input should be taken from the KAAS results webpage for K-code assignment. It should be formatted in two columns with no header. The first column should contain the gene ID while the second column should contain the k-code that was assigned to the gene. Not all genes may have a k-code assignment.
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
After the first time an input file has been processed, all the information extracted from the KEGG servers for that set of genes will be stored on data files (.dat) with the same name as the input files had. These data files can be loaded by GAEV at any time to generate tables without needing to extract the information from the KEGG servers anymore. The option to load a data file instead of entering an input file can be found on the first menu displayed by GAEV after running.

### Output files
GAEV will output tables in HTML files. Users may choose between a gene-centric table, a pathway-centric table, or both. Users may also apply filters to the data, so that tables will be generated with only genes that contain a target term in its name, definition, or list of linked pathways. Each item in the associated pathway column of the table is an embedded url that will take the user to the pathway map on KEGG. Genes on the pathway map that were present in the original input file's genome assembly will be displayed in green. Meanwhile, the target gene is displayed in red to be easily distinguishable.  
An example output file can be found in the gene_annotation_easy_viewer folder. It was produced using the example input file, choosing not to apply any filters, and choosing to display both genes and pathway tables.

### Color Pathways
This new function allows users color any number of gene lists on KEGG pathways with different colors. This is useful if users need to annotate distinct lists of genes of interests (ex. overexpressed and underexpressed genes). It is complete with UI that can be accessed from the first menu by selecting option 5. Requires a dat file from a previous GAEV run to function. This function takes lists of user defined gene ids (gene ids must match with the GAEV dat file) with each gene seperated by new lines. Each list of genes can be colored with a different hue. Common colors can be specified with the color name (ex. "red") or or with hex codes (ex. "#FF0000"). Output is an HTML file with a hyperlink to each pathway similar to the standard GAEV output. Example input and output files for this function have been included [here](https://github.com/UtaDaphniaLab/Gene_Annotation_Easy_Viewer/tree/master/gene_annotation_easy_viewer/Color%20Pathways%20Example%20Inputs%20Outputs).

## Authors
**Trung Huynh** - *Initial work on code*  
**Sen Xu** - *Conceptualization*

## License
This project is licensed under the MIT license and is available for free.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2549592.svg)](https://doi.org/10.5281/zenodo.2549592)



