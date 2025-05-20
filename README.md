## This is a sample molecular biology and bioinformatics workflow that is run from one Jupyter notebook.
- Automating the analysis of multiple Sanger sequencing files (limit of 50 set by NCBI).
- DNA Sanger sequencing data to a CSV file of the top Blast matches, XML files of 50 matches are also saved.

### ab1_blast_workflow_template.ipynb ### 
- can be used with your own sequencing files.

### ab1_blast_workflow_example.ipynb ### 
- is an example that has been run with 10 .ab1 sequencing files. Data for these are in the [seq_data](seq_data/old_ab1s/) folder and results are in the [python_outputs](python_outputs/old_ab1s/) folder. The subfolder for both is /old_ab1s.

## functions working in the background
- these are in [src folder](/src) in a file called ab1_helpers.py 

Tested using Python 3.13 in a Jupyter notebook in VSCode.

Required modules (BioPython, ipykernel for Jupyter notebooks) can be installed from the requirements.txt file. 



-Annette
