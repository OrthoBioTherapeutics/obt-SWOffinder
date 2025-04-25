[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10115908.svg)](https://doi.org/10.5281/zenodo.10115908)

# OBT-SWOffinder
SWOffinder is a method based on Smith-Waterman alignment to find all off-target sites up to some edit distance efficiently.
We implemented a trace-back method to find off-target sites under versatile criteria, such as separate limits on the number of insertions, deletions, and mismatches.


## Usage
1. Populate `guides.txt` with the sequences (**including the PAM sequence**) of each guide of interest.
2. Activate conda env `obt-swoff`
3. Run SWOffinder with target guides and ref genome (`run_swoffinder.sh` + `guides.txt` + ref genome)
4. Calculate gap-aware CFD score calculator R script (score_pam.R)
5. Optional: aggregate CFD scores for each guide (aggregate.py)
6. Use output of step 2 for gene annotation R script
7. Use output of step 4 for disease annotation python script
8. Optional: run combine_excel.py to combine csv files into excel sheets


## Installation
`conda env create -f obt-swoff.yaml`


# Original SWOffinder
1. First, you need to compile and build the SmithWatermanOffTarget package: `javac -d bin SmithWatermanOffTarget/*.java`
2. To search, execute the command line `java -cp bin SmithWatermanOffTarget.SmithWatermanOffTargetSearchAlign <Genome reference path> <sgRNA list path> <Output path prefix> <maxE> <maxM> <maxMB> <maxB> <Threads> <Best in-a-window> <Best in-a-window size> <PAM> <Allow PAM edits>`

The arguments specification:

1. **Genome reference path**: The path of the FASTA file containing the Genome reference.
2. **sgRNA list**: The path of a text file containing a list of sgRNAs with their PAM (see sgRNAs.txt file for example).
3. **Output path prefix**: The path for the output files with file prefix. files are saved in "CSV" format.
4. **maxE**: Max edits allowed (integer).
5. **maxM**: Max mismatches allowed without bulges (integer).
6. **maxMB**: Max mismatches allowed with bulges (integer).
7. **maxB**: Max bulges allowed (integer).
8. **Threads**: The number of threads to use for the run (integer). 
9. **Best in-a-window**: Flag whether to choose the best off-target site in a window or not (true or false).
10. **Best in-a-window size**: The window size for choosing the best in a window (integer). Please insert even if **Best in-a-window** is false.
11. **PAM**: The PAM sequence (for example, NGG).
12. **Allow PAM edits**: Flag whether to allow PAM edits or not (true or false).


