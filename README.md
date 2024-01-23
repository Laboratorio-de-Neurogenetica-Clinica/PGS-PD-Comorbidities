# PGS-PD-Comorbidities
Pipeline with methodology to replicate manuscript "Heterogeneous contribution of Polygenic Scores to comorbidities presentation in Parkinson's disease."
## Files
0. *EPIgwas* (NOT INCLUDED): genotype in PLINK format (mygwas.ped and mygwas.map) and sample-wise microarray intensity files.  
1. *mysamples*: All samples included in the study.
3. *refseq.genes.bed*: All RefSeq genes.
4. *geneset.list*: genes of interest.
5. *repetitive.regions.bed*: list of genomic repetetitive regions (e.g. centromere, telomere).
6. Expected-output/: Folder containg all intermediate files. 

## Dependencies
1. [PLINK 1.9](https://www.cog-genomics.org/plink2)
2. [KING](http://people.virginia.edu/~wc9c/KING/)
3. [EIGENSTRAT](https://github.com/DReichLab/EIG)
4. [PennCNV](http://penncnv.openbioinformatics.org/en/latest/)
5. [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html)
6. [Bedtools](https://bedtools.readthedocs.io/en/latest/)
7. [R](https://www.r-project.org/)
8. [Python](https://www.python.org/)
## Step 1. Quality controls for GWAS summary statistics

