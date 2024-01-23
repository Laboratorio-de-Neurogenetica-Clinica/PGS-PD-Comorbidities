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
1. [PLINK 1.9](https://www.cog-genomics.org/plink/1.9/)
2. [PLINK 2](https://www.cog-genomics.org/plink/2.0/)
3. [EIGENSTRAT](https://github.com/DReichLab/EIG)
4. [PennCNV](http://penncnv.openbioinformatics.org/en/latest/)
5. [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html)
6. [Bedtools](https://bedtools.readthedocs.io/en/latest/)
7. [R](https://www.r-project.org/)
8. [Python](https://www.python.org/)
## Step 1. Quality controls for GWAS summary statistics

### Downlaoad GWAS EPI 
```
wget -c https://www.epigad.org/download/final_sumstats.zip
unzip final_sumstats.zip
```
### Adjust GWAS SS
GWAS SS file need to be formated in this way
1. Space delimieted file
2. Header: CHR BP MarkerName Allele1 Allele2 Freq1 Beta P-value
   a. CHR
   b. BP
   c. MarkerName (If not present, fill with NA)
   d. Allele1 (Effect allele gwas)
   e. Allele2 (Non-effect allele)
   f. Freq1 (If not present, fill with NA)
   g. Beta
   h. P-value
```
{
echo "CHR BP MarkerName Allele1 Allele2 Freq1 Beta P-value"
cat ILAE3_TRANS_all_epilepsy_final.tbl | sed '1d' | cut -f1,2,3,4,5,6,12,10  
} | sed 's/\t/ /g' > GWAS_SS_EPI_all_2022_RAW_hg19.tsv
```
