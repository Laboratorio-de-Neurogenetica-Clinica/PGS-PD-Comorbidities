# PGS-PD-Comorbidities
Pipeline with methodology to replicate manuscript "Heterogeneous contribution of Polygenic Scores to comorbidities presentation in Parkinson's disease."
## Files
0. *EPIgwas* (NOT INCLUDED): genotype in PLINK format (mygwas.ped and mygwas.map) and sample-wise microarray intensity files.
1. *over chain*: chain file to liftover GWAS SS.
3. *refseq.genes.bed*: All RefSeq genes.
4. *geneset.list*: genes of interest.
5. *repetitive.regions.bed*: list of genomic repetetitive regions (e.g. centromere, telomere).
6. Expected-output/: Folder containg all intermediate files. 

## Dependencies
1. [PLINK 1.9](https://www.cog-genomics.org/plink/1.9/)
2. [PLINK 2](https://www.cog-genomics.org/plink/2.0/)
3. [Liftover](https://genome.sph.umich.edu/wiki/LiftOver)
4. 
5. [PennCNV](http://penncnv.openbioinformatics.org/en/latest/)
6. [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html)
7. [Bedtools](https://bedtools.readthedocs.io/en/latest/)
8. [R](https://www.r-project.org/)
9. [Python](https://www.python.org/)

## Directory structure
PGS-PD-Comorbidities
├── PRS_EPI_all_2022
│   ├── 1_SS_GWAS_EPI_all_2022
|   |   ├──     
│   ├── FILTERED.test.bim
│   ├── FILTERED.test.fam
│   └── FILTERED.test.log
├── 1_GWAS_liftover.sh
└── 2_GWAS_QC.sh

## Step 1. Quality controls for GWAS summary statistics

### 1. Downlaoad GWAS EPI 
```
GWAS_SS_PATH=PRS_EPI_all_2022/1_SS_GWAS_EPI_all_2022
wget -P $GWAS_SS_PATH -c https://www.epigad.org/download/final_sumstats.zip
unzip -d $GWAS_SS_PATH $GWAS_SS_PATH/final_sumstats.zip 

```
### 2 .Adjust GWAS SS
GWAS SS file need to be formated in this way
A Space delimieted file with a header like this "CHR BP MarkerName Allele1 Allele2 Freq1 Beta P-value"
1. CHR
2. BP
3. MarkerName (If not present, fill with NA)
4. Allele1 (Effect allele gwas)
5. Allele2 (Non-effect allele)
6. Freq1 (If not present, fill with NA)
7. Beta
8. P-value
```
GWAS_SS_PATH=PRS_EPI_all_2022/1_SS_GWAS_EPI_all_2022
{
echo "CHR BP MarkerName Allele1 Allele2 Freq1 Beta P-value"
cat $GWAS_SS_PATH/ILAE3_TRANS_all_epilepsy_final.tbl | sed '1d' | cut -f1,2,3,4,5,6,12,10  
} | sed 's/\t/ /g' > $GWAS_SS_PATH/GWAS_SS_EPI_all_2022_RAW_hg19.tsv
```
### 3. Liftover GWAS SS to GRCh38
```
./1_GWAS__liftover.sh -c EPI_all_2022 -X true -C hg19ToHg38.over.chain.gz
```
### 4. Perform QC steps in GWAS SS 
```
./2_GWAS_QC.sh -c EPI_all_2022 -X true
```

#### 5. Perform QC in UKBB data 
##### filter BGEN files of UKB by allele frequency
##### run in jupyter lab instance of mem2_ssd2_v2_x32
```
for chr in {1..22} X ; do 
	cat ukb21007_c${chr}_b0_v1.sample | sed 's/0 0 0 0/0 0 0 D/g'> sample_chr${chr}.sample
	plink2  --bgen ukb21007_c${chr}_b0_v1.bgen ref-first \
			--sample sample_chr${chr}.sample \
			--oxford-single-chr ${chr}  \
			--set-all-var-ids @:#:\$r:\$a --new-id-max-allele-len 100 missing \
			--maf 0.01 \
			--make-pfile --out ukb_impTopMed_chr${chr}
	grep -v '#' ukb_impTopMed_chr${chr}.pvar > ukb_impTopMed_chr${chr}_noHeader.pvar
done
```
#### QC UK biobank 
```
./3_QC_UKBB.r EPI_all_2022
```

##### extract common positions between UKBB and Epilepsy GWAS 
```
for chr in {1..22} X ; do
	CMRB=EPI_all_2022
	plink2  --pfile ukb_impTopMed_chr${chr} --exclude ${CMRB}_chr${chr}_snplist_pvar.mismatch --ref-allele force ${CMRB}_chr${chr}_pvar_updated --make-pfile --out ukb_imp_chr${chr}_swapped --memory 120000
	plink2  --pfile ukb_imp_chr${chr}_swapped --set-all-var-ids @:#:\$r:\$a --new-id-max-allele-len 100 missing --make-bed --out ukb_imp_chr${chr}_recoded --memory 120000
	rm ukb_impTopMed_chr${chr}.*
done
```
### 6. PRS Calcualtions Clumping + thresholding
#### Clumping
```
for chr in 1 2 21 22 ; do
	plink   --bfile ukb_imp_chr${chr}_recoded --clump-p1 1 --clump-r2 0.1 --clump-kb 500 --clump final_SS_QC.tsv --clump-snp-field MarkerName --clump-field Pvalue --out clumped_SS_chr$chr --memory 60000
	awk 'NR!=1{print $3}' clumped_SS_chr$chr.clumped | sed '/^$/d' >  SS_TEMP_chr${chr}.valid.snp
	plink2  --bfile ukb_imp_chr${chr}_recoded --extract SS_TEMP_chr${chr}.valid.snp --make-pfile --out clumped_SS_chr$chr --memory 60000
done 
```
##### Merge clumped data by chromosome in a unique file
```
ls | grep ".pvar" | grep "clumped" > TEMP1.txt ; awk -F'.' '{print $1 }' TEMP1.txt > merge_list.chrALL.txt
rm TEMP1.txt
plink2 --pmerge-list merge_list.chrALL.txt --make-bed --out clumped_SS_chr.all
for chr in {1..22} X; do
	mv SS_TEMP_chr${chr}.valid.snp SS_chr${chr}.valid.snp
done   
```
#### Define SNPs and their Pvalue 
```
echo -n "" > SNP.pvalue
for chr in {1..22} X ; do
	cat final_SS_QC.tsv | awk -F' ' '{print $1,$9}' | grep -wFf SS_chr${chr}.valid.snp  >> SNP.pvalue
done
echo "0.00001 0 0.00001" >> range_list
echo "0.0001 0 0.0001" >> range_list
echo "0.001 0 0.001" >> range_list
echo "0.01 0 0.01" >> range_list
echo "0.05 0 0.05" >> range_list
echo "0.1 0 0.1" >> range_list
echo "0.2 0 0.2" >> range_list
echo "0.3 0 0.3" >> range_list
echo "0.4 0 0.4" >> range_list
echo "0.5 0 0.5" >> range_list
echo "0.8 0 0.8" >> range_list
echo "1 0 1" >> range_list 
```
#### Take valid SNPs from GWAS SS to use it in PRS calcualtions
```
echo -n "" > clumped_final_SS_common_SNPS_all.txt
for chr in {1..22} ; do
	grep -wFf SS_chr${chr}.valid.snp final_SS_QC.tsv >> clumped_final_SS_common_SNPS_all.txt
done 
```
#### Calculate PRS for all Pvalue trehshold
```
plink --bfile clumped_SS_chr.all --score clumped_final_SS_common_SNPS_all.txt 1 5 8 header --q-score-range range_list SNP.pvalue --out SS_PLINK1.9_PRS --memory 60000
cat clumped_final_SS_common_SNPS_all.txt | grep -v "nan" > TEMP
plink2 --bfile clumped_SS_chr.all --score TEMP 1 5 8 header --q-score-range range_list SNP.pvalue --out SS_PLINK2_PRS --memory 60000 
```

