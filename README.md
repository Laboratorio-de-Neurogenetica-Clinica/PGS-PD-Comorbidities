# PGS-PD-Comorbidities
Pipeline with methodology to replicate manuscript "Heterogeneous contribution of Polygenic Scores to comorbidities presentation in Parkinson's disease."
This tutorial allow to replicate the analysis of PGS of Epilepsy in Parkinson disease Patients from UK biobank. 
## Files
0. *EPI_GWAS* (NOT INCLUDED): Summary statistics from the last meatanalysis of epilepsy .
1. *over chain*: chain file to liftover GWAS SS.
2. *geneset.list*: genes of interest.
3. *repetitive.regions.bed*: list of genomic repetetitive regions (e.g. centromere, telomere).
4. PRS_EPI_all_2022/: Folder containg all intermediate files. 

## Dependencies
1. [PLINK 1.9](https://www.cog-genomics.org/plink/1.9/)
2. [PLINK 2](https://www.cog-genomics.org/plink/2.0/)
3. [Liftover](https://genome.sph.umich.edu/wiki/LiftOver)
4. [R](https://www.r-project.org/)

## Directory structure
	PGS-PD-Comorbidities
	├── PRS_EPI_all_2022
	│   ├── 1_SS_GWAS_EPI_all_2022
	|   ├── 2_QC_target_data     
	│   ├── 3_Clump_PRS
	│   ├── 4_Tresholding
	├── 1_GWAS_liftover.sh
	├── 2_GWAS_QC.sh
 	└── 3_QC_UKBB.r

## Step 1: Quality controls for GWAS summary statistics
For a conplete detail of Quality controls (QC) procedures refers to [Choi *et al*, 2020](https://www.nature.com/articles/s41596-020-0353-1). Summary statistics from Epilepsy can be downloaded from [here](https://www.epigad.org/download/final_sumstats.zip), or with the code below. GWAS of epilepsy can be found in [ILAE *et al*, 2022](https://www.nature.com/articles/s41588-023-01485-w#article-info).  

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

In order to run this tutorial with other comorbidities, please follow the recommendations above, this is only an example of processing for epilepsy GWAS.
```
GWAS_SS_PATH=PRS_EPI_all_2022/1_SS_GWAS_EPI_all_2022
{
echo "CHR BP MarkerName Allele1 Allele2 Freq1 Beta P-value"
cat $GWAS_SS_PATH/ILAE3_TRANS_all_epilepsy_final.tbl | sed '1d' | cut -f1,2,3,4,5,6,12,10  
} | sed 's/\t/ /g' > $GWAS_SS_PATH/GWAS_SS_EPI_all_2022_RAW_hg19.tsv
```
### 3. Liftover GWAS SS to GRCh38
Most of the GWAS are published in GRCh37 version. And the imputed data from UK Biobank is in GRCh38, so it is necessary to update the GWAS SS potions to GRCh48 version with Liftover. The following bash script can be used to generate the GWAS SS in GRCh38.
```
#!/usr/bin/bash
set -e
#### ARGUMENTS #####
# -c Comorbidity GWAS used
# -d Directory
# -X Indicate if GWAS SS has X chromosome
# -C Path to over cahin file to liftover


while getopts "c:d:X:C:" flag ;
do
    case "${flag}" in
        c) CMRB=${OPTARG};;
        d) WORK_START=${OPTARG};;
        X) xchrms=${OPTARG};;
        C) CHAIN=${OPTARG};;
    esac
done
if [[ -z $WORK_START  ]] ;then WORK_START=$(pwd) ; fi
test() { 
    cd $WORK_START
}
test || {
echo "Enter a valid Directory" 1>&2
false
}
if [[ -z $CMRB ]]; then echo "Indicate some parameter for Comorbilities (-c)" ; exit 0 ; fi
test() { 
    cd PRS_$CMRB
}
test || {
echo "Enter a valid Comorbidity" 1>&2
false
}
if [[ -z $xchrms ]] ;then echo "-X not present: indicate if GWAS SS have X chromosome (true or false)"; exit 0 ; fi
if [[ $xchrms == true ]] ;then 
    CHRS=$(echo 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X   ) ;### variable util para loops
elif [[ $xchrms == false ]] ; then
    CHRS=$(echo 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22  ) ;### variable util para loops
else 
    echo "invalid -X valid options are: true or false"
fi

if [[ -z $CHAIN ]] ;then echo "-C not present: indicate path to hg19ToHg38.over.chain.gz "; exit 0 ; fi


WORKD=$WORK_START
echo Working in $WORKD 

### generate .bed file to liftover
cd $WORKD/1_SS_GWAS_${CMRB}

printf "Formating GWAS SS to match hg19 to liftover: "
cat GWAS_SS_${CMRB}_RAW_hg19.tsv | awk -F' ' '{if (NR>1) print "chr"$1 "\t" ($2-1) "\t" $2 "\t" $3","$4","$5","$6","$7","$8 }' | sed 's/chr23/chrX/g'  > GWAS_SS_hg19.bed
printf "DONE\n"

printf "Liftover to hg38: "
liftOver GWAS_SS_hg19.bed $CHAIN GWAS_SS_hg38.bed GWAS_SS_unlifted.bed
printf "DONE\n"

printf "Formating hg38 to GWAS SS pipeline PRS: "
cat GWAS_SS_hg38.bed | sed 's/,/\t/g' | sed 's/\t/ /g' | awk -F' ' '{printf $1; for (i=3; i <= NF; i++) printf FS$i; print NL }' | sed 's/^chr//g' > GWAS_SS_${CMRB}_RAW_hg38.tsv 
printf "DONE\n"
```
```
./1_GWAS__liftover.sh -c EPI_all_2022 -X true -C hg19ToHg38.over.chain.gz
```
### 4. Perform QC steps in GWAS SS 
The GWAS SS quality controls consist of eliminating duplicate and ambiguous SNPs. The following bash script allows to generate the final GWAS SS.
```
#!/usr/bin/bash
set -e
#### ARGUMENTS #####
# -c Comorbidity GWAS used
# -d Directory
# -X Indicate if GWAS SS has X chromosome

while getopts "c:d:X:" flag ;
do
    case "${flag}" in
        c) CMRB=${OPTARG};;
        d) WORK_START=${OPTARG};;
        X) xchrms=${OPTARG};;

    esac
done
if [[ -z $WORK_START  ]] ;then WORK_START=$(pwd) ; fi
test() { 
    cd $WORK_START
}
test || {
echo "Enter a valid Directory" 1>&2
false
}
if [[ -z $CMRB ]]; then echo "Indicate some parameter for Comorbilities (-c)" ; exit 0 ; fi
test() { 
    cd PRS_$CMRB
}
test || {
echo "Enter a valid Comorbidity" 1>&2
false
}
if [[ -z $xchrms  ]] ;then echo "-X not present: indicate if GWAS SS have X chromosome (true or false)"; exit 0 ; fi


WORKD=$(pwd)
echo Working in $WORKD 

################ set Variables ##############################
if [[ $xchrms == true ]] ;then 
    CHRS=$(echo 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X   ) ;### variable util para loops
elif [[ $xchrms == false ]] ; then
    CHRS=$(echo 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22  ) ;### variable util para loops
else 
    echo "invalid -X valid options are: true or false"
fi


########################################## QC GWAS SS
echo -e "\nQC base data, SS_GWAS_$CMRB\n"
cd $WORKD/1_SS_GWAS_$CMRB

############################ QC SS Comorbidities
##### Format
#header CHR BP MarkerName Allele1 Allele2 Freq1 Beta P-value
# 1 CHR
# 2 BP
# 3 MarkerName (If not present, fill with NA)
# 4 Allele1 (effect allele)
# 5 Allele2 (non-effect allele)
# 6 Freq1 (If not present, fill with NA)
# 7 Beta
# 8 P-value

if [[ ! -f final_${CMRB}_QC.tsv ]] ; then    
    cat GWAS_SS_${CMRB}_RAW_hg38.tsv |  sort  -nk1,1 -k2,2 > ${CMRB}_sorted.tsv
    # filter
    cat ${CMRB}_sorted.tsv | awk '{if ($6 >= 0.01) {print $0}}' > ${CMRB}_filtered.tsv

    #Alleles to upercase
    cat ${CMRB}_filtered.tsv | awk '{print $1,$2,$3,toupper($4) ,toupper($5) ,$6 ,$7 ,$8 }' > ${CMRB}_filtered_touper.tsv

    #recode SNPS
    cat ${CMRB}_filtered_touper.tsv | awk -F' ' '{print " " $1":"$2":"$4":"$5 FS $3 FS $1 FS $2 FS $4 FS $5 FS $6 FS $7 FS $8}' > ${CMRB}_recode.tsv
    ## remove duplciated
    cat ${CMRB}_recode.tsv | awk '{seen[$1$2]++; if(seen[$1$2]==1){ print}}' > ${CMRB}_nodup.tsv

    # Remove ambiguous SNPs
    awk '!( ($5=="A" && $6=="T") || \
        ($5=="T" && $6=="A") || \
        ($5=="G" && $6=="C") || \
        ($5=="C" && $6=="G")) {print}' ${CMRB}_nodup.tsv > ${CMRB}_no_ambiguos.tsv

    cat ${CMRB}_no_ambiguos.tsv | sed 's/ 23:/X:/g' | sed 's/ 23 / X /g' > final_${CMRB}_QC.tsv 
else 
    echo QC of Base data already performed, skipping
fi
```
```
./2_GWAS_QC.sh -c EPI_all_2022 -X true
```

## Step 2: Perform QC in UKBB data 
The reference and alternative alleles of the UK biobank must match those reported in the GWAS SS, so we must generate a quality control of the UK biobank in order to work with the SNPs in a unified way.
### 1. Filter BGEN files of UKB by allele frequency
Due to the magnitude of the imputed UKBB data we first filtered by allelic frequency greater than 1% and then performed the corresponding quality controls on smaller files.
The next step was performed in the UK Biobank Research Analysis Platform (RAP), due to the availability of the data, we cannot provide the BGEN file but it can be replicated with the following command in RAP, using Swiss army knife tool or in an instance of JupyterLab, the requested virtual machine was mem2_ssd2_v2_x32
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
### 2. QC UK biobank 
Once the BGEN files have been filtered, it is necessary to identify the common SNPs between the UK bibank and the GWAS SS. To do this, using the .pvar file (provided in this tutorial) and the GWAS QC file generated in step 1, we identify the common SNPs and generate a list to update the .pvar file with the alleles reported in the GWAS SS. The following R script is the one in charge of what is described above
```
#!/usr/bin/env Rscript
####USAGE######
# ./3_QC_UKB.r [Comorbidity]
#
#
#
#
args <- commandArgs(trailingOnly = TRUE)

library(data.table)
library(magrittr)
library(stringr)


# Function for finding the complementary allele
complement <- function(x) {
    switch (
        x,
        "A" = "T",
        "C" = "G",
        "T" = "A",
        "G" = "C",
        return(NA)
    )
}
CMRB=args[1]
# Read in summary statistic data (require data.table v1.12.0+)
print("Reading SS file")
height <- fread(paste("final_",CMRB,"_QC.tsv",sep="")) %>%
    setnames(., colnames(.), c("ID", "Marker", "CHR", "BP", "A1", "A2","Freq1","Beta","P-value")) %>%
    # And immediately change the alleles to upper cases
    .[,c("A1","A2"):=list(toupper(A1), toupper(A2))] 
str(height)
for (chr in c(1:22,"X")) {
  # Read in bim file  A1=minor A2=major
    print(paste("Reading BIM file of chromosome",chr,sep=" "))
    bim <- fread(paste("ukb_impTopMed_chr",chr,"_noHeader.pvar",sep=""),colClasses=c(V1="character")) %>%
        # Note: . represents the output from previous step
        # The syntax here means, setnames of the data read from
        # the bim file, and replace the original column names by 
        # the new names
       setnames(., colnames(.), c("CHR", "BP", "ID", "A1", "A2")) %>%
        # And immediately change the alleles to upper cases
        .[,c("A1","A2"):=list(toupper(A1), toupper(A2))]
    bim$Position <- str_c(bim$CHR,bim$BP,sep=":")
    colnames(bim) <- c("CHR", "BP", "ID", "B.A1", "B.A2","Position") 
    
    #filter SS to merge
    SS_filtered <- SS[SS$CHR == chr]
    SS_filtered$Position <- str_c(SS_filtered$CHR,SS_filtered$BP,sep=":") 
    # Merge summary statistic with target
    print("Processing")
    info <- merge(bim, SS_filtered, 
              by = c("CHR","BP","Position"),
              suffixes = c("_bim","_SS") )
    
    # Get SNPs that have the same alleles across base and target
    info.match <- subset(info, A1 == B.A1 & A2 == B.A2)
    
    # Identify SNPs that are complementary between base and target
    info$C.A1 <- sapply(info$B.A1, complement)
    info$C.A2 <- sapply(info$B.A2, complement)
    info.complement <- subset(info, A1 == C.A1 & A2 == C.A2)
    # Update the complementary alleles in the bim file
    # This allow us to match the allele in subsequent analysis
    complement.snps <- bim$ID %in% info.complement$ID_bim
    bim[complement.snps,]$B.A1 <-
        sapply(bim[complement.snps,]$B.A1, complement)
    bim[complement.snps,]$B.A2 <-
        sapply(bim[complement.snps,]$B.A2, complement)
    
    # identify SNPs that need recoding
    info.recode <- subset(info, A1 == B.A2 & A2 == B.A1)
    # Update the recode SNPs
    recode.snps <- bim$ID %in% info.recode$ID_bim
    tmp <- bim[recode.snps,]$B.A1
    bim[recode.snps,]$B.A1 <- bim[recode.snps,]$B.A2
    bim[recode.snps,]$B.A2 <- tmp
    
    # identify SNPs that need recoding & complement
    info.crecode <- subset(info, A1 == C.A2 & A2 == C.A1)
    # Update the recode + strand flip SNPs
    com.snps <- bim$ID %in% info.crecode$ID_bim
    tmp <- bim[com.snps,]$B.A1
    bim[com.snps,]$B.A1 <- as.character(sapply(bim[com.snps,]$B.A2, complement))
    bim[com.snps,]$B.A2 <- as.character(sapply(tmp, complement))
    
    # Output updated bim file
    print(paste("Writing pvar File updated of chromosome",chr,sep=" "))
    fwrite(bim[,c("ID", "B.A1")], paste(CMRB,"_chr",chr,"_pvar_updated",sep=""),quote = F,row.names = F , col.names=F, sep="\t")
    
    # Output mismatch file
    print(paste("Writing QCded SNPs of chromosome",chr,sep=" "))
    mismatch <-
        bim$ID[!(bim$ID %in% info.match$ID_bim |
                  bim$ID %in% info.complement$ID_bim | 
                  bim$ID %in% info.recode$ID_bim |
                  bim$ID %in% info.crecode$ID_bim)]
    write.table(mismatch, paste(CMRB,"_chr",chr,"_snplist_pvar.mismatch",sep=""), quote=F, row.names=F, col.names=F)
    
}
```
```
./3_QC_UKBB.r EPI_all_2022
```

### 3. Extract common positions between UKBB and Epilepsy GWAS 
With the files generated in the previous step, we proceed to filter the UKbiobank imputed data and update the reported alleles. Due to the availability of the data, we cannot provide the filtered imputed data, but it can be replicated with the following command in RAP, using Swiss army knife tool or in an instance of JupyterLab, the requested virtual machine was mem2_ssd2_v2_x32
```
for chr in {1..22} X ; do
	CMRB=EPI_all_2022
	plink2  --pfile ukb_impTopMed_chr${chr} --exclude ${CMRB}_chr${chr}_snplist_pvar.mismatch --ref-allele force ${CMRB}_chr${chr}_pvar_updated --make-pfile --out ukb_imp_chr${chr}_swapped --memory 120000
	plink2  --pfile ukb_imp_chr${chr}_swapped --set-all-var-ids @:#:\$r:\$a --new-id-max-allele-len 100 missing --make-bed --out ukb_imp_chr${chr}_recoded --memory 120000
	rm ukb_impTopMed_chr${chr}.*
done
```
## Step 3: PRS Calcualtions Clumping
With the common SNPs we proceed to clumping, which consists of selecting independent variants in the haplotype blocks, considering the GWAS SS pvalue. We will use the clumped data to calculate the PRS using different pvalue thresholds.  
### 1. Clumping
The following code needs to be done in RAP, which consists of taking the QC imputed data and performing the clumping and generating a new pfiles file.
```
for chr in 1 2 21 22 ; do
	plink   --bfile ukb_imp_chr${chr}_recoded --clump-p1 1 --clump-r2 0.1 --clump-kb 500 --clump final_SS_QC.tsv --clump-snp-field MarkerName --clump-field Pvalue --out clumped_SS_chr$chr --memory 60000
	awk 'NR!=1{print $3}' clumped_SS_chr$chr.clumped | sed '/^$/d' >  SS_TEMP_chr${chr}.valid.snp
	plink2  --bfile ukb_imp_chr${chr}_recoded --extract SS_TEMP_chr${chr}.valid.snp --make-pfile --out clumped_SS_chr$chr --memory 60000
done 
```
#### 2. Merge clumped data by chromosome in a unique file
The following code needs to be done in RAP, which allows to merge the data by chromosomes, to generate a single file on which to calculate the PRS.
```
ls | grep ".pvar" | grep "clumped" > TEMP1.txt ; awk -F'.' '{print $1 }' TEMP1.txt > merge_list.chrALL.txt
rm TEMP1.txt
plink2 --pmerge-list merge_list.chrALL.txt --make-bed --out clumped_SS_chr.all
for chr in {1..22} X; do
	mv SS_TEMP_chr${chr}.valid.snp SS_chr${chr}.valid.snp
done   
```
### 3. Define SNPs and their Pvalue 
The following code needs to be done in RAP, it takes the GWAS SS and extracts the SNPs ID along with their pvalue, to be used to define the snps that meet the pvalue treshold. additionally it generates the list of the pvalues included in each treshold AND filters the GWAS SS with the valid SNPs obtained in step 2.
```
#Generate SNP.pvalue
echo -n "" > SNP.pvalue
for chr in {1..22} X ; do
	cat final_SS_QC.tsv | awk -F' ' '{print $1,$9}' | grep -wFf SS_chr${chr}.valid.snp  >> SNP.pvalue
done

#Generate list of treshold used
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

#Filter GWAS SS to use only valid SNPs	
echo -n "" > clumped_final_SS_common_SNPS_all.txt
for chr in {1..22} ; do
	grep -wFf SS_chr${chr}.valid.snp final_SS_QC.tsv >> clumped_final_SS_common_SNPS_all.txt
done 
```
### 5. Calculate PRS for all Pvalue trehshold
The following code needs to be done in RAP, it takes the files generated above to calculate the PRS for each of the pvalue treshold defined above
```
plink --bfile clumped_SS_chr.all --score clumped_final_SS_common_SNPS_all.txt 1 5 8 header --q-score-range range_list SNP.pvalue --out SS_PLINK1.9_PRS --memory 60000
```

