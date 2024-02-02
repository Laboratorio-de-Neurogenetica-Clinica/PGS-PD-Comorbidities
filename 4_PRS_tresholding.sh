#!/usr/bin/bash
eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
set -e

while getopts "c:d:X:t:" flag ;
do
    case "${flag}" in
        c) CMRB=${OPTARG};;
        d) WORK_START=${OPTARG};;
        X) xchrms=${OPTARG};;
        t) TRESH=${OPTARG};;

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

if [[ -z $TRESH ]] ; then 
    
    
    ###### QC Samples
    cd $WORKD/4_Thresholding
    METADATA=$WORKD/Metadata/ 

    #Extract IDS from psam
    cat $WORKD/3_Clump_PRS/clumped_SS_chr22.psam | grep -v '#'| awk -F' ' '{print $1}' | sort  > Phenotipo.Ids     
    
    ##### QC ancestria
    cat $METADATA/CLINICAL_COVAR_QC_participant.tsv | cut -f1,26-45 > $METADATA/22009_PCA.tab
    mkdir -p PCA_results     
    $WORKD/../5_PCA_only_ukb.r EPI_all_2022 ${METADATA}/PCA.tsv
    
    grep -wf keep_ancestry_eur.txt Phenotipo.Ids > Phenotipo_QC_EUR.Ids 
    
    #QC Genetic Sex
    cat $METADATA/CLINICAL_COVAR_QC_participant.tsv | cut -f1,20 > $METADATA/22001_Genetic_sex.tab
        fi 
        cut -f2 $METADATA/Genetic_sex.tab | paste  $METADATA/31_SEX.tab - > TEMP  
        cat TEMP | awk -F'\t' '{if ($2 != $3) print  }' | sed 's/\t$/\tNA/g' | grep -v "NA" | cut -f1 | grep -v "ID" > remove_sex_discordant.ids 
        rm TEMP 
    fi
    
    ##### QC RELACIONADOS
    if [[ ! -f remove_REL.ids ]] ; then
        cp $METADATA/ukb_rel_a81646_s488130.dat rel_IDs_all.txt
        if [[ ! -f $METADATA/COVARS_rel_ids ]] ;then
            cat $METADATA/CLINICAL_COVAR_QC_participant.tsv | awk -F'\t' '{print $1 FS $71 FS $76 FS $2}' | sed 's/\t\t/\t0\t/g' | sed 's/\t\t/\t0\t/g' | sed 's/[0-9]*-[0-9]*-[0-9]*/1/g'  > $METADATA/COVARS_rel_ids
        fi 
        cp $METADATA/COVARS_rel_ids COVARS
    
        cat rel_IDs_all.txt | awk -F' ' '{if ($5 > 0.0442 ) print $1 FS $2}' | grep -v "ID"   >   rel_IDs_all_filter.txt
        cut -d ' ' -f1 rel_IDs_all_filter.txt > rel_ID1.txt 
        cut -d ' ' -f2 rel_IDs_all_filter.txt > rel_ID2.txt 
    
        $WORKD/../Scripts/Rscripts/ukb_imp_rel_IDS.r
    
        cat REL_IDS_to_remove.txt | sort | uniq | grep -v "remove" > TEMP 
        sed -e "s/\r//" TEMP  > remove_REL.ids 
        rm TEMP rel_IDs_all.txt COVARS rel_IDs_all_filter.txt rel_ID1.txt rel_ID2.txt 
    fi
    #QC Overlapping
    #if [[ ! -f REL_overlap_to_remove.txt ]] ; then
    #    cd $WORKD/3_clump_PRS
    #    wget https://personal.broadinstitute.org/sripke/share_links/checksums_download/id_geno_checksum.v2
    #    chmod 755 id_geno_checksum.v2
    #    conda activate plink
    #    ./id_geno_checksum.v2 
    #fi
    
    #Generar fenotipo Final
    if [[ ! -f Phenotipo_QC_EUR_final.Ids ]] ; then
        cat Phenotipo_QC_EUR.Ids | grep -vwFf remove_REL.ids  | grep -vwFf remove_sex_discordant.ids   > Phenotipo_QC_EUR_final.Ids 
    fi
    
    ##### Traer covariables
    if [[ ! -f final.covar.txt ]] ;then 
        if [[ ! -f $METADATA/COVARS_6 ]] ;then
        cat $METADATA/CLINICAL_COVAR_QC_participant.tsv | awk -F'\t' '{print $3 FS $6 FS $4}'| paste $METADATA/COVARS_rel_ids - > $METADATA/COVARS_6 
        fi 
        cat $METADATA/COVARS_6 | cut -f1-5,7 | sed 's/\t\t/\tNA\t/g' | sed 's/\t$/\tNA/g' | grep -wFf Phenotipo_QC_EUR_final.Ids  > final.covar.txt
    fi
    
    
    ##### traer PCA
    if [[ ! -f final.PCA.eigenvec ]] ;then
        cat $METADATA/22009_PCA.tab | grep -wFf Phenotipo_QC_EUR_final.Ids | awk -F'\t' '{print $1 FS $1 FS $3 FS $4 FS $5 FS $6 FS $7 FS $8 FS $9}' > final.PCA.eigenvec
    fi
    
    #### Definir Casos y controles Diferentes configuraciones
    #
    if [[ -z $age ]]; then echo "PHENOTYPE DATA : (-e)  Default value for Age filter (all)" ; age=0 ; fi
    if [[ -z $ctrl ]]; then echo "PHENOTYPE DATA : (-C) Default value for filter only in controls (false)" ; ctrl=false ; fi
    if [[ -z $excl_focalizada ]]; then echo "PHENOTYPE DATA : (-F) Exclude FE in controls, default false" ; excl_focalizada=false ; fi
    if [[ -z $FAM  ]] ;then echo "FAM (-f) patologias o Sin_Filtrar" ; exit 0 ; fi
    echo $fam
    $WORKD/../Scripts/Bash_scripts/Capturar_diagnosticos.sh -c$CMRB -I$ICD10 -e$age -C$ctrl -F$excl_focalizada -f$FAM
    
    mkdir -p $WORKD/4_Thresholding/${CMRB}/${FAM}
    cd $WORKD/4_Thresholding/${CMRB}/${FAM}
    
    if [[ ! -f final.Phenotipo_QC.fam ]] ; then 
        awk '{print $1 FS $1}' ../../Phenotipo_QC_EUR_final.Ids > keep_Ids 
        cp $METADATA/Diagnosticos/${CMRB}/${FAM}/case_controls.txt .
        if [[ $CMRB == Altura_2022 ]] ; then
            sed -i 's/\t$/\t-9/g' case_controls.txt   
        else
            sed -i 's/case/2/g' case_controls.txt 
            sed -i 's/control/1/g' case_controls.txt
        fi
        conda activate plink2

        printf "#FID\tIID\tPHENO1\n" | cat - case_controls.txt > case_controls_header.txt
        plink2 --psam $WORKD/4_Thresholding/clumped_SS_chr1.psam \
             --pheno case_controls_header.txt  \
             --keep keep_Ids \
             --make-just-psam --out final.Phenotipo_QC

        cat final.Phenotipo_QC.psam | sed 's/control/1/g' | 
            sed 's/case/2/g' | sed 's/NONE/-9/g' | 
            awk -F'\t' '{if (NR != 1) print $1 "\t" $2 "\t0\t0\t" $3 "\t" $4  }'  > final.Phenotipo_QC.fam    
        rm case_controls.txt  keep_Ids 
    fi
    echo $SOFTWARE
    
    cd $WORKD/4_Thresholding
    #### find best PRS 
    if [[ -z $SOFTWARE ]]; then echo "Software: (-S) PRS from PLINK1.9 or PLINK2, default false" ; echo $SOFTWARE ; exit 0 ; fi
    echo $SOFTWARE
    echo $CMRB

    ## define thresholds to use
    ls ${CMRB} | grep "SS_${SOFTWARE_PRS}" | grep $SUFFIX_SOFT | sed 's/PLINK1.9/PLINK19/g'| sed 's/0\./0,/g' | awk -F'.' '{print $2}' | sed 's/0,/0\./g'  | sort | sed 's/\n/\t/g' > ${CMRB}/used_threshold_${SOFTWARE}.txt
    
    for i in $(cat ${CMRB}/used_threshold_${SOFTWARE}.txt) ;
    do
        cp ${CMRB}/SS_${SOFTWARE}_PRS.$i.${SUFFIX_SOFT} ${CMRB}/${CMRB}_${SOFTWARE}_PRS.$i.${SUFFIX_SOFT} ;
    done 
        if [[ $CMRB == Altura_2022 ]] ;then
            echo "altura"
            $WORKD/../Scripts/Rscripts/CMRB_Best_PRS_altura.R $FAM $CMRB $SOFTWARE $SUFFIX_SOFT
        else 
            echo "other CMRBs"
            $WORKD/../Scripts/Rscripts/CMRB_Best_PRS.R $FAM $CMRB $SOFTWARE $SUFFIX_SOFT
        fi
        exit 0 ;  #calculate for all treshold find de best treshold
fi
echo "Using Defined Threshold"
echo $TRESH
cd $WORKD/4_Thresholding
echo $CMRB
if [[ ! -f ${CMRB}/${FAM}/${CMRB}_${TRESH}_${SOFTWARE}_PRS.covar.txt ]] ;then
    if [[ $CMRB == Altura_2022 ]] ;then
        $WORKD/../Scripts/Rscripts/CMRB_Best_PRS_altura.R $FAM $CMRB $SOFTWARE $SUFFIX_SOFT $TRESH
    else 
        $WORKD/../Scripts/Rscripts/CMRB_Best_PRS.R $FAM $CMRB $SOFTWARE $SUFFIX_SOFT $TRESH #save data for previous defined treshold
    fi
fi


## Analyze and plot data 
if [[ $CMRB == Altura_2022 ]] ; then
    echo DONE
else
    $WORKD/../Scripts/Rscripts/CMRB_PRS_PD_EPI.r $FAM $CMRB $SOFTWARE $SUFFIX_SOFT $TRESH
fi



exit 0




./Best_PRS.sh -c EPI_all_2022 -I G40 -X true -S PLINK1.9 -f patologias 
./Best_PRS.sh -c EPI_all_2022 -I G40 -X true -S PLINK2 -f patologias 
./Best_PRS.sh -c EPI_all_2022 -I G40 -X true -S PLINK1.9 -f patologias -t 0.05
./Best_PRS.sh -c EPI_all_2022 -I G40 -X true -S PLINK2 -f patologias -t 0.05
./Best_PRS.sh -c EPI_all_2022 -I G40 -X true -S PLINK1.9 -f patologias -t 0.5
./Best_PRS.sh -c EPI_all_2022 -I G40 -X true -S PLINK2 -f patologias -t 0.5

./Best_PRS.sh -c GE_2022 -I G40 -X true -S PLINK1.9 -f patologias 
./Best_PRS.sh -c GE_2022 -I G40 -X true -S PLINK2 -f patologias 
./Best_PRS.sh -c GE_2022 -I G40 -X true -S PLINK1.9 -f patologias -t 0.05
./Best_PRS.sh -c GE_2022 -I G40 -X true -S PLINK2 -f patologias -t 0.05
./Best_PRS.sh -c GE_2022 -I G40 -X true -S PLINK1.9 -f patologias -t 0.5
./Best_PRS.sh -c GE_2022 -I G40 -X true -S PLINK2 -f patologias -t 0.5


./Best_PRS.sh -c T2D_2022 -I E11 -X false -S PLINK1.9 -f patologias 
./Best_PRS.sh -c T2D_2022 -I E11 -X false -S PLINK2 -f patologias 
./Best_PRS.sh -c T2D_2022 -I E11 -X false -S PLINK1.9 -f patologias -t 0.05
./Best_PRS.sh -c T2D_2022 -I E11 -X false -S PLINK2 -f patologias -t 0.05
./Best_PRS.sh -c T2D_2022 -I E11 -X false -S PLINK1.9 -f patologias -t 0.5
./Best_PRS.sh -c T2D_2022 -I E11 -X false -S PLINK2 -f patologias -t 0.5

./Best_PRS.sh -c Altura_2022 -I false -X false -S PLINK1.9 -f patologias 
./Best_PRS.sh -c Altura_2022 -I false -X false -S PLINK2 -f patologias 
./Best_PRS.sh -c Altura_2022 -I false -X false -S PLINK1.9 -f patologias -t 0.05
./Best_PRS.sh -c Altura_2022 -I false -X false -S PLINK2 -f patologias -t 0.05
./Best_PRS.sh -c Altura_2022 -I false -X false -S PLINK1.9 -f patologias -t 0.5
./Best_PRS.sh -c Altura_2022 -I false -X false -S PLINK2 -f patologias -t 0.5







./Best_PRS.sh -c EPI_all_2022 -I G40 -X true -S PLINK1.9 -f Sin_Filtrar 
./Best_PRS.sh -c EPI_all_2022 -I G40 -X true -S PLINK2 -f Sin_Filtrar 
./Best_PRS.sh -c EPI_all_2022 -I G40 -X true -S PLINK1.9 -f Sin_Filtrar -t 0.05
./Best_PRS.sh -c EPI_all_2022 -I G40 -X true -S PLINK2 -f Sin_Filtrar -t 0.05
./Best_PRS.sh -c EPI_all_2022 -I G40 -X true -S PLINK1.9 -f Sin_Filtrar -t 0.5
./Best_PRS.sh -c EPI_all_2022 -I G40 -X true -S PLINK2 -f Sin_Filtrar -t 0.5

./Best_PRS.sh -c GE_2022 -I G40 -X true -S PLINK1.9 -f Sin_Filtrar 
./Best_PRS.sh -c GE_2022 -I G40 -X true -S PLINK2 -f Sin_Filtrar 
./Best_PRS.sh -c GE_2022 -I G40 -X true -S PLINK1.9 -f Sin_Filtrar -t 0.05
./Best_PRS.sh -c GE_2022 -I G40 -X true -S PLINK2 -f Sin_Filtrar -t 0.05
./Best_PRS.sh -c GE_2022 -I G40 -X true -S PLINK1.9 -f Sin_Filtrar -t 0.5
./Best_PRS.sh -c GE_2022 -I G40 -X true -S PLINK2 -f Sin_Filtrar -t 0.5


./Best_PRS.sh -c T2D_2022 -I E11 -X false -S PLINK1.9 -f Sin_Filtrar 
./Best_PRS.sh -c T2D_2022 -I E11 -X false -S PLINK2 -f Sin_Filtrar 
./Best_PRS.sh -c T2D_2022 -I E11 -X false -S PLINK1.9 -f Sin_Filtrar -t 0.05
./Best_PRS.sh -c T2D_2022 -I E11 -X false -S PLINK2 -f Sin_Filtrar -t 0.05
./Best_PRS.sh -c T2D_2022 -I E11 -X false -S PLINK1.9 -f Sin_Filtrar -t 0.5
./Best_PRS.sh -c T2D_2022 -I E11 -X false -S PLINK2 -f Sin_Filtrar -t 0.5

./Best_PRS.sh -c Altura_2022 -I false -X false -S PLINK1.9 -f Sin_Filtrar 
./Best_PRS.sh -c Altura_2022 -I false -X false -S PLINK2 -f Sin_Filtrar 
./Best_PRS.sh -c Altura_2022 -I false -X false -S PLINK1.9 -f Sin_Filtrar -t 0.05
./Best_PRS.sh -c Altura_2022 -I false -X false -S PLINK2 -f Sin_Filtrar -t 0.05
./Best_PRS.sh -c Altura_2022 -I false -X false -S PLINK1.9 -f Sin_Filtrar -t 0.5
./Best_PRS.sh -c Altura_2022 -I false -X false -S PLINK2 -f Sin_Filtrar -t 0.5








echo "DONE"