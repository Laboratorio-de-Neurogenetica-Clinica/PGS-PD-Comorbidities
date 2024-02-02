#!/home/carlos/miniconda3/envs/r4.1/bin/Rscript --vanilla
#
#
#
#
library(data.table)
library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(ggpubr)
library(rstatix)
library(stringr)
library(magrittr)
library(survival)
library(finalfit)
library(survminer)


############################### Figure sup 4 ,PRS in age at onset

##define functions

#Load data
load_data <- function (PRS) {
  data <- fread(paste0(PRS,"_PRS_pheno_covar.tsv"))
  names(data)[3] <- "zSCORE"
  data
}
#Label Plot
label_plot <- function (plot,label,size) {
    plot + labs(tags=label) + theme(plot.tag = element_text(size=size, face= "bold") )
}
#compare AAO of PD between CMRBs
violin_plot_AAO_CMRB <- function(PRS,Disease,subset,age=0,onset=0){
    DATA <- load_data(PRS)
    if (subset == "Disease") {
        data_test <- DATA[DATA[[Disease]] == 1 ]
    } else if (subset== "PD") {
        data_test <- DATA[DATA$PD == 1 ]
    }
    ##Filter by Age
    data <- data_test[data_test$AGE >= age]
    
    
    ## cargar AAO PD
    #AAO <-  fread ("../Data/CLINICAL_COVAR_QC_participant.tsv",sep="\t")
    #which(str_detect(colnames(AAO), regex("Date G40|Date G20|Date E11|Date F32|Date G470|Date G43|Date G47")))
    AAO <-  load_data(PRS)[,c(1,25,26,27)]  
    colnames(AAO) <- c("FID","Birth","Date_PD",paste0("Date_",PRS))
    
    #AAO PD
    AAO$Year_PD <- as.numeric(format(AAO$Date_PD,'%Y'))
    AAO$AAO_PD <- AAO$Year_PD - AAO$Birth
    #AAO CMRB
    AAO$Date_CMRB <- ifelse(str_detect(AAO[[paste0("Date_",Disease)]],"Code"),
                            "",
                            ifelse(is.na(AAO[[paste0("Date_",Disease)]]),
                                   "",
                                   AAO[[paste0("Date_",Disease)]]))
    AAO$Date_CMRB <- as.Date(AAO$Date_CMRB)
    AAO$Year_CMRB <- as.numeric(format(AAO$Date_CMRB,'%Y'))
    AAO$AAO_CMRB <- AAO$Year_CMRB - AAO$Birth
    
    ## merge data
    data_PD_AAO_all <- merge(data,AAO, by="FID")
    if (subset == "PD") {
        data_PD_AAO_all <- data_PD_AAO_all[!is.na(data_PD_AAO_all$AAO_PD)]
        data_onset <- data_PD_AAO_all[data_PD_AAO_all$AAO_PD >= as.numeric(onset)]
    } else if (subset == "Disease") {
        data_PD_AAO_all <- data_PD_AAO_all[!is.na(data_PD_AAO_all$AAO_CMRB)]
        data_onset <- data_PD_AAO_all[data_PD_AAO_all$AAO_CMRB >= as.numeric(onset)]
        
    }
    
    if (subset == "Disease") {
        Age_mean <- data_onset %>% 
            group_by(as.factor(PD)) %>% 
            summarise(count = n(),MEAN=mean(AGE),SD=sd(AGE), MAX=max(AGE))
        
        stat.test <- data_onset %>% t_test(AAO_CMRB ~ PD)
        stat.test <- stat.test %>% add_xy_position(x = "PD")
        stat.test$p.adj <- stat.test$p
        
    } else if (subset == "PD") {
        Age_mean <- data_onset %>% 
            group_by(as.factor( data_onset[[Disease]] ) ) %>% 
            summarise(count = n(),MEAN=mean(AGE),SD=sd(AGE), MAX=max(AGE))
        
        formula <- as.formula(paste0("AAO_PD ~ " , Disease))
        stat.test <- data_onset %>% t_test(formula)
        stat.test <- stat.test %>% add_xy_position(x = Disease)
        stat.test$p.adj <- stat.test$p
    }
    
    
    if (subset == "Disease") {
        ggplot(data_onset, aes(x = as.factor(PD), y = AAO_CMRB)) +
            geom_violin(aes(fill = as.factor(PD)),
                        show.legend = FALSE) +
            geom_boxplot(width=0.15, 
                         fill="white",
                         outlier.shape = NA) +
            stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.005, 
                               step.increase = 0.02,
                               size = 6) +
            #annotate(geom = "text", 
            #         label =  paste0(round(Age_mean$MEAN[1],0),"±", round(Age_mean$SD[1],0) ) ,
            #         x= 0.5, y= 50) +
            #annotate(geom = "text", label =  paste0(round(Age_mean$MEAN[2],0),"±", round(Age_mean$SD[2],0) ),
            #         x= 2.5, y= 50) +
            
            theme_minimal() +
            theme_bw() +
            ylim(onset,90)+
            scale_x_discrete(labels=c("Controls","PD")) +
            scale_fill_manual(values=c("orange","#BEAED4"),
            ) +
            ylab(paste0("Age at Onset of ",Disease ))  +
            xlab(paste0("Individuals with ", Disease))  + 
            
            theme(axis.text.x = element_text(size=20),
                  axis.title = element_text(size=20),
                  axis.text.y= element_text(size=15))
        
        
    }else if (subset == "PD") {
        ggplot(data_onset, aes(x = as.factor(data_onset[[Disease]]), y = AAO_PD)) +
            geom_violin(aes(fill = as.factor(data_onset[[Disease]])),
                        show.legend = FALSE) +
            geom_boxplot(width=0.15, 
                         fill="white",
                         outlier.shape = NA) +
            stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.005, 
                               step.increase = 0.02,
                               size = 6) +
            #annotate(geom = "text", 
            #         label =  paste0(round(Age_mean$MEAN[1],0),"±", round(Age_mean$SD[1],0) ) ,
            #         x= 0.5, y= 50) +
            #annotate(geom = "text", label =  paste0(round(Age_mean$MEAN[2],0),"±", round(Age_mean$SD[2],0) ),
            #         x= 2.5, y= 50) +
            
            theme_minimal() +
            theme_bw() + 
            ylim(onset,90)+
            scale_x_discrete(labels=c("Control",Disease)) +
            scale_fill_manual(values=c("#7FC97F","#BEAED4")) +
            ylab(paste0("Age at Onset of PD "))  +
            xlab("PD subset") + 
            theme(axis.text.x = element_text(size=20),
                  axis.title = element_text(size=20),
                  axis.text.y= element_text(size=15))
    }
    
}

#make plots
p1 <- violin_plot_AAO_CMRB("T2D","T2D",subset = "PD", age=0, onset=50) %>% label_plot(.,"A",25)
p2 <- violin_plot_AAO_CMRB("MDD","MDD",subset = "PD", age=0, onset=50) %>% label_plot(.,"B",25)
p3 <- violin_plot_AAO_CMRB("MH","MH",subset = "PD", age=0, onset=50) %>% label_plot(.,"C",25)
p4 <- violin_plot_AAO_CMRB("EPI","EPI",subset = "PD", age=50, onset=50) %>% label_plot(.,"D",25)

##grid FIG SUP 5 
lay <- rbind(c(1,2),
             c(3,4))
grid_landscape <- arrangeGrob(grobs = list(p1,p2,p3,p4), layout_matrix = lay)
#landscape 15*15
ggsave(file="All_figures/FigureSup5.pdf",grid_landscape,width = 15,height = 15)

