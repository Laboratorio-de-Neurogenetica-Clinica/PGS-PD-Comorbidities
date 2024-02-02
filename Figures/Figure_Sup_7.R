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


############################### Figure sup 7 ,PRS and age at onset of CMRB disease

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
#Plot Prevalence PD Vs onset of CMRB in PRS deciles
Decile_AGE <- function(PRS,Disease,subset="PD",age=0,onset=0){
    DATA <- load_data(PRS)
    subset <- as.character(subset)
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
        data_PD_AAO_all$AGE_plot <-  data_PD_AAO_all$AAO_PD
    }else if (subset == "Disease"){
        data_PD_AAO_all$AGE_plot <-  data_PD_AAO_all$AAO_CMRB
    }
    data_PD_AAO <- subset(data_PD_AAO_all, data_PD_AAO_all$AGE_plot >= onset )
    
    
    #Deciles
    data_PD_AAO$dec1 <- ifelse(data_PD_AAO$zSCORE <= quantile(data_PD_AAO$zSCORE,0.20), 1, 0)
    data_PD_AAO$dec2_9 <- ifelse(data_PD_AAO$zSCORE > quantile(data_PD_AAO$zSCORE,0.33) & data_PD_AAO$zSCORE <= quantile(data_PD_AAO$zSCORE,0.66), 1, 0)
    data_PD_AAO$dec10 <- ifelse(data_PD_AAO$zSCORE > quantile(data_PD_AAO$zSCORE,0.80), 1, 0)
    
    #agrupate in one varaible
    data_PD_AAO$dec[data_PD_AAO$dec1 == 1] <- "low"
    #data_PD_AAO$dec[data_PD_AAO$dec2_9 == 1] <- "medium"
    data_PD_AAO$dec[data_PD_AAO$dec10 == 1] <- "high"
    
    #phantom_data <- data_PD_AAO[1,]
    #phantom_data$MH <- 0
    #phantom_data$AAO_PD <- 30
    #phantom_data_1 <- phantom_data
    #phantom_data_2 <- phantom_data
    #phantom_data_3 <- phantom_data
    #phantom_data_1$dec <- "high"
    #phantom_data_2$dec <- "medium"
    #phantom_data_3$dec <- "low"
    
    #data_PD_AAO_with0 <- rbind(data_PD_AAO,phantom_data_1,phantom_data_3) 
    data_PD_AAO_with0 <- data_PD_AAO
    data_PD_AAO_with0 <- data_PD_AAO_with0[!is.na(data_PD_AAO_with0$dec)]
    total_low <- sum(data_PD_AAO_with0$dec == "low")
    total_high <-sum(data_PD_AAO_with0$dec == "high")
    
    if (subset== "PD"){
        data_plot <- data_PD_AAO_with0 %>% 
            select(FID,zSCORE,Disease,AGE,PD,AGE_plot,dec) %>% 
            group_by(dec,AGE_plot) %>% 
            summarise(cases=sum(!!as.symbol(  Disease) )  )
    } else if (subset== "Disease" ) {
        data_plot <- data_PD_AAO_with0 %>% 
            select(FID,zSCORE,PD,AGE,AGE_plot,dec) %>% 
            group_by(dec,AGE_plot) %>% 
            summarise(cases=sum(PD )  )
    }
    data_plot$N_total <- ifelse(data_plot$dec == "low", total_low, total_high ) 
    
    data_plot <- data_plot %>%
        mutate (prevalence=(cases/N_total)*100,
                accumulate = cumsum(prevalence))
    
    
    data_plot
    
    
    #plot
    plot1 <- ggplot(data_plot,aes(AGE_plot, accumulate, color=dec)) +
        geom_point() +
        geom_line() +
        theme_bw() + 
        labs(colour=paste("PGS Groups"))  +
        scale_color_manual(values = c("high" = "#D95F02", "low" = "#1B9E77"),
                           limits = c("high", "low"),
                           breaks = c("high", "low"),
                           labels = c("High", "Low")) +
        theme_pubclean()
    plot1
    plot2 <- plot1 + theme(legend.title = element_text(color = "black", size = 20),
                           legend.text = element_text(color = "black", size = 15),
                           legend.background = element_rect(fill = "white", colour = "gray30"),
                           axis.text.y = element_text(size=20),
                           axis.text.x = element_text(size=20) ,
                           axis.title = element_text(size=25))
    plot2
    if (subset == "PD") {
        plot3 <- plot2 + 
            xlab("Age At Onset of PD (Years)")  +
            ylab(paste0("% Prevalence of ", Disease))     
    } else if (subset == "Disease") {
        plot3 <- plot2 + 
            xlab(paste0("Age At Onset of ",Disease ,"(Years)") ) +
            ylab(paste0("% Prevalence of PD"))   
    }
    plot3
}
#OR for PD before CMRB onset
PRS_OR_event_time <- function(PRS,Disease,subset="PD",age=0,onset=0){
    DATA <- load_data(PRS)
    subset_carga <- as.character(subset)
    if (subset == FALSE) {
        data_test <- DATA
    } else {
        data_test <- DATA[DATA[[subset_carga]] == 1 ]
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
    data_PD_AAO_all <- data_PD_AAO_all[!is.na(data_PD_AAO_all$AAO_CMRB)]
    data_PD_AAO_all$onsetPD <-  as.numeric(ifelse(data_PD_AAO_all$AAO_PD >= data_PD_AAO_all$AAO_CMRB,
                                                  "1",
                                                  "0"))
    data_PD_AAO <- subset(data_PD_AAO_all, data_PD_AAO_all$AAO_PD >= 50)
    
    
    data_onset <- data_PD_AAO
    
    
    data_onset$SCORE_20_80 <- ifelse(data_onset$zSCORE <= quantile(data_onset$zSCORE,0.8), 0, 1) 
    values_tresh <- data.frame(NULL)
    values_tresh <- cbind(Percentil=c(20))
    
    formula <- as.formula(paste0(" onsetPD ~ SCORE_20_80 + AGE + SEX +PC1+PC2+PC3+PC4 "))
    glm_tresh <- glm(formula,family = "binomial",data=data_onset )
    BETA <-  summary(glm_tresh)$coeff[2,1]
    SE <- summary(glm_tresh)$coeff[2,2]
    p.value <- summary(glm_tresh)$coeff[2,4]
    values_tresh <- cbind(values_tresh,BETA)
    values_tresh <- cbind(values_tresh,SE)
    values_tresh <- as.data.frame(cbind(values_tresh,p.value))
    
    
    values_tresh$BETA_low <- values_tresh$BETA - 1.96 * values_tresh$SE  
    values_tresh$BETA_high <- values_tresh$BETA + 1.96 * values_tresh$SE  
    values_tresh$OR <- signif(exp(values_tresh$BETA),3)  
    values_tresh$OR_low <-  signif(exp(values_tresh$BETA_low),3)
    values_tresh$OR_high <-  signif(exp(values_tresh$BETA_high),3)
    
    values_tresh$OR_signif <- signif(values_tresh$OR,2)
    values_tresh$CMRB <- Disease
    values_tresh$Population <- ifelse(subset==FALSE,"GP",subset)
    values_tresh
}
#plot prevous OR
plot_PRS_OR_event_time <- function(data) {
    ###compare ODDS ratio Age
    plot1 <- ggplot(data,aes(y=as.factor(CMRB),x=as.numeric(OR_signif ))) + 
        geom_pointrange(aes(xmin=as.numeric(OR_low ),xmax=as.numeric(OR_high )),
                        position=position_dodge2(width=0.5,reverse = TRUE),
                        size=1.3)+
        geom_vline(xintercept = 1.0) + 
        theme_bw() + 
        xlab("OR for future PD diagnosis")  +
        ylab("")  +
        scale_y_discrete(limits=c("EPI","MH","T2D","MDD") )+
        theme_pubclean()
    plot1
    plot2 <- plot1 + theme(
        axis.text.y = element_text(size=20),
        axis.text.x = element_text(size=20) ,
        axis.title = element_text(size=25))
    plot2
}

#make plots
q1 <- violin_plot_AAO_CMRB("T2D","T2D",subset = "Disease", age=50, onset=0) %>% label_plot(.,"A",25)
q2 <- violin_plot_AAO_CMRB("MDD","MDD",subset = "Disease", age=50, onset=0) %>% label_plot(.,"B",25)
q3 <- violin_plot_AAO_CMRB("MH","MH",subset = "Disease", age=50, onset=0) %>% label_plot(.,"C",25)
q4 <- violin_plot_AAO_CMRB("EPI","EPI",subset = "Disease", age=50, onset=0) %>% label_plot(.,"D",25)
q5 <- Decile_AGE("T2D","T2D",subset="Disease",age=50,onset=0) %>% label_plot(., "E",25)
q6 <- Decile_AGE("MDD","MDD",subset="Disease",age=50,onset=0) %>% label_plot(., "F",25)
q7 <- Decile_AGE("MH","MH",subset="Disease",age=50,onset=0) %>% label_plot(., "G",25)
q8 <- Decile_AGE("EPI","EPI",subset="Disease",age=50,onset=0)  %>% label_plot(., "H",25)
{
OR_Pvalue_event_time <- data.frame(NULL)
for (i in list(c("EPI","EPI"),c("T2D","T2D"),c("MH","MH"),c("MDD","MDD"))) { #c("EPI","Insomnia","MDD","T2D")
    OR_Pvalue_event_time <- rbind(OR_Pvalue_event_time,PRS_OR_event_time(i[[1]],i[[2]],"PD",age=0,onset = 50))
}
OR_Pvalue_event_time
} # Run PRS_OR_event_time() in all cmrbs
p1 <- plot_PRS_OR_event_time(OR_Pvalue_event_time) %>% label_plot(., "I",25)

#####grid
lay <- rbind(c(1,2,9,9),
             c(3,4,9,9),
             c(5,6,7,8))
grid_landscape <- arrangeGrob(grobs = list(q1,q2,q3,q4,q5,q6,q7,q8,p1), layout_matrix = lay)
#landscape 20*15
ggsave(file="All_figures/FigureSup7.pdf",grid_landscape,width = 20,height = 15)        





