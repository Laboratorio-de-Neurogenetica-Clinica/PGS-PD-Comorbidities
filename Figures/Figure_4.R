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

############################### Figure 4 , PRS asocaition with age at onset

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
## Plot Prevalence CMRB Vs Time in PRS deciles
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
#OR for Cmrb Vs AGE at onset Categoric
PRS_OR_onset <- function(PRS,Disease,subset="PD",age="0",onset="0",groups=c(0,50,70,Inf)){
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
    AAO <-  load_data(PRS)[,c(1,25,26)]  
    colnames(AAO) <- c("FID","Birth","Year_PD")
    AAO$Year_PD <- as.numeric(format(AAO$Year_PD,'%Y'))
    AAO$AAO_PD <- AAO$Year_PD - AAO$Birth
    
    
    
    ## merge data
    data_PD_AAO_all <- merge(data,AAO, by="FID")
    data_PD_AAO_all$AGE_plot <-  data_PD_AAO_all$AAO_PD - data_PD_AAO_all$AGE
    min_age <- abs(min( data_PD_AAO_all$AGE_plot,na.rm=TRUE))
    data_PD_AAO_all$AGE_plot_Positive <-  data_PD_AAO_all$AGE_plot + min_age
    data_PD_AAO_all$AGE_plot_Positive[!is.na(data_PD_AAO_all$AGE_plot_Positive)]
    
    #filter by onset
    data_PD_AAO <- subset(data_PD_AAO_all, data_PD_AAO_all$AAO_PD >= onset)
    
    #Deciles
    data_PD_AAO$dec1 <- ifelse(data_PD_AAO$zSCORE <= quantile(data_PD_AAO$zSCORE,0.33), 1, 0)
    data_PD_AAO$dec2_9 <- ifelse(data_PD_AAO$zSCORE > quantile(data_PD_AAO$zSCORE,0.33) & data_PD_AAO$zSCORE <= quantile(data_PD_AAO$zSCORE,0.66), 1, 0)
    data_PD_AAO$dec10 <- ifelse(data_PD_AAO$zSCORE > quantile(data_PD_AAO$zSCORE,0.66), 1, 0)
    
    #agrupate in one varaible
    data_PD_AAO$dec <- "min"
    data_PD_AAO$dec[data_PD_AAO$dec2_9 == 1] <- "medium"
    data_PD_AAO$dec[data_PD_AAO$dec10 == 1] <- "low"
    
    
    #groups
    labels_1 <- ifelse(groups[1] == 0, paste0("<",groups[2]),paste0(groups[1],"-",groups[2])) 
    labels_2 <- paste0(groups[2],"-",groups[3])
    labels_3 <- ifelse(groups[4] == Inf, paste0(">",groups[3]),paste0(groups[3],"-",groups[4])) 
    labels_all <- c(labels_1,labels_2,labels_3)
    
    
    
    data_PD_AAO$AAO_PD_CAT <- cut(data_PD_AAO$AAO_PD, breaks =groups,
                                  labels = labels_all,
                                  right = FALSE)
    data_PD_AAO <- data_PD_AAO[!is.na(data_PD_AAO$AAO_PD_CAT)]
    data_menor <- data_PD_AAO[data_PD_AAO$AAO_PD_CAT == labels_1] 
    data_media <- data_PD_AAO[data_PD_AAO$AAO_PD_CAT == labels_2] 
    data_mayor <- data_PD_AAO[data_PD_AAO$AAO_PD_CAT == labels_3] 
    data_EDADES <- list(data_menor,
                        data_media,
                        data_mayor)
    names(data_EDADES) <- labels_all
    VALUES_TRESH <- data.frame(NULL)
    for (i in 1:length(data_EDADES)) {
        i <- as.numeric(i)
        
        data_EDADES[[i]]$SCORE_20_80 <- ifelse(data_EDADES[[i]]$zSCORE <= quantile(data_EDADES[[i]]$zSCORE,0.8), 0, 1) 
        values_tresh <- data.frame(NULL)
        values_tresh <- cbind(Percentil=c(20))
        
        formula <- as.formula(paste0(Disease, " ~ SCORE_20_80 + SEX+ AGE +PC1+PC2+PC3+PC4 "))
        glm_tresh <- glm(formula,family = "binomial",data=data_EDADES[[i]] )
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
        values_tresh$AAO_PD_CAT <- names(data_EDADES)[i]
        
        
        #numero casos controles
        var <- "SCORE_20_80"
        case_control <- data_EDADES[[i]] %>% 
            group_by(.data[[var]]) %>% 
            summarise(Total=n(),Cases=sum( .data[[Disease]] ),Controls= (Total - Cases))
        
        Table <- data.frame(NULL)
        Table <- rbind(Table,c(ifelse(subset==FALSE,"GP",subset),
                               Disease,
                               paste0(case_control[1,3],"/",case_control[1,4]),
                               paste0(case_control[2,3],"/",case_control[2,4])) )
        
        
        names(Table) <- c("Population","CMRB","Low_PRS","High_PRS")
        Table
        row <- cbind(values_tresh,Table[,c(3,4)])
        
        VALUES_TRESH <- rbind(VALUES_TRESH,row)
        
        
    }
    VALUES_TRESH
    data_PD_AAO$SCORE_20_80 <- ifelse(data_PD_AAO$zSCORE <= quantile(data_PD_AAO$zSCORE,0.8), 0, 1) 
    values_tresh <- data.frame(NULL)
    values_tresh <- cbind(Percentil=c(20))
    
    formula <- as.formula(paste0(Disease, " ~ SCORE_20_80 * AAO_PD_CAT + SEX+ AGE +PC1+PC2+PC3+PC4 "))
    glm_tresh <- glm(formula,family = "binomial",data=data_PD_AAO )
    BETA <-  summary(glm_tresh)$coeff[10,1]
    SE <- summary(glm_tresh)$coeff[10,2]
    p.value <- summary(glm_tresh)$coeff[10,4]
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
    values_tresh$Population <- "PD"
    values_tresh$AAO_PD_CAT <- "All"
    
    #numero casos controles
    var <- "SCORE_20_80"
    case_control <- data_PD_AAO %>% 
        group_by(.data[[var]]) %>% 
        summarise(Total=n(),Cases=sum( .data[[Disease]] ),Controls= (Total - Cases))
    
    Table <- data.frame(NULL)
    Table <- rbind(Table,c(ifelse(subset==FALSE,"GP",subset),
                           Disease,
                           paste0(case_control[1,3],"/",case_control[1,4]),
                           paste0(case_control[2,3],"/",case_control[2,4])) )
    
    
    names(Table) <- c("Population","CMRB","Low_PRS","High_PRS")
    Table
    values_tresh <- cbind(values_tresh,Table[,c(3,4)])
    VALUES_TRESH <- rbind(VALUES_TRESH,values_tresh)
}
#plot Previous function
plot_PRS_OR_onset <- function(data,groups=c("50-70",">70")) {
    data$AAO_PD_CAT <- factor(data$AAO_PD_CAT,levels=c("<50",">70","50-70","All"))
    data_to_plot <- data[data$AAO_PD_CAT %in% groups,]
    
    Pvalues <- data[data$AAO_PD_CAT == "All",][,c("p.value","CMRB","AAO_PD_CAT")]
    Pvalues$OR <- data[data$AAO_PD_CAT == ">70",]$OR_signif
    #data <- data[data$AAO_PD_CAT != "<50"]
    
    ###compare ODDS ratio Age
    plot1 <- ggplot(data_to_plot,aes(y=as.factor(CMRB),x=as.numeric(OR_signif ), color=AAO_PD_CAT)) + 
        geom_pointrange(aes(xmin=as.numeric(OR_low ),xmax=as.numeric(OR_high )),
                        position=position_dodge2(width=0.5,reverse = TRUE),
                        size=1.3)+
        geom_vline(xintercept = 1.0) + 
        theme_bw() + 
        xlab("OR")  +
        ylab("") +
        labs(colour=paste("AAO PD")) +
        scale_color_manual(values = c(">70" = "#D95F02", "50-70" = "#7570B3"),
                           limits = c(">70", "50-70"),
                           breaks = c(">70", "50-70"),
                           labels = c(">70", "50-70")) +    
        scale_y_discrete(limits=c("EPI","MH","MDD","T2D")) +
        theme_pubclean()
    plot2 <- plot1 + theme(legend.title = element_text(color = "black", size = 30),
                           legend.text = element_text(color = "black", size = 25),
                           legend.background = element_rect(fill = "white", colour = "gray30"),
                           axis.text.y = element_text(size=30),
                           axis.text.x = element_text(size=30) ,
                           axis.title = element_text(size=30)) +
        scale_x_continuous(limits=c(0,4.5))
    
    plot2
    
}

#Make plots
p1 <- Decile_AGE("T2D","T2D",subset="PD",onset=50) %>% label_plot(., "B",25)
p2 <- Decile_AGE("MDD","MDD",subset="PD",onset=50) %>% label_plot(., "C",25)
p3 <- Decile_AGE("MH","MH",subset="PD",onset=50) %>% label_plot(., "D",25)
p4 <- Decile_AGE("EPI","EPI",subset="PD",onset=50)  %>% label_plot(., "E",25)
{
    age=50
    OR_Pvalue <- data.frame(NULL)
    for (i in list(c("EPI","EPI"),c("T2D","T2D"),c("MDD","MDD"),c("MH","MH"))) { #c("EPI","Insomnia","MDD","T2D")
        OR_Pvalue <- rbind(OR_Pvalue,PRS_OR_onset(i[[1]],i[[2]],"PD",age, onset= 0))
    }
    
    
} #Run PRS_OR_onset in all cmrbs    
p5 <- plot_PRS_OR_onset(OR_Pvalue,groups = c("50-70",">70") )  %>% label_plot(., "A",20)

##Arrenge in grid
lay <- rbind(c(6,6,1,2),
             c(6,6,3,4))
grid_landscape <- arrangeGrob(grobs = list(p1,p2,p3,p4,p5), layout_matrix = lay)
#landscape 12*20
ggsave(file="All_figures/Figure4.pdf",grid_landscape,width = 20,height = 12)
    
#Extract info to make table sup 3
OR_Pvalue$Ci_95 <- paste0(OR_Pvalue$OR_low," - ",OR_Pvalue$OR_high )
TABLE <- OR_Pvalue[,c("Population","CMRB","AAO_PD_CAT","OR","Ci_95","p.value","High_PRS","Low_PRS")]
TABLE$p.value <- signif(TABLE$p.value, 3)
    
Table_SUP3 <- TABLE[TABLE$Population == "PD",]
fwrite( Table_SUP3 , "Tables/TableSup3.tsv",sep="\t" )
    
    
   