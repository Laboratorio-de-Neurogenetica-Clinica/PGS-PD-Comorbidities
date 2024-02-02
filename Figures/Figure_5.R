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

############################### Figure 4 , PRS association with sex

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
#OR for Cmrbs Vs SEX
PRS_OR_sex <- function(PRS,Disease,subset="PD",age="0",group="SEX"){
    DATA <- load_data(PRS)
    subset_carga <- as.character(subset)
    if (subset == FALSE) {
        data_test <- DATA[DATA$PD == 0 ]
    } else {
        data_test <- DATA[DATA[[subset_carga]] == 1 ]
    }
    
    ##Filter by Age
    data <- data_test[data_test$AGE >= age]
    
    
    if (subset == "PD") {
        
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
        data_PD_AAO_all <- data_PD_AAO_all[!is.na(data_PD_AAO_all$AGE_plot_Positive)]
        data_PD_AAO <- subset(data_PD_AAO_all, data_PD_AAO_all$AAO_PD >= 50 )
        
        data <- data_PD_AAO
    }
    
    
    #groups
    
    data_male <- data[data$SEX == "Male"] 
    data_female <- data[data$SEX == "Female"] 
    
    data_SEX <- list(data_male,
                     data_female)
    labels_all <- c("Male","Female")
    names(data_SEX) <- labels_all
    VALUES_TRESH <- data.frame(NULL)
    for (i in 1:length(data_SEX)) {
        i <- as.numeric(i)
        
        data_SEX[[i]]$SCORE_20_80 <- ifelse(data_SEX[[i]]$zSCORE <= quantile(data_SEX[[i]]$zSCORE,0.8), 0, 1) 
        values_tresh <- data.frame(NULL)
        values_tresh <- cbind(Percentil=c(20))
        
        formula <- as.formula(paste0(Disease, " ~ SCORE_20_80 + AGE +PC1+PC2+PC3+PC4 "))
        glm_tresh <- glm(formula,family = "binomial",data=data_SEX[[i]] )
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
        values_tresh$SEX <- names(data_SEX)[i]
        
        #numero casos controles
        var <- "SCORE_20_80"
        case_control <- data_SEX[[i]] %>% 
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
    
#   data$SCORE_20_80 <- ifelse(data$zSCORE <= quantile(data$zSCORE,0.8), 0, 1) 
#    values_tresh <- data.frame(NULL)
#    values_tresh <- cbind(Percentil=c(20))
#    
#    formula <- as.formula(paste0(Disease, " ~ SCORE_20_80 * SEX  + AGE +PC1+PC2+PC3+PC4 "))
#    glm_tresh <- glm(formula,family = "binomial",data=data )
#    BETA <-  summary(glm_tresh)$coeff[9,1]
#    SE <- summary(glm_tresh)$coeff[9,2]
#    p.value <- summary(glm_tresh)$coeff[9,4]
#    values_tresh <- cbind(values_tresh,BETA)
#    values_tresh <- cbind(values_tresh,SE)
#    values_tresh <- as.data.frame(cbind(values_tresh,p.value))
#    
#    
#    values_tresh$BETA_low <- values_tresh$BETA - 1.96 * values_tresh$SE  
#    values_tresh$BETA_high <- values_tresh$BETA + 1.96 * values_tresh$SE  
#    values_tresh$OR <- signif(exp(values_tresh$BETA),3)  
#    values_tresh$OR_low <-  signif(exp(values_tresh$BETA_low),3)
#    values_tresh$OR_high <-  signif(exp(values_tresh$BETA_high),3)
#    
#    values_tresh$OR_signif <- signif(values_tresh$OR,2)
#    values_tresh$CMRB <- Disease
#    values_tresh$Population <- ifelse(subset==FALSE,"GP",subset)
#    values_tresh$SEX <- "All"
    
    #numero casos controles
    #var <- "SCORE_20_80"
    #case_control <- data %>% 
    #    group_by(.data[[var]]) %>% 
    #    summarise(Total=n(),Cases=sum( .data[[Disease]] ),Controls= (Total - Cases))
    #
#    Table <- data.frame(NULL)
#    Table <- rbind(Table,c(ifelse(subset==FALSE,"GP",subset),
#                           Disease,
#                           paste0(case_control[1,3],"/",case_control[1,4]),
#                           paste0(case_control[2,3],"/",case_control[2,4])) )
#    
#    
#    names(Table) <- c("Population","CMRB","Low_PRS","High_PRS")
#    Table
    #values_tresh <- cbind(values_tresh,Table[,c(3,4)])
#    VALUES_TRESH <- rbind(VALUES_TRESH,values_tresh)
    
    
}
#plot Previous function
plot_PRS_OR_sex <- function(data ) {
    data$SEX <- factor(data$SEX,levels=c("Male","Female"))
    
    data_to_plot <- data 
    
    
    
    ###compare ODDS ratio Age
    plot1 <- ggplot(data_to_plot,aes(y=as.factor(CMRB),x=as.numeric(OR_signif ), shape= SEX , color = Population)) + 
        geom_pointrange(aes(xmin=as.numeric(OR_low ),xmax=as.numeric(OR_high )),
                        position= position_dodge(width = 0.6) ,
                        size=1.3)+
        geom_vline(xintercept = 1.0) + 
        theme_bw() + 
        xlab("OR")  +
        ylab("") +
        labs(colour=paste("Population"),
             shape= "Sex") +
        scale_y_discrete(limits=c("EPI","MH","MDD","T2D"))  +
        theme_pubclean()
    plot1
    plot2 <- plot1 + theme(legend.title = element_text(color = "black", size = 30),
                           legend.text = element_text(color = "black", size = 25),
                           legend.background = element_rect(fill = "white", colour = "gray30"),
                           axis.text.y = element_text(size=30),
                           axis.text.x = element_text(size=30) ,
                           axis.title = element_text(size=30)) +
        scale_x_continuous(limits=c(0,4)) +
        scale_color_manual(values = c(ifelse(subset == "PD" , "#1B9E77", "#1B9E77"),
                                      ifelse(subset == "All" , "#D95F02","#D95F02" )
        )
        )
    
    
    plot2
    
}

#Make plots
{
age=50
OR_Pvalue_sex <- data.frame(NULL)
for (i in list(c("EPI","EPI"),c("T2D","T2D"),c("MH","MH"),c("MDD","MDD"))) { #c("EPI","Insomnia","MDD","T2D")
    OR_Pvalue_sex <- rbind(OR_Pvalue_sex,PRS_OR_sex(i[[1]],i[[2]],"PD",age,group="SEX"))
    OR_Pvalue_sex <- rbind(OR_Pvalue_sex,PRS_OR_sex(i[[1]],i[[2]],FALSE,age,group="SEX"))
}

} # Run OR_Pvalue_sex() in all CMRBs 
p1 <- plot_PRS_OR_sex(OR_Pvalue_sex) %>% label_plot(., "A",25)


#landscape 30*20
ggsave(file="All_figures/Figure5.pdf",p1,width = 11,height = 11)

#Extract info to make table sup 4
OR_Pvalue_sex$Ci_95 <- paste0(OR_Pvalue_sex$OR_low," - ",OR_Pvalue_sex$OR_high )
TABLE <- OR_Pvalue_sex[,c("Population","CMRB","SEX","OR","Ci_95","p.value","High_PRS","Low_PRS")]
TABLE$p.value <- signif(TABLE$p.value, 3)

Table_SUP4 <- TABLE
fwrite( Table_SUP4 , "Tables/Table_SUP4.tsv",sep="\t" )
