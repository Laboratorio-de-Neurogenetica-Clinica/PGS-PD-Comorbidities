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
library(GGally)
library(cowplot)
library(magrittr)

############################### Figure 3 , OR for PD and GP / regression between PRSs

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
###### Plot PRS Violin Plots
PRS_OR <- function(PRS,Disease,subset=FALSE,age=50){
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
  data$SCORE_20_80 <- ifelse(data$zSCORE <= quantile(data$zSCORE,0.8), 0, 1) 
  data$SCORE_10_90 <- ifelse(data$zSCORE <= quantile(data$zSCORE,0.9), 0, 1) 
  data$SCORE_5_95 <- ifelse(data$zSCORE <= quantile(data$zSCORE,0.95), 0, 1)
  data$SCORE_1_99 <- ifelse(data$zSCORE <= quantile(data$zSCORE,0.99), 0, 1)
  
  values_tresh <- data.frame(NULL)
  values_tresh <- cbind(Percentil=c(20,10,5,1))
  tresh <- c()
  tresh <-c(tresh, as.numeric(min(data[data$SCORE_20_80 == 1 ]$zSCORE )))
  tresh <-c(tresh, as.numeric(min(data[data$SCORE_10_90 == 1 ]$zSCORE )))
  tresh <-c(tresh, as.numeric(min(data[data$SCORE_5_95 == 1 ]$zSCORE )))
  tresh <-c(tresh, as.numeric(min(data[data$SCORE_1_99 == 1 ]$zSCORE )))
  values_tresh <- cbind(values_tresh,tresh)
  COEFFS <- data.frame(NULL)
  for (i in c(20,10,5,1)) {
    # Data pvalue
    formula <- as.formula(paste0(Disease, " ~ ", paste0("SCORE_",i,"_",100-i), " + SEX+ AGE +PC1+PC2+PC3+PC4 "))
    glm_tresh <- glm(formula,family = "binomial",data=data )
    BETA <-  summary(glm_tresh)$coeff[2,1]
    SE <- summary(glm_tresh)$coeff[2,2]
    p.value <- summary(glm_tresh)$coeff[2,4]
    valores <- data.frame(BETA,SE,p.value)
    
    #numero casos controles
    var <- paste0("SCORE_",i,"_",100-i)
    case_control <- data %>% 
    group_by(.data[[var]]) %>% 
        summarise(Total=n(),Cases=sum( .data[[Disease]] ),Controls= (Total - Cases))

    Table <- data.frame(NULL)
    Table <- rbind(Table,c(ifelse(subset==FALSE,"GP",subset),
                               Disease,
                               paste0(case_control[1,3],"/",case_control[1,4]),
                               paste0(case_control[2,3],"/",case_control[2,4])) )
        
    
    names(Table) <- c("Population","CMRB","Low_PRS","High_PRS")
    Table
    Row <- cbind(valores,Table)
    COEFFS <- rbind(COEFFS,Row)
  }
  values_tresh <- cbind(values_tresh,COEFFS)
  names(values_tresh) <- c("Percentil","Tresh","BETA","SE","p.value","Population","CMRB","Low_PRS","High_PRS")
  values_tresh$BETA_low <- values_tresh$BETA - 1.96 * values_tresh$SE  
  values_tresh$BETA_high <- values_tresh$BETA + 1.96 * values_tresh$SE  
  values_tresh$OR <- signif(exp(values_tresh$BETA),3)  
  values_tresh$OR_low <-  signif(exp(values_tresh$BETA_low),3)
  values_tresh$OR_high <-  signif(exp(values_tresh$BETA_high),3)
  
  values_tresh$OR_signif <- signif(values_tresh$OR,2)
  values_tresh$CMRB <- rep(Disease,4)
  values_tresh$Population <- rep(ifelse(subset==FALSE,"GP",subset),4)
  
  values_tresh
}
###Plot PRSs
plot_PRS_OR <- function(data) {
    plot <- ggplot(data,aes(y=as.factor(CMRB),x=as.numeric(OR_signif ), color=Population), group=Population) + 
        geom_pointrange(aes(xmin=as.numeric(OR_low ),xmax=as.numeric(OR_high )),
                    position=position_dodge(width=0.2),
                    size = 0.8) +
        geom_vline(xintercept = 1.0) + 
        theme_bw() + 
        xlab("OR")  +
        ylab("") +
        labs(colour=paste("Subset cohort")) +
        scale_y_discrete(limits=c("EPI","MH","MDD","T2D"),
                     breaks=c("EPI","MH","MDD","T2D"))  + theme_pubclean()
    plot

    plot1 <- plot + theme(
        axis.text.y  = element_text(size=30),
        axis.text.x  = element_text(size=30), 
        axis.title = element_text(size=30)) +
        theme(legend.position = "none") + 
        scale_color_brewer(palette = "Dark2") 
    p1 <- label_plot(plot1,"A",25) 
    p1
}
#plot PRS VS PRS 
PRS_PRS <- function(){
    data <- fread("ALL_PRS_pheno_covar.tsv")[,c(1,27,4,8,10,12,14,16,34)]
    data <- data[data$AGE >= 50]
    
    data_PD <- data[data$PD == 1 ]
    colnames(data_PD) <- c("FID","AGE","PGS EPI (zSCORE)","PGS T2D (zSCORE)","PGS Height (zSCORE)","PGS MDD (zSCORE)","PGS MH (zSCORE)","PGS Insomnia (zSCORE)","PD")
    
    
    plot2 <- ggpairs(data_PD, columns= c("PGS T2D (zSCORE)", "PGS MDD (zSCORE)","PGS MH (zSCORE)","PGS EPI (zSCORE)"), 
                     title = "",
                     upper = list(continuous = wrap("cor",color="black", size = 5)),
                     lower = list(continuous = wrap("smooth",color="gray", size = 2.5))) +
        theme_bw() 
    p2 <- label_plot(plot2,"B",25)
    p2
}
#make plots
{
OR_Pvalue_all <- data.frame(NULL)
age=50
for (i in list(c("T2D","T2D"),c("EPI","EPI"),c("MDD","MDD"),c("MH","MH"))) {
  OR_Pvalue_all <- rbind(OR_Pvalue_all,PRS_OR(i[1],i[2],FALSE,age))
  OR_Pvalue_all <- rbind(OR_Pvalue_all,PRS_OR(i[1],i[2],"PD",age))
}
OR_Pvalue_all
OR_Pvalue <- filter(OR_Pvalue_all,Percentil==20) } # Run OR_Pvalue_all() in all cmrbs
p1 <- plot_PRS_OR(OR_Pvalue)
p2 <- PRS_PRS()
  
#Save Plots
ggsave(file="All_figures/Figure3_A.pdf",p1,width = 10,height = 7)
ggsave(file="All_figures/Figure3_B.pdf",p2,width = 7,height = 6.5)
  

#Extract info to make table 2 and sup2
OR_Pvalue_all$Ci_95 <- paste0(OR_Pvalue_all$OR_low," - ",OR_Pvalue_all$OR_high )
TABLE <- OR_Pvalue_all[,c("Population","CMRB","Percentil","OR","Ci_95","p.value","High_PRS","Low_PRS")]
TABLE$p.value <- signif(TABLE$p.value, 3)

Table_2 <- TABLE[TABLE$Population == "PD",]
fwrite( Table_2 , "Tables/Table_2.tsv",sep="\t" )

Table_SUP2 <- TABLE[TABLE$Population == "GP",]
fwrite( Table_SUP2 , "Tables/TableSup2.tsv",sep="\t" )


