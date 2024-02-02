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
library(pROC)
library(magrittr)


############################### Figure sup 8 ,ROC curves

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
###### ROC curves
ROC_AUC <- function (PRS,Disease,subset,age=0){
  DATA <- load_data(PRS)
  subset_carga <- as.character(subset)
  if (subset == FALSE) {
    data <- DATA
  } else {
    data <- DATA[DATA[[subset_carga]] == 1 ]
  } 
  
  ## divide cohort in training and validation
  all_rows <- nrow(data)
  training <- all_rows * 0.8
  set.seed(123456)
  pheno <- sample_n(data, training )
  pheno_validation <- subset(data, !(FID %in% pheno$FID))
  pheno_validation$PHENO <- pheno_validation[[Disease]] 
  
  # We can then calculate the  model (with PRS) using a logistic regression 
  formula <- as.formula(paste0(Disease," ~ zSCORE + SEX + AGE + PC1 + PC2 + PC3 + PC4 +TD "))
  model_glm <- glm(formula, data=pheno, family = "binomial")
  # predict CMRB with model
  predicted <- predict(model_glm, pheno_validation, type="response")
  # Generate curva roc
  auc <- paste0("AUC=",round(auc(pheno_validation[[Disease]], predicted),2))
  rocobj <- roc(pheno_validation[[Disease]], predicted)
  rocobj
} 
#plot previous
PLOT_ROC <- function(PRS,Disease,subset,age=0) {
rocobj <- list(Data1 = ROC_AUC(PRS,Disease,FALSE),
               Data2 = ROC_AUC(PRS,Disease,subset))

auc_1 <- paste0("AUC=",signif(auc(rocobj[[1]]),2))
auc_2 <- paste0("AUC=",signif(auc(rocobj[[2]]),2))

   ##############################################drawing plot
  p <- ggroc(rocobj,aes=c("linetype", "color"),linetype = 1,size = 0.6,alpha=1)+
    labs(x = "1-Specificity", y = "Sensitivity") + 
    geom_abline(intercept = 1, slope = 1, color='grey',size = 0.5,linetype = "dashed") + 
    scale_y_continuous(expand = c(0.01, 0.001),breaks= seq(from = 0 , to = 1 ,by = 0.2))+
    scale_x_reverse(expand = c(0.01, 0.001),breaks= seq(from = 0 , to = 1 ,by = 0.2),limits = c(1,-0.01),
                    labels = c("1.0","0.8","0.6","0.4","0.2","0.0")) +
    scale_color_brewer(palette = "Dark2",
                       labels = c(paste0("GP - ",auc_1), paste0("PD - ",auc_2))) +
    theme_bw() + 
    theme(legend.position = c(0.6, 0.15),
          legend.background = element_rect(fill="white",
                                           size=0.4, linetype="solid", 
                                           colour ="gray30"),
          legend.title = element_blank(),
          legend.text=element_text(size=23,family = "sans"),
          #axis.text.x = element_blank(),
          axis.text.x = element_text(colour="black",size=20,face="plain",family = "sans"),
          axis.text.y = element_text(colour="black",size=20,face="plain",family = "sans"),
          axis.title.x = element_text(colour="black",size=25,face="plain",family = "sans"),
          axis.title.y = element_text(colour="black",size=25,face="plain",family = "sans"),
          plot.title = element_text(size=25))+
    ggtitle(paste0("Roc Curves for ",Disease))
  p
}

#make plots
p1 <- PLOT_ROC("T2D","T2D","PD") %>% label_plot(.,"A",25)
p2 <- PLOT_ROC("MDD","MDD","PD") %>% label_plot(.,"B",25)
p3 <- PLOT_ROC("MH","MH","PD") %>% label_plot(.,"C",25)
p4 <- PLOT_ROC("EPI","EPI","PD") %>% label_plot(.,"D",25)

##grid
lay <- rbind(c(1,2),
             c(3,4))
grid_landscape <- arrangeGrob(grobs = list(p1,p2,p3,p4), layout_matrix = lay)
#   15 *15
ggsave(file="All_figures/FigureSup8.pdf",grid_landscape,width = 15,height = 15)



