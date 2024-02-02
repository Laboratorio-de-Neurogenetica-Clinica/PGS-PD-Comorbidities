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


############################### Figure sup 6 ,PRS in age at onset linear asociation

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
#Plot PRS Vs AGE at onset continuos
LM_age_zSCORE <- function (PRS,Disease,subset,age=0) {
  data_test <- load_data(PRS)
  subset_carga <- as.character(subset)
  Disease <- as.character(Disease)
  
  
  ## cargar AAO PD
  AAO <-  load_data(PRS)[,c(1,25,26)]  
  colnames(AAO) <- c("FID","Birth","Year_PD")
    AAO$Year_PD <- as.numeric(format(AAO$Year_PD,'%Y'))
  AAO$AAO_PD <- AAO$Year_PD - AAO$Birth
  
  ##Filter by Age
  DATA <- data_test[data_test$AGE >= age]
  
  ### Filter PD 
  data_PD <- DATA[DATA[[subset_carga]] == 1 ]
  
  #Merge data
  data_PD_AAO <- merge(data_PD,AAO, by="FID")
  # rename
  data <- data_PD_AAO
  ##Filter by Age
  data_test <- data[data$AGE >= age]

    ggplotRegression <- function (fit) {
    
    require(ggplot2)
    
    ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
      geom_jitter(size= 1,
                 ) +
      stat_smooth(method = "lm", col = "red") +
      labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 3),
                         "P =",signif(summary(fit)$coef[2,4], 3))) 
  }
  formula <- as.formula("zSCORE ~ AAO_PD + SEX + AGE + TD + PC1 + PC2 +PC3 +PC4 ")
  regression <- lm (formula , data_test)
  p <- ggplotRegression(regression) + theme_bw() +
    theme(plot.title = element_text(size=15),
          axis.title = element_text(size = 20),
          axis.text = element_text(size=15)) +
    xlab("AAO PD") +
    ylab(paste0(PRS," PGS"))
  p
}


#make plots
p1 <- LM_age_zSCORE("T2D","T2D","PD",age) %>% label_plot(., "A",25)
p2 <- LM_age_zSCORE("MDD","MDD","PD",age) %>% label_plot(., "B",25)
p3 <- LM_age_zSCORE("MH","MH","PD",age) %>% label_plot(., "C",25)
p4 <- LM_age_zSCORE("EPI","EPI","PD",age) %>% label_plot(., "D",25)

##grid
lay <- rbind(c(1,2),
             c(3,4))
grid_landscape <- arrangeGrob(grobs = list(p1,p2,p3,p4), layout_matrix = lay)
#landscape 20*10
ggsave(file="All_figures/FigureSup6.pdf",grid_landscape,width = 20,height = 10)



