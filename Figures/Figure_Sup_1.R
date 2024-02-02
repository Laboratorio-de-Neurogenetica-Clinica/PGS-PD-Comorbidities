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
library(magrittr)

############################### Figure Sup 1 , DECILE PLOTS

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
decil_plot <- function(PRS,Disease,subset=FALSE,age=0) {
  DATA <- load_data(PRS)
  subset_carga <- as.character(subset)
  if (subset == FALSE) {
      data_test <- DATA
  } else {
      data_test <- DATA[DATA[[subset_carga]] == 1 ]
  }
  ##Filter by Age
  data <- data_test[data_test$AGE >= age]
  
  
  data$dec1 <- ifelse(data$zSCORE <= quantile(data$zSCORE,0.1), 1, 0)
  data$dec2 <- ifelse(data$zSCORE > quantile(data$zSCORE,0.1) & data$zSCORE <= quantile(data$zSCORE,0.2), 1, 0)
  data$dec3 <- ifelse(data$zSCORE > quantile(data$zSCORE,0.2) & data$zSCORE <= quantile(data$zSCORE,0.3), 1, 0)
  data$dec4 <- ifelse(data$zSCORE > quantile(data$zSCORE,0.3) & data$zSCORE <= quantile(data$zSCORE,0.4), 1, 0)
  data$dec5 <- ifelse(data$zSCORE > quantile(data$zSCORE,0.4) & data$zSCORE <= quantile(data$zSCORE,0.5), 1, 0)
  data$dec6 <- ifelse(data$zSCORE > quantile(data$zSCORE,0.5) & data$zSCORE <= quantile(data$zSCORE,0.6), 1, 0)
  data$dec7 <- ifelse(data$zSCORE > quantile(data$zSCORE,0.6) & data$zSCORE <= quantile(data$zSCORE,0.7), 1, 0)
  data$dec8 <- ifelse(data$zSCORE > quantile(data$zSCORE,0.7) & data$zSCORE <= quantile(data$zSCORE,0.8), 1, 0)
  data$dec9 <- ifelse(data$zSCORE > quantile(data$zSCORE,0.8) & data$zSCORE <= quantile(data$zSCORE,0.9), 1, 0)
  data$dec10 <- ifelse(data$zSCORE > quantile(data$zSCORE,0.9), 1, 0)
  
  data$dec <- 1
  data$dec[data$dec2 == 1] <- 2
  data$dec[data$dec3 == 1] <- 3
  data$dec[data$dec4 == 1] <- 4
  data$dec[data$dec5 == 1] <- 5
  data$dec[data$dec6 == 1] <- 6
  data$dec[data$dec7 == 1] <- 7
  data$dec[data$dec8 == 1] <- 8
  data$dec[data$dec9 == 1] <- 9
  data$dec[data$dec10 == 1] <- 10
  
  
  
  
  formula <- paste0(Disease, " ~ ","as.factor(data$dec) + SEX + PC1 + PC2 + PC3 + PC4")
  decTests <- glm(as.formula(formula), family="binomial", data=data)
  
  
  summary(decTests)
  summary_stats <- data.frame(summary(decTests)$coef[2:10,1:2])
  names(summary_stats) <- c("BETA","SE")
  summary_stats$DECILE <- c("2nd","3rd","4th","5th","6th","7th","8th","9th","10th")
  summary_stats$Pvalue <- summary(decTests)$coef[2:10,4]
  summary_stats$Pvalue <- signif(summary_stats$Pvalue, 2)
  summary_stats[10,] <- c(0,0,"1st",NA)
  summary_stats$decil_sort <- c(2,3,4,5,6,7,8,9,10,1)
  summary_stats_sorted <- summary_stats[order(summary_stats$decil_sort),]
  
  write.table(summary_stats_sorted, "TEMP.csv", quote = F, row.names = F, sep = ",")
  
  
  to_plot <- read.table("TEMP.csv", header = T, sep = ",")
  to_plot$low <- to_plot$BETA - 1.96*to_plot$SE
  to_plot$high <- to_plot$BETA + 1.96*to_plot$SE
  to_plot$low_OR <- exp(to_plot$low) 
  to_plot$high_OR <- exp(to_plot$high)
  to_plot$OR <- exp(to_plot$BETA)
  p <- ggplot(to_plot, aes(as.factor(decil_sort), OR)) + 
    geom_errorbar(aes(ymin = low_OR, ymax = high_OR), width=0.3) +
    geom_point(size=3.5) +
    scale_x_discrete(labels=c("1"="1th","2"="2th","3"="3th","4"="4th","5"="5th","6"="6th","7"="7th","8"="8th","9"="9th","10"="10th")) +
    geom_text(aes(x=decil_sort, y=OR,label=Pvalue),
              angle = 45,
              vjust= -0.8,
              hjust= 0.9)
  p  
  p2 <- p +   
    ylab(paste0("OR for ",Disease," in ",ifelse(subset==FALSE,"General Population",subset)))  +
    xlab(paste0(PRS," PGS decile")) +
    theme_bw() + 
    theme(axis.text.y  = element_text(size=20),
          axis.text.x  = element_text(size=20), 
          axis.title = element_text(size=25))
  
  p2
}

#Make plots
p1 <- decil_plot("T2D","T2D","PD",age=0)  %>% label_plot(., "A", 25)
p2 <- decil_plot("MDD","MDD","PD")  %>% label_plot(., "B", 25)
p3 <- decil_plot("MH","MH","PD")  %>% label_plot(., "C", 25)
p4 <- decil_plot("EPI","EPI","PD",age=0) %>% label_plot(., "D", 25)

##Arrenge in grid
lay <- rbind(c(1,2),
             c(3,4))
grid_landscape <- arrangeGrob(grobs = list(p1,p2,p3,p4), layout_matrix = lay)
#landscape 17*12
ggsave(file="All_figures/FigureSup1.pdf",grid_landscape,width = 17,height = 12)


