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

############################### Figure sup 3 , PD risk in Disease, VIOLINPLOTS, DENSITY PLOT y DECILE PLOT

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
#Plot PRS Violin Plots
Violin_plots <- function (PRS,Disease,subset,age=0,test="t.test") {
  DATA <- load_data(PRS)
  subset_carga <- as.character(subset)
  Disease <- as.character(Disease)
  
  if (subset == "PD") {
      ### Filter PD 
      data_subset <- DATA[DATA[[subset_carga]] == 1 ]      
  } else {
      ### Filter CMRB 
      data_subset <- DATA[DATA[[Disease]] == 1 ]
      
  }
  #add cmrb statues without PD
  #data_cmrb <- DATA[DATA[[subset_carga]] == 0 & DATA[[Disease]] ==1  ]
  #data <- rbind(data_cmrb, data_PD) 
  
  

  ##Filter by Age
  data_test <- data_subset[data_subset$AGE >= age]
  ## Generate 3 cohorts and count
  if (subset == "PD") {
      ### with and withou cmrb 
      data_test$PHENO_PLOT <- as.factor(ifelse(data_test[[Disease]] == 0 , 0, 1))
  } else {
      ### without PD 
      data_test$PHENO_PLOT <- as.factor(ifelse(data_test$PD == 1 , 0, 1))
      
  }
  N_records <- data_test %>% 
    group_by(PHENO_PLOT) %>% 
    summarise(count = n(),MEAN=mean(zSCORE), MAX=max(zSCORE))
  
  ## pvalues 
  if (test == "t.test") {
  stat.test <- data_test %>% t_test(zSCORE ~ PHENO_PLOT)
  stat.test <- stat.test %>% add_xy_position(x = "PHENO_PLOT")
  stat.test$p.adj <- stat.test$p
  
  } else if (test == "anova"){
    stat.test <- data_test %>% tukey_hsd(zSCORE ~ PHENO_PLOT)
    stat.test <- stat.test %>% add_xy_position(x = "PHENO_PLOT")
  }
  

  ##plot Violin Plot
  p <- ggplot(data_test,aes_string(x= reorder(as.factor(data_test$PHENO_PLOT),data_test$zSCORE),
                                   y="zSCORE",
                                   fill=as.factor(Disease))) +
      scale_x_discrete(
          #limits=c("1","0"),
          #breaks=c("0","1"),
          labels= c(ifelse(subset=="PD",
                           subset,
                           paste0(Disease,"\nwith PD")),
                    ifelse(subset=="PD",
                           paste0(Disease,"\nwith PD"),
                           paste0("Only\n",Disease)
                           ))
      )
  
  p
  p2 <- p+
      geom_text(data= N_records, aes(x=PHENO_PLOT ,y=MEAN,label=paste0("N=",count)),
                vjust= -6, hjust=- 0.4 ,size=6 ) +
    geom_violin(aes_string(fill = as.factor(data_test$PHENO_PLOT)),
                width=0.5) +
    geom_boxplot(width=0.15, 
                 fill="white",
                 outlier.shape = NA) +
    theme_minimal() +
    stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.005, 
                       step.increase = 0.02,
                       size = 6) +
      ylim(min(data_test$zSCORE),
           max(data_test$zSCORE) + 0.5)
  p2
  p3 <- p2 + scale_fill_manual(values=c(ifelse(subset=="PD",
                                             "#7FC97F",
                                             "#BEAED4"),
                                        ifelse(subset=="PD",
                                               "#BEAED4",
                                               "orange"),
                                        "black")) +
    
    ylab(paste0(PRS," PGS (zSCORE)"))  +
    xlab(paste0(ifelse(age != 0 , paste0("Patients Over ",age," Years") ,""))) + theme(legend.position = "none") + 
    theme(axis.text.x = element_text(size=20),
          axis.title = element_text(size=20),
          axis.text.y= element_text(size=15))
  
  p3
}
#Histogram
PRS_histogram <- function(PRS,Disease,subset=FALSE,age=0){
    DATA <- load_data(PRS)
    subset <- as.character(subset)
    Disease <- as.character(Disease)
  
    ### Filter  
    data_subset <- DATA[DATA[[subset]] ==1 ]
      
  
  
  ##Filter by Age
  data <- data_subset[data_subset$AGE >= age]
  data[[Disease]] <- as.character(data[[Disease]])
  data[[subset]] <-as.character(data[[subset]])
  
  
  TEST <-  ks.test(x= data$zSCORE[data[[Disease]] == 0]  ,
                   y= data$zSCORE[data[[Disease]] == 1] )
  stat.test <- c(pvalue=signif(TEST[[2]],3))
  histogram1 <- density(data$zSCORE[data[[Disease]] == 0])
  histogram2 <- density(data$zSCORE[data[[Disease]] == 1])
  
  max_height1 <- max(histogram1$y)
  max_height2 <- max(histogram2$y)
  maxvalue <- ifelse( max_height1 > max_height2,
                      max_height1,
                      max_height2)
  
  
  #plot
  p <- ggplot(data,aes_string("zSCORE", fill=Disease  )) + 
      geom_density(weight=10 ,alpha= 0.4 ) +
      geom_vline(xintercept = mean(data$zSCORE[data[[Disease]]==0]),
                 alpha=1, linewidth=1,
                 linetype="longdash",
                 color=ifelse(subset=="PD",
                            "#7FC97F",
                            "orange") ) +
      geom_vline(xintercept = mean(data$zSCORE[data[[Disease]]==1]),
                 alpha=1, linewidth=1,
                 linetype="longdash",color="#BEAED4" ) +
      annotate("text", x = 2, y = maxvalue, label = stat.test, size=6)
  p2 <-  p + 
    ylab(ifelse(subset==Disease, "General Population (%)",paste0(subset," Subset (%)")))  +
    xlab(paste0(PRS," PGS (zSCORE)")) +
    scale_fill_manual(name = "Status",
                      values=c(ifelse(subset=="PD",
                                      "#7FC97F",
                                      "orange"),
                               "#BEAED4"),
                      labels= c(ifelse(subset=="PD",
                                       subset,
                                       paste0("Only\n",subset)),
                                paste0(subset,"\nwith PD"))
                      )  +
    labs(fill="Statues") +
    theme_bw() + 
    theme(legend.position = c(0.15,0.85),
          legend.title = element_text(color = "black", size = 20,),
          legend.text = element_text(color = "black", size = 15),
          legend.background = element_rect(fill = "white", colour = "gray30"),
          legend.key.height=unit(1, "cm"),
          legend.key.width = unit(1,"cm"),
          axis.text = element_text(size=20),
          axis.title = element_text(size=20),) 
    
  p2
}
#DecilePlot
decil_plot <- function(PRS,Disease,subset=FALSE,age=0) {
    DATA <- load_data(PRS)
    subset <- as.character(subset)
    Disease <- as.character(Disease)
    
    ### Filter  
    data_subset <- DATA[DATA[[subset]] ==1 ]
    
    
    
    ##Filter by Age
    data  <- data_subset[data_subset$AGE >= age]
    
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
        theme(axis.text.y  = element_text(size=30),
              axis.text.x  = element_text(size=25), 
              axis.title = element_text(size=30))
    
    p2
}

#make plots
p1 <- Violin_plots("T2D","T2D",subset=FALSE,test ="anova",age=0) %>% label_plot(.,"A",25)
p2 <- Violin_plots("MDD","MDD",subset=FALSE,test = "t.test",age=0) %>% label_plot(.,"B",25)
out <- grid.rect(gp=gpar(col="white"))
p3 <- Violin_plots("MH","MH",subset=FALSE,test = "t.test",age=0) %>% label_plot(.,"C",25)
p4 <- Violin_plots("EPI","EPI",subset=FALSE,test = "t.test",age=0) %>% label_plot(.,"D",25)
p5 <- PRS_histogram("T2D","PD","T2D",age=0) %>% label_plot(.,"E",25)
p6 <- PRS_histogram("MDD","PD","MDD",age=0) %>% label_plot(.,"F",25)
out <- grid.rect(gp=gpar(col="white"))
p7 <- PRS_histogram("MH","PD","MH",age=0) %>% label_plot(.,"G",25)
p8 <-PRS_histogram("EPI","PD","EPI",age=0) %>% label_plot(.,"H",25)
p9 <- decil_plot("T2D","PD","T2D",age=0) %>% label_plot(.,"I",25)
p10 <- decil_plot("MDD","PD","MDD",age=0) %>% label_plot(.,"J",25)
out <- grid.rect(gp=gpar(col="white"))
p11 <- decil_plot("MH","PD","MH",age=0) %>% label_plot(.,"K",25)
p12 <- decil_plot("EPI","PD","EPI",age=0) %>% label_plot(.,"L",25)



######grid ######
lay_vertical <- rbind(c(1,5,9),
                      c(2,6,10),
                      c(3,7,11),
                      c(4,8,12))
grid_vertical <- arrangeGrob(grobs = list(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12), layout_matrix = lay_vertical)
### PDF Horizontal 22.5*30
ggsave(file="All_figures/FigureSup4.pdf",grid_vertical,width = 22.5,height = 30)

