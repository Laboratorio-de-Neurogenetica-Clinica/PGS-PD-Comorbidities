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
 
############################### Figure 2 , Violin plot and density plots

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
###### Plot PRS Violin Plots ####
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
                             paste0(Disease,"\n with PD"),
                             paste0("Only\n",Disease)
                      ))
        )
    
    p
    p2 <- p+
        geom_text(data= N_records, aes(x=PHENO_PLOT ,y=MEAN,label=paste0("N=",count)),
                  vjust= -6, hjust=- 0.3 ,size=6 ) +
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
        theme_bw() + 
        ylab(paste0(PRS," PGS (zSCORE)"))  +
        xlab(paste0(ifelse(age != 0 , paste0("Patients Over ",age," Years") ,""))) + theme(legend.position = "none") + 
        theme(axis.text.x = element_text(size=20),
              axis.title = element_text(size=20),
              axis.text.y= element_text(size=15))
    
    p3
}
###### Plot PRS Histogram ####
PRS_histogram <- function(PRS,Disease,subset=FALSE,age=0,test="ks.test"){
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
                                    paste0(Disease,"\nwith PD"))
        )  +
        labs(fill="Statues") +
        theme_bw() + 
        theme(legend.position = c(0.20,0.80),
              legend.title = element_text(color = "black", size = 20,),
              legend.text = element_text(color = "black", size = 15),
              legend.background = element_rect(fill = "white", colour = "gray30"),
              legend.key.height=unit(1, "cm"),
              legend.key.width = unit(1,"cm"),
              axis.text = element_text(size=20),
              axis.title = element_text(size=20),) 
    
    p2
}

#make Plots
p1 <- Violin_plots("T2D","T2D","PD",test ="t.test",age=50) %>% label_plot(.,"A",25)
p2 <- Violin_plots("MDD","MDD","PD",test ="t.test",age=50) %>% label_plot(.,"B",25)
p3 <- Violin_plots("MH","MH","PD",test ="t.test",age=50)  %>% label_plot(.,"C",25)
p4 <- Violin_plots("EPI","EPI","PD",test = "t.test",age=50)  %>% label_plot(.,"D",25)
p5 <- PRS_histogram("T2D","T2D","PD",age=50)  %>% label_plot(.,"E",25) 
p6 <- PRS_histogram("MDD","MDD","PD",age=50)  %>% label_plot(.,"F",25)
p7 <- PRS_histogram("MH","MH","PD",age=50)  %>% label_plot(.,"G",25)
p8 <- PRS_histogram("EPI","EPI","PD",age=50) %>% label_plot(.,"H",25)


##Arrenge in grid
lay_vertical <- rbind(c(1,5),
             c(2,6),
             c(3,7),
             c(4,8))

grid_vertical <- arrangeGrob(grobs = list(p1,p2,p3,p4,p5,p6,p7,p8), layout_matrix = lay_vertical)

### PDF Horizontal 20*11
ggsave(file="All_figures/Figure2.pdf",grid_vertical,width = 11,height = 20)
