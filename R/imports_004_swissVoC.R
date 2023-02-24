#' ---
#' title: Imports and the Swiss SARS-CoV-2 Epidemic
#' subtitle: model spread of new variants in Switzerland
#' author: "Martina Reichmuth"
#' date: "'01/12/2021"
#' ---

### Load data
#Rows that have no entry for 'Clade' nor 'HeadOfClade' (end of file) are 'separate' nodes like b) below - they don't have any Swiss sequences above them. These are very likely introductions.
#Rows that have an entry for 'Clade' and 'HeadOfClade' are the first node in a chain, like node 1554 in a) below. They don't have any Swiss sequences above them. These are also likely introductions.
#Rows that have an entry for 'Clade' but not 'HeadOfClade' are the other nodes in a chain, like the other nodes in a) below. They can either be interpreted as further introductions, or not, depending on how you want to count.
# Alpha imports from Emma Hodcroft
imports_conservative_alpha <- read.csv("./data/conservativeClusters-Alpha.csv")[,-1]
imports_alpha <-read.csv("./data/liberalClusters-Alpha.csv")[,-1]
imports_alpha$MinDateSwissChildren[grepl("2020-01-05", imports_alpha$MinDateSwissChildren)] <- ("2020-12-26")
imports_conservative_alpha$MinDateSwissChildren[grepl("2020-01-05", imports_conservative_alpha$MinDateSwissChildren)] <- ("2020-12-26")
imports_alpha$approach <- ifelse(imports_alpha$HeadOfClade!="True" & is.na(imports_alpha$Clade)| imports_alpha$HeadOfClade=="True" & !is.na(imports_alpha$Clade)  & imports_alpha$Clade>-1,"conservative", "liberal" )
imports_alpha$MinDateNonSwissChildren <- as_date(imports_alpha$MinDateNonSwissChildren)

# Delta imports from Emma Hodcroft
imports_conservative_delta <- read.csv("./data/conservativeClusters-Delta.csv")[,-1]
imports_delta <-read.csv("./data/liberalClusters-Delta.csv")[,-1]
imports_delta <- imports_delta[!grepl("2021-02-13", imports_delta$MinDateSwissChildren),]
imports_delta$approach <- ifelse(imports_delta$HeadOfClade!="True" & is.na(imports_delta$Clade)| imports_delta$HeadOfClade=="True" & !is.na(imports_delta$Clade)  & imports_delta$Clade>-1,"conservative", "liberal" )
imports_delta$MinDateNonSwissChildren <- as_date(imports_delta$MinDateNonSwissChildren)

imports_files <- list(imports_conservative_alpha, imports_alpha,imports_conservative_delta, imports_delta)

for (n in 1:length(imports_files)) {
  imports <- imports_files[[n]]
  if(n %in% c(1,2)){
    variant_data <- alpha_ch
    variant <- "alpha"
  }
  if(n %in% c(3,4)){
    variant_data <- delta_ch
    variant <- "delta"
  }
  cov_ch <- subset(swiss_cov, as_date(date) %in% seq(min(variant_data$date),max(variant_data$date),1))
  time_window <- c(min(variant_data$date),max(variant_data$date))
  period <- seq(min(variant_data$date),max(variant_data$date),1)
  times_length <- length(period)
  
  imports <- (as.numeric(as_date(imports$MinDateSwissChildren)) - as.numeric(min(variant_data$date)))
  imports <-  imports-7 # delay to getting tested
  imports <- hist(imports, 0:times_length, plot = FALSE)$counts
  
  sim <- imports_003_SEIR(imports,variant)
  if(n ==1){
    sim_1 <- sim
    sim_1$approach <- "conservative"
  }
  else if(n ==2){
    sim_2 <- sim
    sim_2$approach <- "liberal"
  }
  else if(n ==3){
    sim_3 <- sim
    sim_3$approach <- "conservative"
  }
  else if(n ==4){
    sim_4 <- sim
    sim_4$approach <- "liberal"
  }
}

for (n in 1:2) {
  if(n %in% c(1)){
    labels_variant <- c("All cases","Alpha cases", "other cases")
    titel_variant <- bquote("Proportion of Alpha\n")
    variant_name <- "Alpha cases" 
    variant_data <- alpha_ch
    cov_ch <- subset(swiss_cov, as_date(date) %in% seq(min(variant_data$date),max(variant_data$date),1))
    time_window <- c(min(variant_data$date),max(variant_data$date))
    period <- seq(min(variant_data$date),max(variant_data$date),1)
    sim<- rbind(sim_1,sim_2)
    sim_alpha <- sim
    imports_all <- imports_alpha
  }
  if(n %in% c(2)){
    labels_variant <- c("All cases","Delta cases", "other cases")
    titel_variant <- bquote("Proportion of Delta\n")
    variant_name <- "Delta cases" 
    variant_data <- delta_ch
    cov_ch <- subset(swiss_cov, as_date(date) %in% seq(min(variant_data$date),max(variant_data$date),1))
    time_window <- c(min(variant_data$date),max(variant_data$date))
    period <- seq(min(variant_data$date),max(variant_data$date),1)
    sim<- rbind(sim_3,sim_4)
    sim_delta <- sim
    imports_all <- imports_delta
  }
  imports_all$MinDateSwissChildren <- as_date(imports_all$MinDateSwissChildren)-7
  imports_all$approach <- factor(imports_all$approach, levels =  c("liberal","conservative"))
  sim$approach <- factor(sim$approach, levels =  c("liberal","conservative"))
  
  # Plot
  imports_hist_plot <- ggplot()+
    geom_bar(imports_all, mapping=aes(x=MinDateSwissChildren, fill=approach,color=approach), position="stack",stat = "count")+
    scale_x_date(date_breaks = "1 month", 
                 date_labels = "%b", limits = c(time_window[1]-1,time_window[2]+1))+
    scale_y_continuous(limits=c(0,100))+
    theme_minimal()+
    scale_fill_manual(values= col_9[c(2,4)],name=variant_name, labels= c("liberal","conservative")) +
    scale_color_manual(values= col_9[c(2,4)],name=variant_name, labels= c("liberal","conservative")) +
    theme(legend.position=c(.2,.7),
          plot.subtitle = element_text(hjust = 0.5),
          axis.title.y = element_text(size = 10),
          title = element_text(size = 10))+
    labs(tag="",subtitle = bquote(), x = "", y ="Number of imports \n")
  
  
  proportion_plot_Re_t <- ggplot(data= sim)+
    geom_errorbar(data= variant_data, aes(x = date, ymin=lower, ymax=upper), color= col_9[9], alpha=0.5,width=.1) +
    geom_point(data= variant_data,aes(x = date, y=proportion), color= col_9[9])+
    geom_line(aes(x=date, y=T2/T_total, color=approach), size=1)+
    scale_color_manual(values= col_9[c(2,4)],name=variant_name, labels= c("liberal","conservative")) +
    theme_minimal()+
    theme(legend.position = "none")+
    scale_y_continuous(limits = c(0,1))+
    theme(plot.subtitle = element_text(hjust = 0.5),
          axis.title.y = element_text(size = 10),
          title = element_text(size = 10))+
    scale_x_date(date_breaks = "1 month", 
                 date_labels = "%b", limits = c(time_window[1]-1,time_window[2]+1))+
    labs(tag="",subtitle = bquote(), x = "", y =paste0(titel_variant,"\n"))#mu_wt=1d^-1 #.(imports)~" import per day"
  
  
  cases_day_plot_Re_t <-  ggplot(data= sim)+
    geom_line(data = cov_ch, aes(x = date, y = weigthed_cases),color=col_9[9], size=1)+
    geom_line(aes(x=date, y=T_total,group=approach), alpha=0.6,size=1)+
    geom_line(aes(x=date, y=T2, color=approach), alpha=0.6,size=1)+
    scale_color_manual(values= col_9[c(2,4)],name=variant_name, labels= c("liberal","conservative")) +
    theme_minimal()+
    scale_y_continuous(limits = c(0,8100))+
    scale_x_date(date_breaks = "1 month", 
                 date_labels = "%b", limits = c(time_window[1]-1,time_window[2]+1))+
    theme_minimal()+
    theme(legend.position = "none",
          plot.subtitle = element_text(hjust = 0.5),
          axis.title.y = element_text(size = 10),
          title = element_text(size = 10))+
    labs(tag="",subtitle = "", x = "", y =bquote("Number of reported \ncases per day"))
  
  if(n %in% c(1)){
    plot_model_ch_alpha <- ggarrange(imports_hist_plot, cases_day_plot_Re_t,proportion_plot_Re_t, 
                                     labels = c("A", "C", "E"),
                                     ncol = 1, nrow = 3)
    plot_model_ch_alpha1 <- ggarrange(imports_hist_plot,proportion_plot_Re_t, 
                                      labels = c("A", "C"),
                                      ncol = 1, nrow = 2)
    proportion_plot_Re_t_alpha <- ggplot(data= sim)+
      geom_errorbar(data= variant_data, aes(x = date, ymin=lower, ymax=upper), color= col_9[9], alpha=0.5,width=.1) +
      geom_point(data= variant_data,aes(x = date, y=proportion), color= col_9[9])+
      geom_line(aes(x=date, y=proportion, color=approach), size=1)+
      scale_color_manual(values= col_9[c(2,4)],name="", labels= c("liberal","conservative")) +
      theme_minimal()+
      scale_y_continuous(limits = c(0,1))+
      theme(legend.position=c(.2,.9),
            plot.subtitle = element_text(hjust = 0.5),
            axis.title.y = element_text(size = 10),
            title = element_text(size = 10))+
      scale_x_date(date_breaks = "1 month", 
                   date_labels = "%b", limits = c(time_window[1]-1,time_window[2]+1))+
      labs(tag="",subtitle = bquote(), x = "", y =paste0(titel_variant,"\n"))
  }
  if(n %in% c(2)){
    plot_model_ch_delta <- ggarrange(imports_hist_plot, cases_day_plot_Re_t,proportion_plot_Re_t, 
                                     labels = c("B", "D", "F"),
                                     ncol = 1, nrow = 3)
    plot_model_ch_delta1 <- ggarrange(imports_hist_plot,proportion_plot_Re_t, 
                                      labels = c("B", "D"),
                                      ncol = 1, nrow = 2)
    proportion_plot_Re_t_delta <- ggplot(data= sim)+
      geom_errorbar(data= variant_data, aes(x = date, ymin=lower, ymax=upper), color= col_9[9], alpha=0.5,width=.1) +
      geom_point(data= variant_data,aes(x = date, y=proportion), color= col_9[9])+
      geom_line(aes(x=date, y=proportion, color=approach), size=1)+
      scale_color_manual(values= col_9[c(2,4)],name="", labels= c("liberal","conservative")) +
      theme_minimal()+
      theme(legend.position = "none")+
      scale_y_continuous(limits = c(0,1))+
      theme(plot.subtitle = element_text(hjust = 0.5),
            axis.title.y = element_text(size = 10),
            title = element_text(size = 10))+
      scale_x_date(date_breaks = "1 month", 
                   date_labels = "%b", limits = c(time_window[1]-1,time_window[2]+1))+
      labs(tag="",subtitle = bquote(), x = "", y =paste0(titel_variant,"\n"))#mu_wt=1d^-1 #.(imports)~" import per day"
    
  }
}
plot_model_ch <- ggarrange(plot_model_ch_alpha, plot_model_ch_delta,
                           ncol = 2, nrow = 1)
ggsave(plot_model_ch, filename = paste0("./data/figures/Figure4_import_variants_",format(Sys.time(), "%Y-%m-%d"), ".pdf"), height = 7, width = 7,  bg = "transparent")
ggsave(plot_model_ch, filename = paste0("./data/figures/Figure4_import_variants_",format(Sys.time(), "%Y-%m-%d"), ".png"), height = 7, width = 7,  bg = "transparent")
ggsave(plot_model_ch, filename = paste0("./data/figures/Figure4_import_variants_",format(Sys.time(), "%Y-%m-%d")), height =9, width = 8,  bg = "transparent",device='tiff', dpi=700)

