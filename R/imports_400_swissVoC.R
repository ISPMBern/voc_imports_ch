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
alpha_files <- list.files(paste0(path_name,"data/Alpha/"))
conservative_alpha <- alpha_files[grepl("conservativeClustersDate", alpha_files)]
conservative_alpha <- conservative_alpha[order(as.numeric(gsub("[^0-9]+","\\1",conservative_alpha)))]
liberal_alpha <- alpha_files[grepl("liberalClustersDate", alpha_files)]
liberal_alpha <- liberal_alpha[order(as.numeric(gsub("[^0-9]+","\\1",liberal_alpha)))]
# Delta imports from Emma Hodcroft
delta_files <- list.files(paste0(path_name,"data/Delta/"))
conservative_delta <- delta_files[grepl("conservativeClustersDate", delta_files)]
conservative_delta <- conservative_delta[order(as.numeric(gsub("[^0-9]+","\\1",conservative_delta)))]
liberal_delta <- delta_files[grepl("liberalClustersDate", delta_files)]
liberal_delta <- liberal_delta[order(as.numeric(gsub("[^0-9]+","\\1",liberal_delta)))]

subsampling_i <- length(liberal_alpha)
lagdetection_l <- 7#c(3:13)
simulations <-c()
imports_files <- c()
imports_sum <- c()

for(l in length(lagdetection_l)){
  for(i in subsampling_i){
    imports_conservative_alpha <- read.csv(paste0(path_name,"data/Alpha/",conservative_alpha[i]))[,-1]
    imports_alpha <- read.csv(paste0(path_name,"data/Alpha/",liberal_alpha[i]))[,-1]
    imports_alpha$approach <- ifelse(imports_alpha$HeadOfClade!="True" & is.na(imports_alpha$Clade)| imports_alpha$HeadOfClade=="True" & !is.na(imports_alpha$Clade)  & imports_alpha$Clade>-1,"conservative", "liberal" )
    imports_alpha$MinDateSwissChildren[grepl("2020-01-05", imports_alpha$MinDateSwissChildren)] <- ("2020-12-26")
    imports_conservative_alpha$MinDateSwissChildren[grepl("2020-01-05", imports_conservative_alpha$MinDateSwissChildren)] <- ("2020-12-26")
    imports_alpha$MinDateNonSwissChildren <- as_date(imports_alpha$MinDateNonSwissChildren)
    
    imports_conservative_delta <- read.csv(paste0(path_name,"data/Delta/",conservative_delta[i]))[,-1]
    imports_delta <- read.csv(paste0(path_name,"data/Delta/",liberal_delta[i]))[,-1]
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
     
      clusters <- c()
      imports$TotalDatesSwissChildren <- gsub("[',]","",imports$TotalDatesSwissChildren)
      imports$TotalDatesSwissChildren <- gsub("\\[|\\]","",imports$TotalDatesSwissChildren)
      
      clusters$name <- (1:length(imports$TotalDatesSwissChildren))
      clusters$MinDateSwissChildren <- as_date(imports$MinDateSwissChildren)
      clusters <- as.data.frame(clusters)
      
      cluster <- melt(strsplit(imports$TotalDatesSwissChildren, split=" "))
      colnames(cluster) <- c("cluster_dates", "name")
      cluster$cluster_dates <- as_date(cluster$cluster_dates)
      
      cluster <- cluster[ cluster$cluster_dates>=min(clusters$MinDateSwissChildren),]
      clusters <- merge(cluster, clusters, by="name", all.y = T)
      clusters <- clusters[order(clusters$cluster_dates),]
      cluster <- clusters %>% group_by(name,cluster_dates) %>% summarise(MinDateSwissChildren=unique(MinDateSwissChildren), cluster_size = length(cluster_dates)) %>% mutate(cluster_size = cumsum(cluster_size))
      cluster <- cluster[!is.na(cluster$cluster_dates),]
                         
      cluster_size_plot <- ggplot()+
        theme_cowplot()+
        geom_line(data=cluster, aes(x=cluster_dates,y =cluster_size, color= as.character(name)))+
        scale_x_date(date_breaks = "1 month", date_labels = "%b")+ 
        theme(legend.position ="none")+
        coord_cartesian(ylim=c(0,max(cluster$cluster_size[cluster$cluster_dates<=min(cluster$cluster_dates)+60])),
                        xlim=c(min(cluster$cluster_dates),min(cluster$cluster_dates)+60))+
        labs(tag="",x = "", y ="Cluster size \n ")
        
      
    
      imports <- (as.numeric(as_date(imports$MinDateSwissChildren)) - as.numeric(min(variant_data$date)))
      imports_sum <- c(imports_sum,length(imports))
      imports <-  imports-lagdetection_l[l] # delay to getting tested
      imports <- hist(imports, 0:times_length, plot = FALSE)$counts
      
      sim <- imports_300_SEIR(imports,variant)
      
      sim$variant <- variant
      sim$lagdetection <- lagdetection_l[l]
      
      if(n==1){
        sim$approach <- "conservative"
        cluster_size_plot_conserv_alpha <- cluster_size_plot +labs(tag="A",y ="Alpha variant\nCluster size", subtitle = "conservative")
        cluster_size_plot_conserv_alpha <-ggdraw() +draw_plot(cluster_size_plot_conserv_alpha)+
          draw_plot(ggplot()+
                      theme_cowplot()+
                      geom_line(data=cluster, aes(x=cluster_dates,y =cluster_size, group=name), color= col_9[9])+
                      scale_x_date(date_breaks = "1 month", date_labels = "%b")+ # limits = c(min(variant_data$date),min(variant_data$date)+90)
                      theme(legend.position ="none", axis.text.y = element_text(size=rel(0.4)), axis.text.x = element_text(size=rel(0.3)))+
                      labs(tag="",x = "", y =paste0(" \n")), x = 0.1, y = .4, width = .55, height = .55)
        
         }
      if(n==2){
        sim$approach <- "liberal"
        sim_alpha_time_dominance <- min(sim$date[sim$proportion>=0.5])
        cluster_size_plot_liberal_alpha <- cluster_size_plot +labs(tag="B", subtitle = "liberal")
        cluster_size_plot_liberal_alpha <-ggdraw() +draw_plot(cluster_size_plot_liberal_alpha)+
          draw_plot(ggplot()+
                      theme_cowplot()+
                      geom_line(data=cluster, aes(x=cluster_dates,y =cluster_size, group=name), color= col_9[9])+
                      scale_x_date(date_breaks = "1 month", date_labels = "%b")+ # limits = c(min(variant_data$date),min(variant_data$date)+90)
                      theme(legend.position ="none", axis.text.y = element_text(size=rel(0.4)), axis.text.x = element_text(size=rel(0.3)))+
                      labs(tag="",x = "", y =paste0(" \n")), x = 0.1, y = .4, width = .55, height = .55)
      }
      if(n==3){
        sim$approach <- "conservative"
        cluster_size_plot_conserv_delta <- cluster_size_plot +labs(tag="C",y ="Delta variant\nCluster size")
        cluster_size_plot_conserv_delta <-ggdraw() +draw_plot(cluster_size_plot_conserv_delta)+
          draw_plot(ggplot()+
                      theme_cowplot()+
                      geom_line(data=cluster, aes(x=cluster_dates,y =cluster_size, group=name), color= col_9[9])+
                      scale_x_date(date_breaks = "1 month", date_labels = "%b")+ # limits = c(min(variant_data$date),min(variant_data$date)+90)
                      theme(legend.position ="none", axis.text.y = element_text(size=rel(0.4)), axis.text.x = element_text(size=rel(0.3)))+
                      labs(tag="",x = "", y =paste0(" \n")), x = 0.1, y = .4, width = .55, height = .55)
        
      }
      if(n==4){
        sim$approach <- "liberal"
        sim_delta_time_dominance <- min(sim$date[sim$proportion>=0.5])
        cluster_size_plot_liberal_delta <- cluster_size_plot +labs(tag="D")
        cluster_size_plot_liberal_delta <-ggdraw() +draw_plot(cluster_size_plot_liberal_delta)+
          draw_plot(ggplot()+
                      theme_cowplot()+
                      geom_line(data=cluster, aes(x=cluster_dates,y =cluster_size, group=name), color= col_9[9])+
                      scale_x_date(date_breaks = "1 month", date_labels = "%b")+ # limits = c(min(variant_data$date),min(variant_data$date)+90)
                      theme(legend.position ="none", axis.text.y = element_text(size=rel(0.4)), axis.text.x = element_text(size=rel(0.3)))+
                      labs(tag="",x = "", y =paste0(" \n")), x = 0.1, y = .4, width = .55, height = .55)
        
      }
       
       simulations[[length(simulations)+1]] <- sim
       sim <- NULL
      }
      
  }
  simulations_baseline <-simulations

for(m in 1:2) {
  if(m %in% c(1)){
    labels_variant <- c("All cases","Alpha cases", "other cases")
    titel_variant <- bquote("Proportion of Alpha\n")
    variant_name <- "Alpha cases" 
    variant_data <- alpha_ch
    cov_ch <- subset(swiss_cov, as_date(date) %in% seq(min(variant_data$date),max(variant_data$date),1))
    time_window <- c(min(variant_data$date),max(variant_data$date))
    period <- seq(min(variant_data$date),max(variant_data$date),1)
    variant <- "alpha"
  }
  if(m %in% c(2)){
    labels_variant <- c("All cases","Delta cases", "other cases")
    titel_variant <- bquote("Proportion of Delta\n")
    variant_name <- "Delta cases" 
    variant_data <- delta_ch
    cov_ch <- subset(swiss_cov, as_date(date) %in% seq(min(variant_data$date),max(variant_data$date),1))
    time_window <- c(min(variant_data$date),max(variant_data$date))
    period <- seq(min(variant_data$date),max(variant_data$date),1)
    variant <- "delta"
  }
  
  sim <- as.data.frame(do.call(rbind, lapply(simulations_baseline, `length<-`, max(lengths(simulations_baseline)))))
  sim <- sim[sim$variant==variant,]
  sim$Imports[is.na(sim$Imports)] <- 0
  time_window_imports <- max(sim$date[sim$Imports>0])
  sim$Imports[sim$Imports==0] <- NA
  sim$approach <- factor(sim$approach, levels =  c("liberal","conservative"))
  
  approach <- NULL
  imports_hist_plot <- ggplot()+
    geom_col(sim, mapping=aes(x=date, y=Imports, fill=approach,color=approach),position="identity")+
    geom_col(sim[sim$approach=="conservative",], mapping=aes(x=date, y=Imports, fill=approach,color=approach),position="identity")+
    scale_x_date(date_breaks = "1 month", date_labels = "%b", expand = c(0, 0), limits = c(time_window[1]-1,time_window_imports+1))+
    coord_cartesian(xlim=c(time_window[1]-1,time_window_imports+1))+
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
    geom_line(data= sim, aes(x = date, y = proportion, color=approach))+
    scale_color_manual(values= col_9[c(2,4)],name=variant_name, labels= c("liberal","conservative")) +
    theme_minimal()+
    theme(legend.position = "none")+
    scale_y_continuous(limits = c(0,1))+
    theme(plot.subtitle = element_text(hjust = 0.5),
          axis.title.y = element_text(size = 10),
          title = element_text(size = 10))+
    scale_x_date(date_breaks = "1 month", 
                 date_labels = "%b", limits = c(time_window[1]-1,time_window[2]+1))+
    labs(tag="",subtitle = bquote(), x = "", y =paste0(titel_variant,"\n"))
  
  
  cases_day_plot_Re_t <-  ggplot(data= sim)+
    geom_line(data = cov_ch, aes(x = date, y = weigthed_cases),color=col_9[9], size=1)+
    geom_line(data= sim, aes(x = date, y = T_total, group=approach), color="black", alpha =0.5)+
    geom_line(data= sim, aes(x = date, y = T2, color=approach))+
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
  
  if(m %in% c(1)){
    plot_model_ch_alpha <- ggarrange(imports_hist_plot, cases_day_plot_Re_t,proportion_plot_Re_t, 
                                     labels = c("A", "C", "E"),
                                     ncol = 1, nrow = 3)
    simu_alpha <-   sim
    
  }
  if(m %in% c(2)){
    plot_model_ch_delta <- ggarrange(imports_hist_plot, cases_day_plot_Re_t,proportion_plot_Re_t, 
                                     labels = c("B", "D", "F"),
                                     ncol = 1, nrow = 3)
    simu_delta <- sim
    
  }
}
plot_model_ch <- ggarrange(plot_model_ch_alpha, plot_model_ch_delta,
                           ncol = 2, nrow = 1)

ggsave(plot_model_ch, filename = paste0("./data/figures/Figure4", ".pdf"), height = 7, width = 7,  bg = "transparent")
ggsave(plot_model_ch, filename = paste0("./data/figures/Figure4", ".png"), height = 7, width = 7,  bg = "transparent")
ggsave(ggarrange(cluster_size_plot_conserv_alpha,cluster_size_plot_liberal_alpha, cluster_size_plot_conserv_delta, cluster_size_plot_liberal_delta,ncol = 2, nrow = 2), filename = paste0("./data/figures/SF5.png"), height = 5, width = 7,  bg = "transparent")
ggsave(ggarrange(cluster_size_plot_conserv_alpha,cluster_size_plot_liberal_alpha, cluster_size_plot_conserv_delta, cluster_size_plot_liberal_delta,ncol = 2, nrow = 2), filename = paste0("./data/figures/SF5.pdf"), height = 5, width = 7,  bg = "transparent")
}


