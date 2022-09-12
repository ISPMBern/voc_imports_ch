#' ---
#' title: Imports and the Swiss SARS-CoV-2 Epidemic
#' subtitle: model spread of new variants in Switzerland
#' author: "Martina Reichmuth"
#' date: "'01/12/2021"
#' ---

# grep time were 50% is reached 1) closing borders, 2) strength of surveillance
scenarios_closure <- as.data.frame(matrix(nrow=(length(seq(as_date(min(imports_alpha$MinDateSwissChildren)),as_date(max(imports_alpha$MinDateSwissChildren)),1))+4+
                                                  length(seq(as_date(min(imports_delta$MinDateSwissChildren)),as_date(max(imports_delta$MinDateSwissChildren)),1))), ncol=3))
colnames(scenarios_closure) <-  c("variant", "time_50", "closed_borders")

# Alert about Alpha: 16 Dec 2020 (https://www.bmj.com/content/371/bmj.m4857) widespread on 23 Dec (https://www.theguardian.com/world/2020/dec/23/what-do-we-know-about-the-two-new-covid-19-variants-in-the-uk); https://www.ecdc.europa.eu/en/publications-data/covid-19-risk-assessment-spread-new-sars-cov-2-variants-eueea
concern_alpha <- as_date("2020-12-16")
# Alert about Delta: 24 May (https://www.ecdc.europa.eu/en/publications-data/threat-assessment-emergence-sars-cov-2-b1617-variants)
concern_delta <- as_date("2021-05-24")
scenarios_closure_afterconcern <- as.data.frame(matrix(nrow=(length(seq(concern_alpha,as_date(max(imports_alpha$MinDateSwissChildren)),1))+
                                                               length(seq(concern_delta,as_date(max(imports_delta$MinDateSwissChildren)),1))), ncol=3))
colnames(scenarios_closure_afterconcern) <-  c("variant", "time_50", "closed_borders")

scenarios_surveillance <- as.data.frame(matrix(nrow=(length(imports_alpha$MinDateSwissChildren)+length(imports_delta$MinDateSwissChildren)), ncol=3))
colnames(scenarios_surveillance) <-  c("variant", "time_50", "surveillance")

probs = c(.025,.5,.975)
dates_surveillances<- c()
for (n in 1:2) {
  if(n %in% c(1)){
    # Alpha imports from Emma Hodcroft
    import <- imports_alpha
    import$approach <- ifelse(import$HeadOfClade!="True" & is.na(import$Clade)| import$HeadOfClade=="True" & !is.na(import$Clade)  & import$Clade>-1,"conservative", "liberal" )
    variant_data <- alpha_ch
    variant <- "alpha"
    concern_voc <- concern_alpha 
  
  }
  if(n %in% c(2)){
    # Delta imports from Emma Hodcroft
    import <- imports_delta
     import$approach <- ifelse(import$HeadOfClade!="True" & is.na(import$Clade)| import$HeadOfClade=="True" & !is.na(import$Clade)  & import$Clade>-1,"conservative", "liberal" )
    variant_data <- delta_ch
    variant <- "delta"
    concern_voc <- concern_delta
  }
  
  first_import_date <- as_date(min(import$MinDateSwissChildren))-7
  last_import_date <- max(as_date(import$MinDateSwissChildren))-7
  length_imports <- length(first_import_date:last_import_date)

  time_window <- c(min(variant_data$date),max(variant_data$date))
period <- seq(min(variant_data$date),max(variant_data$date),1)
times_length <- length(period)
length_beforeconcern <- length(seq(min(variant_data$date),concern_voc,1))-1
length_concern_tolastimport <- length(seq(concern_voc,last_import_date,1))
print(length_imports)
print(n)
import <- as.numeric(as_date(import$MinDateSwissChildren) - as.numeric(min(variant_data$date)))
import <-  import-7 # delay of getting tested
length_beforeimport <-  min(import)-1

for (j in 0:length_imports) {
  imports <- import
  imports <- hist(imports, 0:times_length, plot = FALSE)$counts
  imports[length_beforeimport:(length_beforeimport+j)] <- 0
  sim <- imports_003_SEIR(imports,variant )
  
if(n ==1){
  o <- length_imports+1
  j <- 1+j
  scenarios_closure$variant[j] <- "Alpha"
  scenarios_closure$time_50[j] <- suppressWarnings(min(sim$date[sim$proportion >= .5]))
  scenarios_closure$closed_borders[j] <- j-1
}
else if(n ==2){
  scenarios_closure$variant[o+j] <- "Delta"
  scenarios_closure$time_50[o+j] <- suppressWarnings(min(sim$date[sim$proportion>= .5]))
  scenarios_closure$closed_borders[o+j] <- j

}
}

for (i in 0:length_concern_tolastimport) {
  imports <- import
  imports <- hist(imports, 0:times_length, plot = FALSE)$counts
  imports[length_beforeconcern:(length_beforeconcern+i)] <- 0
  sim <- imports_003_SEIR(imports,variant)
  if(n ==1){
    i <- 1+i
    l <- length_concern_tolastimport+1
    scenarios_closure_afterconcern$variant[i] <- "Alpha"
    scenarios_closure_afterconcern$time_50[i] <- suppressWarnings(min(sim$date[sim$proportion >= .5]))
    scenarios_closure_afterconcern$closed_borders[i] <- i-1
  }
  else if(n ==2){
    scenarios_closure_afterconcern$variant[l+i] <- "Delta"
    scenarios_closure_afterconcern$time_50[l+i] <- suppressWarnings(min(sim$date[sim$proportion>= .5]))
    scenarios_closure_afterconcern$closed_borders[l+i] <- i
    
  }
}
#surveillance strategies:

imports_fun <- function(x) sample(x,k,replace = FALSE)
for (k in 1:(length(import))) {
  #print(k)
        #dates_surveillance <- c()
        imports_all <- replicate(100, imports_fun(import), simplify=FALSE)
        imports <- c()
        #for(m in 1:100){
        m <-1
        #imports <- hist(imports_all[[m]], 0:times_length, plot = FALSE)$counts
          imports <- as.numeric(table(cut(imports_all[[m]], breaks=0:times_length)))
        #print(m)
        sim <- imports_003_SEIR(imports,variant)
        sim <- sim[sim$date >=first_import_date,]
        
       # date_surveillance <- suppressWarnings(min(sim$date[sim$proportion >= .5]))

       # dates_surveillance <- c(date_surveillance, dates_surveillance)
       # }
        
      if(n ==1){
        q <- length(import)
        scenarios_surveillance$variant[k] <- "Alpha"
        scenarios_surveillance$time_50[k] <- suppressWarnings(min(sim$date[sim$proportion >= .5]))#paste0(quantile(as.numeric(dates_surveillance), probs),collapse=";")
        #scenarios_surveillance$num_sims[k] <- sum(!as.character(dates_surveillance)%in% NA)
        scenarios_surveillance$surveillance[k] <- 1-(k/length(import))
      }
      else if(n ==2){
        scenarios_surveillance$variant[q+k] <- "Delta"
        scenarios_surveillance$time_50[q+k] <- suppressWarnings(min(sim$date[sim$proportion >= .5]))#paste0(quantile(as.numeric(dates_surveillance), probs),collapse=";")
        #scenarios_surveillance$num_sims[q+k] <- sum(!as.character(dates_surveillance)%in% NA)
        scenarios_surveillance$surveillance[q+k] <- 1-(k/length(import))
      }
}

for(p in c(0.25, 0.5, 0.75)){
  k<- round(length(import)*(1-p))
  imports_all <- replicate(100, imports_fun(import), simplify=FALSE)
  imports <- c()
  date_surveillance <- c()
  for(m in 1:100){
  imports <- hist(imports_all[[m]], 0:times_length, plot = FALSE)$counts
  imports <- as.numeric(table(cut(imports_all[[m]], breaks=0:times_length)))
  #print(m)
  sim <- imports_003_SEIR(imports,variant)
  sim <- sim[sim$date >=first_import_date,]
  date_surveillance <- c(suppressWarnings(min(sim$date[sim$proportion >= .5])), date_surveillance)
  }
  dates_surveillance<- as.data.frame(t(quantile(as.numeric(date_surveillance), probs)))
  dates_surveillance$precent <-  p
  dates_surveillance$variant <- variant
  dates_surveillances<-rbind(dates_surveillance,dates_surveillances)

  }
}

scenarios_closure_backup <- scenarios_closure
scenarios_closure <- scenarios_closure[scenarios_closure$time_50!=Inf,]
scenarios_closure$date <- as_date(scenarios_closure$time_50)
scenarios_closure <- scenarios_closure[!is.na(scenarios_closure$time_50),]
scenarios_closure$delay_50[scenarios_closure$variant=="Alpha"] <- (scenarios_closure$time_50[scenarios_closure$variant=="Alpha"] - min(as.numeric(sim_alpha$date[sim_alpha$proportion>=0.5&sim_alpha$approach=="liberal"]))) #- (as.numeric(as_date(min(imports_alpha$MinDateSwissChildren))) - as.numeric(min(alpha_ch$date)))
scenarios_closure$delay_50[scenarios_closure$variant=="Delta"] <- (scenarios_closure$time_50[scenarios_closure$variant=="Delta"] - min(as.numeric(sim_delta$date[sim_delta$proportion>=0.5&sim_delta$approach=="liberal"]))) #- (as.numeric(as_date(min(imports_delta$MinDateSwissChildren))) - as.numeric(min(delta_ch$date)))
#scenarios_closure <- scenarios_closure[(scenarios_closure$closed_borders > (as.numeric(as_date(min(imports_alpha$MinDateSwissChildren)))-7-1 - as.numeric(min(alpha_ch$date))) &scenarios_closure$variant=="Alpha")|(scenarios_closure$closed_borders > (as.numeric(as_date(min(imports_delta$MinDateSwissChildren)))-7-1 - as.numeric(min(delta_ch$date))) &scenarios_closure$variant=="Delta"),]
#scenarios_closure$closed_borders[scenarios_closure$variant=="Alpha"] <- scenarios_closure$closed_borders[scenarios_closure$variant=="Alpha"] - min(scenarios_closure$closed_borders[scenarios_closure$variant=="Alpha"])
#scenarios_closure$closed_borders[scenarios_closure$variant=="Delta"] <- scenarios_closure$closed_borders[scenarios_closure$variant=="Delta"] - min(scenarios_closure$closed_borders[scenarios_closure$variant=="Delta"])

scenarios_closure_afterconcern<- scenarios_closure_afterconcern[scenarios_closure_afterconcern$time_50!=Inf,]
scenarios_closure_afterconcern$date <- as_date(scenarios_closure_afterconcern$time_50)
scenarios_closure_afterconcern <- scenarios_closure_afterconcern[!is.na(scenarios_closure_afterconcern$time_50),]
scenarios_closure_afterconcern$delay_50[scenarios_closure_afterconcern$variant=="Alpha"] <- scenarios_closure_afterconcern$time_50[scenarios_closure_afterconcern$variant=="Alpha"]  - min(as.numeric(sim_alpha$date[sim_alpha$proportion>=0.5&sim_alpha$approach=="liberal"])) 
scenarios_closure_afterconcern$delay_50[scenarios_closure_afterconcern$variant=="Delta"] <- scenarios_closure_afterconcern$time_50[scenarios_closure_afterconcern$variant=="Delta"] - min(as.numeric(sim_delta$date[sim_delta$proportion>=0.5&sim_delta$approach=="liberal"])) 
#scenarios_closure_afterconcern <- scenarios_closure_afterconcern[(scenarios_closure_afterconcern$closed_borders > (as.numeric(as_date(min(imports_alpha$MinDateSwissChildren)))-7-1 - as.numeric(min(alpha_ch$date))) &scenarios_closure_afterconcern$variant=="Alpha")|(scenarios_closure_afterconcern$closed_borders > (as.numeric(as_date(min(imports_delta$MinDateSwissChildren)))-7-1 - as.numeric(min(delta_ch$date))) &scenarios_closure_afterconcern$variant=="Delta"),]
#scenarios_closure_afterconcern$closed_borders[scenarios_closure_afterconcern$variant=="Alpha"] <- scenarios_closure_afterconcern$closed_borders[scenarios_closure_afterconcern$variant=="Alpha"] - min(scenarios_closure_afterconcern$closed_borders[scenarios_closure_afterconcern$variant=="Alpha"])
#scenarios_closure_afterconcern$closed_borders[scenarios_closure_afterconcern$variant=="Delta"] <- scenarios_closure_afterconcern$closed_borders[scenarios_closure_afterconcern$variant=="Delta"] - min(scenarios_closure_afterconcern$closed_borders[scenarios_closure_afterconcern$variant=="Delta"])

scenarios_surveillance_backup <- scenarios_surveillance
scenarios_surveillance <- scenarios_surveillance[scenarios_surveillance$time_50!=Inf,]
scenarios_surveillance$date <- as_date(as.numeric(scenarios_surveillance$time_50))
scenarios_surveillance <- scenarios_surveillance[!is.na(scenarios_surveillance$time_50),]
scenarios_surveillance$delay_50[scenarios_surveillance$variant=="Alpha"] <- as.numeric(scenarios_surveillance$time_50[scenarios_surveillance$variant=="Alpha"])  - min(as.numeric(sim_alpha$date[sim_alpha$proportion>=0.5&sim_alpha$approach=="liberal"]))
scenarios_surveillance$delay_50[scenarios_surveillance$variant=="Delta"] <- as.numeric(scenarios_surveillance$time_50[scenarios_surveillance$variant=="Delta"]) - min(as.numeric(sim_delta$date[sim_delta$proportion>=0.5&sim_delta$approach=="liberal"]))

dates_surveillance_CI <- dates_surveillances
dates_surveillance_CI[grepl("lpha",dates_surveillance_CI$variant),c(1:3)] <- round(dates_surveillance_CI[grepl("lpha",dates_surveillance_CI$variant),c(1:3)] - min(as.numeric(sim_alpha$date[sim_alpha$proportion>=0.5&sim_alpha$approach=="liberal"])))
dates_surveillance_CI[grepl("elta",dates_surveillance_CI$variant),c(1:3)] <- round(dates_surveillance_CI[grepl("elta",dates_surveillance_CI$variant),c(1:3)] - min(as.numeric(sim_delta$date[sim_delta$proportion>=0.5&sim_delta$approach=="liberal"])))



scenarios_closure$variant <- factor(scenarios_closure$variant , levels = c("Alpha","Delta"))
scenarios_closure_afterconcern$variant <- factor(scenarios_closure_afterconcern$variant , levels = c("Alpha","Delta"))
scenarios_surveillance$variant <- factor(scenarios_surveillance$variant , levels = c("Alpha","Delta"))

scenarios_closure <- scenarios_closure[scenarios_closure$closed_borders<=100,]
lm(formula = scenarios_closure$delay_50 ~ scenarios_closure$closed_borders)
closing_borders_plot <- ggplot(data= scenarios_closure)+
  geom_line(aes(x=closed_borders,y=delay_50, color=variant), size=2)+
  scale_color_manual(values= col_variants[c(1,2)],name=" ", labels= c("Alpha","Delta")) +
  theme_minimal()+
  scale_y_continuous(limits=c(0,60))+
  theme(plot.subtitle = element_text(hjust = 0.5),
        #legend.position="none",
        axis.title.y = element_text(size = 10),
        title = element_text(size = 10))+
  labs(tag="",subtitle = bquote(), x = "Closed borders (in days)", y="Delay to dominance (in days)")
ggsave(closing_borders_plot, filename = paste0("./data/figures/closing_borders_overtake_VoC_",format(Sys.time(), "%b"), ".pdf"), height = 2, width = 4,  bg = "transparent")
ggsave(closing_borders_plot, filename = paste0("./data/figures/closing_borders_overtake_VoC_",format(Sys.time(), "%b"), ".png"), height = 2, width = 4,  bg = "transparent")


lm(formula = scenarios_closure_afterconcern$delay_50[scenarios_closure_afterconcern$closed_borders<=60] ~ scenarios_closure_afterconcern$closed_borders[scenarios_closure_afterconcern$closed_borders<=60])
scenarios_closure_afterconcern_plot <- ggplot(data= scenarios_closure_afterconcern[scenarios_closure_afterconcern$closed_borders<=60,])+
  geom_line(aes(x=closed_borders,y=delay_50, color=variant), size=2)+
  scale_color_manual(values= col_variants[c(1,2)],name=" ", labels= c("Alpha","Delta")) +
  theme_minimal()+
  scale_y_continuous(limits=c(0,60))+
  theme(plot.subtitle = element_text(hjust = 0.5),
        legend.position="none",
        axis.title.y = element_text(size = 10),
        title = element_text(size = 10))+
  labs(tag="",subtitle = bquote(), x = "Closed borders after awareness (in days)", y="Delay to dominance (in days)")
ggsave(scenarios_closure_afterconcern_plot, filename = paste0("./data/figures/scenarios_closure_afterconcern_VoC",format(Sys.time(), "%b"), ".pdf"), height = 2, width = 4,  bg = "transparent")
ggsave(scenarios_closure_afterconcern_plot, filename = paste0("./data/figures/scenarios_closure_afterconcern_VoC",format(Sys.time(), "%b"), ".png"), height = 2, width = 4,  bg = "transparent")

scenarios_surveillance$surveillance100<- scenarios_surveillance$surveillance*100
lm(formula = scenarios_surveillance$delay_50 ~ scenarios_surveillance$surveillance100)
change_surveillance_plot$variant <- factor(change_surveillance_plot$variant, levels =  c("Alpha","Delta"))
change_surveillance_plot <- ggplot(data= scenarios_surveillance)+
  geom_line(aes(x=surveillance,y=delay_50, color=variant), size=2)+
  scale_color_manual(values= col_variants[c(1,2)],name=" ", labels= c("Alpha","Delta")) +
  theme_minimal()+
  scale_y_continuous(limits=c(0,60))+
  scale_x_continuous(breaks=c(0,0.25,0.5,0.75,0.99))+#(0,0.25,0.5,0.75,0.99))+
  theme(plot.subtitle = element_text(hjust = 0.5),
        legend.position="none",
        axis.title.y = element_text(size = 10),
        title = element_text(size = 10))+
  labs(tag="",subtitle = bquote(), x = "Level of surveillance (in %)", y="Delay to dominance (in days)")
ggsave(change_surveillance_plot, filename = paste0("./data/figures/change_surveillance_overtake_VoC_",format(Sys.time(), "%b"), ".pdf"), height = 2, width = 4,  bg = "transparent")
ggsave(change_surveillance_plot, filename = paste0("./data/figures/change_surveillance_overtake_VoC_",format(Sys.time(), "%b"), ".png"), height = 2, width = 4,  bg = "transparent")





g_legend <- function(a.gplot,num){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]] 
  if(num==1){
    legend$grobs[[1]]$grobs[[1]] <-  editGrob(legend$grobs[[1]]$grobs[[1]], gp=gpar(fill ="transparent",col="transparent"))
  }
  if(num==2){
    legend$grobs[[1]]$grobs[[1]] <-  editGrob(legend$grobs[[1]]$grobs[[1]], gp=gpar(fill ="transparent",shape="transparent"))
    legend$grobs[[2]]$grobs[[1]] <-  editGrob(legend$grobs[[2]]$grobs[[1]], gp=gpar(fill ="transparent",shape="transparent"))
  }
  legend
} 

legend_voc <- g_legend(closing_borders_plot,1)
closing_borders_plot <- closing_borders_plot+theme(legend.position="none")


plot_model_scenarios <- ggarrange(closing_borders_plot,scenarios_closure_afterconcern_plot, change_surveillance_plot,legend_voc,
                                  labels = c("A","B","C", ""),
                                 ncol = 2, nrow = 2)
ggsave(plot_model_scenarios, filename = paste0("./data/figures/Figure3_counterfactual_",format(Sys.time(), "%Y-%m-%d"), ".pdf"), height = 6, width = 8,  bg = "transparent")
ggsave(plot_model_scenarios, filename = paste0("./data/figures/Figure3_counterfactual_",format(Sys.time(), "%Y-%m-%d"), ".png"), height = 6, width = 8,  bg = "transparent")

