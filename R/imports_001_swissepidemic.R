#' ---
#' title: Imports and the Swiss SARS-CoV-2 Epidemic
#' subtitle: data loading and decriptive part
#' author: "Martina Reichmuth"
#' date: "'01/12/2021"
#' ---




## Swiss surveillance of SARS-CoV-2 meta-data (reported by Federal Office of Public Health):
url <- readLines("https://www.covid19.admin.ch/api/data/context", warn=FALSE)
url <- gsub("^(.*)https:", "https:", url[19])
url <- gsub("\"", "",url)
url <- gsub(",", "",url)
download(url, dest=paste0("covid19_bag_", Sys.Date(),".zip"), mode="wb") 
unzip(paste0("covid19_bag_", Sys.Date(),".zip"), exdir = "./temp_data")

swiss_cov_metadata <- read.csv("./temp_data/data/COVID19Cases_geoRegion.csv")
colnames(swiss_cov_metadata)[2] <- "date"
colnames(swiss_cov_metadata)[3] <- "cases_num"
cantons_ch <- c("CH", "AG","AI","AR","BE","BL","BS","FR","GE",
                "GL","GR", "JU","LU","NE","NW","OW","SG","SH",
                "SO","SZ","TG",  "TI","UR","VD","VS","ZG","ZH")
swiss_cov_metadata <- swiss_cov_metadata[swiss_cov_metadata$geoRegio %in% cantons_ch,]

swiss_cov_tests <- read.csv("./temp_data/data/COVID19Test_geoRegion_PCR_Antigen.csv")
colnames(swiss_cov_tests)[1] <- "date"
swiss_cov_tests <- swiss_cov_tests[swiss_cov_tests$geoRegio %in% cantons_ch,]

swiss_cov_tests_pcr <- swiss_cov_tests[swiss_cov_tests$nachweismethode=="PCR",]
colnames(swiss_cov_tests_pcr)[2] <- "pcrtests_num"
colnames(swiss_cov_tests_pcr)[3] <- "pcrtests_pos_num"

swiss_cov_tests_antig <- swiss_cov_tests[swiss_cov_tests$nachweismethode=="Antigen_Schnelltest",]
colnames(swiss_cov_tests_antig)[2] <- "antigtests_num"
colnames(swiss_cov_tests_antig)[3] <- "antigtests_pos_num"
geoR <- grep("geoRegion",colnames(swiss_cov_tests))
swiss_cov_tests <- merge(swiss_cov_tests_antig[,c(1,2,3,geoR)], swiss_cov_tests_pcr[,c(1,2,3,geoR)],by=c("date","geoRegion"), all=TRUE )
swiss_cov_tests<- swiss_cov_tests %>% 
  rowwise() %>% #rowwise will make sure the sum operation will occur on each row
  mutate(tests_num = sum(antigtests_num,pcrtests_num, na.rm=TRUE))%>% 
  mutate(tests_pos_num = sum(antigtests_pos_num,pcrtests_pos_num, na.rm=TRUE))
swiss_cov_tests[is.na(swiss_cov_tests)] <- 0

swiss_cov_Re <- read.csv("./temp_data/data/COVID19Re_geoRegion.csv")
swiss_cov_Re <- swiss_cov_Re[swiss_cov_Re$geoRegio %in% cantons_ch,]
swiss_cov <- c()
swiss_cov <- merge(swiss_cov_metadata[,c("date","geoRegion","cases_num", "pop")], swiss_cov_tests[,c("date","geoRegion","tests_num", "tests_pos_num")],by=c("date","geoRegion"), all=TRUE )
swiss_cov <- merge(swiss_cov, swiss_cov_Re[,c("date","geoRegion", "median_R_mean", "median_R_highHPD", "median_R_lowHPD")],by=c("date","geoRegion"), all=TRUE )
swiss_cov <- swiss_cov[!is.na(swiss_cov$cases_num),]
swiss_cov$date <- as_date(swiss_cov$date)
swiss_cov <- as.data.frame(swiss_cov)
swiss_cov <- swiss_cov[swiss_cov$geoRegion=="CH",]

#list.files("./temp_data/data/")
#"COVID19FullyVaccPersons_indication_w_v2.csv" 
#"COVID19FullyVaccPersons_vaccine_v2.csv"
swiss_cov_vaccine <- read.csv("./temp_data/data/COVID19FullyVaccPersons_indication_w_v2.csv")
colnames(swiss_cov_vaccine)[4] <- "fulvacc_num"
swiss_cov_vaccine$week <- swiss_cov_vaccine$date
swiss_cov_vaccine <- swiss_cov_vaccine[swiss_cov_vaccine$geoRegio %in% cantons_ch,]
swiss_cov$week <- week(swiss_cov$date)
swiss_cov$week <-sprintf("%02d",swiss_cov$week)
swiss_cov$week <- paste0(year(swiss_cov$date), swiss_cov$week)
swiss_cov_week <- as.data.frame(swiss_cov  %>% group_by(week, geoRegion)%>% 
  summarise(pop = unique(pop),
            tests_num = sum(tests_num),
            tests_pos_num = sum(tests_pos_num),
            median_R_mean = median(median_R_mean),
            median_R_highHPD = median(median_R_highHPD),
            median_R_lowHPD = median(median_R_lowHPD)))

swiss_cov_week <- merge(swiss_cov_week, swiss_cov_vaccine[,c("week","geoRegion", "fulvacc_num")],by=c("week","geoRegion"), all=TRUE )

unlink("temp_data", recursive = TRUE)
unlink(paste0("covid19_bag_", Sys.Date(),".zip"), recursive = TRUE)
remove(swiss_cov_metadata)
remove(swiss_cov_tests)
remove(swiss_cov_Re)
remove(swiss_cov_vaccine)
remove(swiss_cov_tests_antig)
remove(swiss_cov_tests_pcr)

## KOF stringency, gives information about the measures in place
KOF <- read.csv(paste0("https://datenservice.kof.ethz.ch/api/v1/public/sets/stringency_plus_web?mime=csv&df=Y-m-d.csv"))
KOF[,"date"] <- seq(as_date("2020-01-01"),(as_date("2020-01-01")+length(KOF[,"date"])-1),1)
KOF$date <- as_date(KOF$date)

##SEROCoV-POP study: https://www.corona-immunitas.ch/en/program/results/
seroprevalence_ch <- data.frame(date=c("2020-10-06","2021-02-01"),p=c(0.15,0.25))#https://www.corona-immunitas.ch/en/program/results/
seroprevalence_ch$date <- as_date(seroprevalence_ch$date)


## Swiss surveillance of SARS-CoV-2 sequencing meta-data (reported by ETH Zurich):
#url <- GET("https://cov-spectrum.ethz.ch/gisaid/api/v1/sample/aggregated?country=Switzerland&fields=date,division,pangoLineage")#Used since 17 Nov 2021
url <- GET("https://cov-spectrum.ethz.ch/gisaid/api/v1/sample/aggregated?country=Switzerland&fields=date,division,pangoLineage&accessKey=9Cb3CqmrFnVjO3XCxQLO6gUnKPd")#Used since 17 Nov 2021
jsonRespParsed<- content(url,as="parsed", encoding="UTF-8") 
seq_ch <- jsonRespParsed%>%bind_rows#%>%select(date,division,pangoLineage)# %>%subset(.,country %in% "Switzerland") #%>%
seq_ch <- seq_ch[,c("date","division","pangoLineage","count")]
seq_ch <- seq_ch[seq_ch$count>0&!is.na(seq_ch$count),]
remove(url)
remove(jsonRespParsed)
seq_ch$country <- "CH"
seq_ch <- seq_ch[rep(row.names(seq_ch), seq_ch$count), c(1,2,3,5)]
seq_ch <- seq_ch[!is.na(seq_ch$pangoLineage),]



who_variant_names <- function(x){
  if(is.na(x)){return("undetermined")}
  else if(grepl("Alpha|alpha|Q.1|Q.2|Q.3|Q.4|Q.6",x)){return("Alpha")}
  else if(grepl("B.1.1.7",x, fixed  = TRUE)){return("Alpha")}
  else if(grepl("Beta|beta|B.1.351|B.1.351.1|B.1.351.2",x,useBytes = TRUE)){return("Beta")}
  else if(grepl("Gamma|gamma|P.1|P.1.1|P.1.2|P.1.3|P.1.4|P.1.5|P.1.6",x,useBytes = TRUE)){return("Gamma")}
  else if(grepl("Delta|delta|B.1.617.2|AY.1.1|AY.2|AY.3|AY.3.1|AY.4|AY.5|AY.5.1|AY.5.2|AY.6|AY.7|AY.7.1|AY.7.2|AY.8|AY.9|AY.10|AY.11|AY.12|AY.13|AY.14|AY.15|AY.16|AY.17|AY.18|AY.19|AY.20|AY.21|AY.22|AY.23|AY.24|AY.25|AY.26|AY.27|AY.28|AY.29|AY.30|AY.31|AY.32|AY.33|AY.34|AY.35|AY.36|AY.37",x,useBytes = TRUE)){return("Delta")}
  else if(grepl("C.37|Lambda|lambda",x)){return("Lambda")}
  else if(grepl("C.36",x)){return("C.36*")}
  else if(grepl("Mu|mu|B.1.621|B.1.621.1|B.1.621.2|B.1.621.3",x,useBytes = TRUE)){return("Mu")}#useBytes = FALSE
  else if(grepl("B.1.1.318|AZ.2|AZ.",x)){return("B.1.1.318")}#,useBytes = FALSE
  else if(grepl("B.1.1.529",x, useBytes = TRUE)){return("Omicron")}#,useBytes = FALSE
  else{return("others")}
}
#seq_ch$who_variants <- sapply(seq_ch$pangoLineage, who_variant_names)
seq_ch$who_variants <- "others"
seq_ch$who_variants[seq_ch$pangoLineage %in% c("Alpha","alpha","B.1.1.7","Q.1","Q.2","Q.3","Q.4","Q.6")] <- "Alpha"
seq_ch$who_variants[seq_ch$pangoLineage %in% c("Beta","beta","B.1.351","B.1.351.1","B.1.351.2")] <- "Beta"
seq_ch$who_variants[seq_ch$pangoLineage %in% c("Gamma","gamma","P.1","P.1.1","P.1.2","P.1.3","P.1.4","P.1.5","P.1.6")] <- "Gamma"
#seq_ch$who_variants[seq_ch$pangoLineage %in% c("Delta","delta","B.1.617.2","AY.1.1","AY.2","AY.3","AY.3.1","AY.4","AY.5","AY.5.1","AY.5.2","AY.6","AY.7","AY.7.1","AY.7.2","AY.8","AY.9","AY.10","AY.11","AY.12","AY.13","AY.14","AY.15","AY.16","AY.17","AY.18","AY.19","AY.20","AY.21","AY.22","AY.23","AY.24","AY.25","AY.26","AY.27","AY.28","AY.29","AY.30","AY.31","AY.32","AY.33","AY.34","AY.35","AY.36","AY.37")] <- "Delta"
seq_ch$who_variants[grepl(c("Delta|delta|B.1.617.2|AY.1.1|AY.2|AY.3|AY.3.1|AY.4|AY.5|AY.5.1|AY.5.2|AY.6|AY.7|AY.7.1|AY.7.2|AY.8|AY.9|AY.10|AY.11|AY.12|AY.13|AY.14|AY.15|AY.16|AY.17|AY.18|AY.19|AY.20|AY.21|AY.22|AY.23|AY.24|AY.25|AY.26|AY.27|AY.28|AY.29|AY.30|AY.31|AY.32|AY.33|AY.34|AY.35|AY.36|AY.37"),seq_ch$pangoLineage)] <- "Delta"
seq_ch$who_variants[seq_ch$pangoLineage %in% c("C.37","Lambda","lambda")] <- "Lambda"
#seq_ch$who_variants[seq_ch$pangoLineage %in% c("C.36")] <- "C.36*"
seq_ch$who_variants[seq_ch$pangoLineage %in% c("Mu","mu","B.1.621","B.1.621.1","B.1.621.2","B.1.621.3")] <- "Mu"
seq_ch$who_variants[seq_ch$pangoLineage %in% c("B.1.1.318","AZ.2","AZ.")] <- "B.1.1.318"
seq_ch$who_variants[seq_ch$pangoLineage %in% c("B.1.1.529")] <- "Omicron"
lev <- c("Alpha",  "Beta",  "Gamma", "Delta","Lambda","B.1.1.318","Mu", "others")
seq_ch$who_variants <- factor(seq_ch$who_variants, levels = lev)

seq_ch$date <- as_date(seq_ch$date)
seq_ch <- seq_ch[!is.na(seq_ch$date),]
seq_ch <- seq_ch[order(seq_ch$date),]
seq_ch <- seq_ch[!(grepl("Delta", seq_ch$who_variants)&grepl(as_date("2021-02-13"), seq_ch$date)),]#misclassified Delta


seq_ch_date  <- seq_ch %>% group_by(date) %>% summarise(seq_num= length(na.omit(date)))
swiss_cov <- merge(swiss_cov, seq_ch_date, by="date", all.x = TRUE)
remove(seq_ch_date)
for(i in 1:length(swiss_cov$date)){
swiss_cov$prop_seq_cases[i] <- ifelse(!is.na(swiss_cov$seq_num[i]),swiss_cov$seq_num[i]/swiss_cov$cases_num[i],NA)
}

# weighted cases num
before <- 3
after <- 3
swiss_cov$weigthed_cases <- NA
weighted_casenum_fun <- function(x){
  set <- subset(swiss_cov, swiss_cov$date >= (as_date(x) - before) & swiss_cov$date <= (as_date(x) + after))
  weigthed<- round(mean(set$cases_num))
  return(weigthed)
}
swiss_cov$weigthed_cases[c(4:(length(swiss_cov$date)-4))] <- sapply(swiss_cov$date[c(4:(length(swiss_cov$date)-4))], weighted_casenum_fun)


head(seq_ch$date[seq_ch$who_variants=="Alpha"]) 
seq_ch$date[seq_ch$who_variants=="Alpha" &seq_ch$date== "2020-01-10"] <- as_date("2020-10-01")
seq_ch$date[seq_ch$who_variants=="Alpha" &seq_ch$date== "2020-05-11"] <- as_date("2020-11-05")
head(seq_ch$date[seq_ch$who_variants=="Delta"]) 
head(seq_ch$date[seq_ch$who_variants=="Omicron"]) 


for (n in c("Alpha","Delta")) {
  if(n =="Alpha"){
    variant_ch <- subset(seq_ch, as_date(date) %in% seq(as_date("2020-10-01"),as_date("2021-05-01"),1))
    cov_ch <- subset(swiss_cov, as_date(date) %in% seq(min(variant_ch$date),max(variant_ch$date),1))
    kof_ch <- subset(KOF, as_date(date) %in% seq(min(variant_ch$date),max(variant_ch$date),1))
    time_window <- c(min(variant_ch$date),max(variant_ch$date))
    period <- seq(min(variant_ch$date),max(variant_ch$date),1)
    
    variant_data <- table(variant_ch$date, variant_ch$who_variants=="Alpha")
    labels_variant <- c("Alpha", "Other variants")
  }
  else if(n =="Delta"){
    variant_ch <- subset(seq_ch, as_date(date) %in% seq(as_date("2021-02-01"),as_date("2021-09-01"),1))
    cov_ch <- subset(swiss_cov, as_date(date) %in% seq(min(variant_ch$date),max(variant_ch$date),1))
    kof_ch <- subset(KOF, as_date(date) %in% seq(min(variant_ch$date),max(variant_ch$date),1))
    time_window <- c(min(variant_ch$date),max(variant_ch$date))
    period <- seq(min(variant_ch$date),max(variant_ch$date),1)
    
    variant_data <- table(variant_ch$date, variant_ch$who_variants=="Delta")
    labels_variant <- c("Delta", "Other variants")
  }
variant_data<- cbind(variant_data[,2],variant_data[,1],rowSums(variant_data))
colnames(variant_data) <- c( "variant", "others", "total")
variant_data <- as.data.frame(variant_data)
variant_data$date <- as_date(rownames(variant_data))
lower <- numeric()
upper <- numeric()
proportion <- numeric()
for(i in 1:length(variant_data$variant)) { #[variant_data$country_restrict == country]
  interval <- binom.test(variant_data$variant[i], variant_data$total[i], conf.level = 0.95)
  lower[i] <- interval$conf.int[1]
  upper[i] <- interval$conf.int[2]
  proportion[i] <- interval$estimate
}
variant_data <- cbind(variant_data, proportion=proportion, lower = lower, upper = upper)

if(n =="Alpha"){
  #ggsave(plot_overview_ch, filename = paste0("./figures/Figure1_Alpha_",format(Sys.time(), "%Y-%m-%d"), ".pdf"), height = 9, width = 8,  bg = "transparent")
  alpha_ch <- variant_data
  #cov_ch <- subset(swiss_cov, as_date(date) %in% seq(min(variant_ch$date),max(variant_ch$date),1))
  #kof_ch <- subset(KOF, as_date(date) %in% seq(min(variant_ch$date),max(variant_ch$date),1))
  #time_window <- c(min(variant_ch$date),max(variant_ch$date))
  #period <- seq(min(variant_ch$date),max(variant_ch$date),1)

}
else if(n =="Delta"){
  #ggsave(plot_overview_ch, filename = paste0("./figures/Figure1_Delta_",format(Sys.time(), "%Y-%m-%d"), ".pdf"), height = 9, width = 8,  bg = "transparent")
  delta_ch <- variant_data
  #cov_ch <- subset(swiss_cov, as_date(date) %in% seq(min(variant_ch$date),max(variant_ch$date),1))
  #kof_ch <- subset(KOF, as_date(date) %in% seq(min(variant_ch$date),max(variant_ch$date),1))
  #time_window <- c(min(variant_ch$date),max(variant_ch$date))
  #period <- seq(min(variant_ch$date),max(variant_ch$date),1)
  
}
}
write.csv(swiss_cov, "swiss_cov.csv")

cov_ch <- subset(swiss_cov, as_date(date) %in% seq(min(alpha_ch$date),max(delta_ch$date),1))
variant_data <- subset(seq_ch, as_date(date) %in% seq(min(alpha_ch$date),max(delta_ch$date),1))
variants_melt <- variant_data %>%  group_by(date) %>% 
  summarise(num_seq = n())%>% ungroup()
variants_melt <- variant_data %>%  group_by(date, who_variants) %>% 
  summarise(num_seq = n())%>% ungroup()
variants_melt$who_variants <- factor(variants_melt$who_variants, levels = lev)

kof_ch <- subset(KOF, as_date(date) %in% seq(min(variant_data$date),max(variant_data$date),1))
time_window <- c(min(variant_data$date),max(variant_data$date))
time_window_plots <- as_date(paste0(gsub('.{2}$', '', time_window),"01"))
period <- seq(min(variant_data$date),max(variant_data$date),1)


# Figure 1 (overview of the Swiss epidemic)
col_9 <- (brewer.pal(9,"Set1"))
#Alpha Beta Gamma Delta Lambda B.1.1.318 C.36* Mu, others # Omicron= #9c6fbf
#col_variants <- c("#690c0c", "#c94a36", "#f28fa1","#305c23","#5e294e", "#999945","#b1cc8b", "#5e5e5d")
col_variants <- c("#fc8d62", "#8da0cb", "#5e5e5d")
col_alpha <- c("#fcbda4", "#fc8d62", "#94533a")
col_delta <- c("#a9b5cf", "#8da0cb", "#5f6b87")
variants_melt$who_alpha_delta <- variants_melt$who_variants
variants_melt$who_alpha_delta[!variants_melt$who_variants %in% c("Alpha", "Delta")] <- "others"
variants_melt$who_alpha_delta <- factor(variants_melt$who_alpha_delta, levels = c("Alpha", "Delta","others"))


yscaling <- function(l) {
  l <- format(l, scientific = TRUE)
  l <- gsub("^(.*)e", "e", l)
  l <- gsub("\\+", "", l)
  l <- gsub("e", "10^", l)
  parse(text=l)
}
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

cov_ch_plot <- ggplot(data = cov_ch, aes(x = date, y = weigthed_cases,fill=col_9[9],color=col_9[9]))+ 
  geom_density(stat = 'identity')+
  scale_x_date(date_breaks = "1 month", 
               date_labels = "%b", limits = c(time_window_plots[1],time_window_plots[2]))+
  scale_color_manual(values= col_9[9],label=c("", ""), name="") +
  scale_fill_manual(values= col_9[9],label=c("", ""), name="") +
  theme_minimal()+
  theme(legend.position = "none")+
  #coord_cartesian(expand = FALSE)+
  labs(tag=bquote(.("")),subtitle = "", x = "", y =bquote("Number of reported \nSARS-CoV-2 cases"))
ggsave(cov_ch_plot, filename = paste0("./data/figures/cases_num_",format(Sys.time(), "%Y-%m-%d"), ".png"), height = 2, width = 5,  bg = "transparent")

plot_serop <- ggplot(seroprevalence_ch) +
  geom_point(aes(x=date,
                 y=p),
             #ymin=qbeta(.025,k,n-k),
             #ymax=qbeta(.975,k,n-k)),
             #col=col_9[c(9,9,9,9,9,7)]) +
             col=col_9[9]) +
  theme_minimal()+
  scale_y_continuous(limits = c(0,1))+
  scale_x_date(date_breaks = "1 month", 
               date_labels = "%b", limits = c(time_window_plots[1],time_window_plots[2]))+
  labs(tag="",x = "", y =bquote("Seroprevalence\n "))

#key data 
#NPI
#https://de.wikipedia.org/wiki/COVID-19-Pandemie_in_der_Schweiz (Access: 2021-05-04)
#https://www.bag.admin.ch/bag/en/home/krankheiten/ausbrueche-epidemien-pandemien/aktuelle-ausbrueche-epidemien/novel-cov/massnahmen-des-bundes.html (Access: 2021-07-19)
kof_plot <- ggplot(data= kof_ch)+
  theme_minimal()+
  geom_line(data= kof_ch,aes(x= date, y=ch.kof.stringency.ch.stringency_plus), col=col_9[9],size = 2) +
  scale_x_date(date_breaks = "1 month", 
               date_labels = "%b",
               limits = as_date(c(time_window_plots[1],time_window_plots[2])))+
 
  annotate("pointrange", x = as_date("2020-10-19"), y = kof_ch$ch.kof.stringency.ch.stringency_plus[kof_ch$date==as_date("2020-10-19")], ymin = kof_ch$ch.kof.stringency.ch.stringency_plus[kof_ch$date==as_date("2020-10-19")], ymax = kof_ch$ch.kof.stringency.ch.stringency_plus[kof_ch$date==as_date("2020-10-19")]+17,colour = col_9[9], size = 0.8)+
  annotate("text",hjust =0.5,vjust =0, x=as_date("2020-10-19"), y=kof_ch$ch.kof.stringency.ch.stringency_plus[kof_ch$date==as_date("2020-10-19")]+20, label= "max 15 people",size=1.8,colour = "black") + 
  
  annotate("pointrange", x = as_date("2020-10-29"), y = kof_ch$ch.kof.stringency.ch.stringency_plus[kof_ch$date==as_date("2020-10-29")], ymin = kof_ch$ch.kof.stringency.ch.stringency_plus[kof_ch$date==as_date("2020-10-29")], ymax = kof_ch$ch.kof.stringency.ch.stringency_plus[kof_ch$date==as_date("2020-10-29")]-42,colour = col_9[9], size = 0.8)+
  annotate("text",hjust =0.5,vjust =1, x=as_date("2020-10-29"), y=kof_ch$ch.kof.stringency.ch.stringency_plus[kof_ch$date==as_date("2020-10-29")]-45, label= "max 10 people if private \notherwise 15 with mask",size=1.8,colour = "black") + 
  
  annotate("pointrange", x = as_date("2020-11-02"), y = kof_ch$ch.kof.stringency.ch.stringency_plus[kof_ch$date==as_date("2020-11-02")], ymin = kof_ch$ch.kof.stringency.ch.stringency_plus[kof_ch$date==as_date("2020-11-02")], ymax = kof_ch$ch.kof.stringency.ch.stringency_plus[kof_ch$date==as_date("2020-11-02")]+32,colour = col_9[9], size = 0.8)+
  annotate("text",hjust =0.5,vjust =0, x=as_date("2020-11-02"), y=kof_ch$ch.kof.stringency.ch.stringency_plus[kof_ch$date==as_date("2020-11-02")]+35, label= "Online learning \nat higher educational institutes",size=1.8,colour = "black") + 
  
  annotate("pointrange", x = as_date("2020-12-12"), y = kof_ch$ch.kof.stringency.ch.stringency_plus[kof_ch$date==as_date("2020-12-12")], ymin = kof_ch$ch.kof.stringency.ch.stringency_plus[kof_ch$date==as_date("2020-12-12")], ymax = kof_ch$ch.kof.stringency.ch.stringency_plus[kof_ch$date==as_date("2020-12-12")]+12,colour = col_9[9], size = 0.8)+
  annotate("text",hjust =0.5,vjust =0, x=as_date("2020-12-12"), y=kof_ch$ch.kof.stringency.ch.stringency_plus[kof_ch$date==as_date("2020-12-12")]+15, label= "Sperrstunde",size=1.8,colour = "black") + 
  
  annotate("pointrange", x = as_date("2020-12-21"), y = kof_ch$ch.kof.stringency.ch.stringency_plus[kof_ch$date==as_date("2020-12-21")], ymin = kof_ch$ch.kof.stringency.ch.stringency_plus[kof_ch$date==as_date("2020-12-21")], ymax = kof_ch$ch.kof.stringency.ch.stringency_plus[kof_ch$date==as_date("2020-12-21")]-32,colour = col_9[9], size = 0.8)+
  annotate("text",hjust =0.5,vjust =1, x=as_date("2020-12-21"), y=kof_ch$ch.kof.stringency.ch.stringency_plus[kof_ch$date==as_date("2020-12-21")]-35, label= "quarantine restrictions \non countries with Alpha",size=1.8,colour = "black") + 
  
  annotate("pointrange", x = as_date("2021-01-18"), y = kof_ch$ch.kof.stringency.ch.stringency_plus[kof_ch$date==as_date("2021-01-18")], ymin = kof_ch$ch.kof.stringency.ch.stringency_plus[kof_ch$date==as_date("2021-01-18")], ymax = kof_ch$ch.kof.stringency.ch.stringency_plus[kof_ch$date==as_date("2021-01-18")]+17,colour = col_9[9], size = 0.8)+
  annotate("text",hjust =0.5,vjust =0, x=as_date("2021-01-18"), y=kof_ch$ch.kof.stringency.ch.stringency_plus[kof_ch$date==as_date("2021-01-18")]+20, label= "mandatory home-office, \nmax 5 people, closures, etc.",size=1.8,colour = "black") + 
  
  annotate("pointrange", x = as_date("2021-01-28"), y = kof_ch$ch.kof.stringency.ch.stringency_plus[kof_ch$date==as_date("2021-01-28")], ymin = kof_ch$ch.kof.stringency.ch.stringency_plus[kof_ch$date==as_date("2021-01-28")], ymax = kof_ch$ch.kof.stringency.ch.stringency_plus[kof_ch$date==as_date("2021-01-28")]-7,colour = col_9[9], size = 0.8)+
  annotate("text",hjust =0.5,vjust =1, x=as_date("2021-01-28"), y=kof_ch$ch.kof.stringency.ch.stringency_plus[kof_ch$date==as_date("2021-01-28")]-10, label= "free testing",size=1.8,colour = "black") + 
  
  annotate("pointrange", x = as_date("2021-03-01"), y = kof_ch$ch.kof.stringency.ch.stringency_plus[kof_ch$date==as_date("2021-03-01")], ymin = kof_ch$ch.kof.stringency.ch.stringency_plus[kof_ch$date==as_date("2021-03-01")], ymax = kof_ch$ch.kof.stringency.ch.stringency_plus[kof_ch$date==as_date("2021-03-01")]-12,colour = col_9[9], size = 0.8)+
  annotate("text",hjust =0.5,vjust =1, x=as_date("2021-03-01"), y=kof_ch$ch.kof.stringency.ch.stringency_plus[kof_ch$date==as_date("2021-03-01")]-15, label= "max of 15 people outside, \nre-opening of all shops, \nleisure and sport activities \nbut mandatory mask wearing",size=1.8,colour = "black") + 
  
  annotate("pointrange", x = as_date("2021-04-07"), y = kof_ch$ch.kof.stringency.ch.stringency_plus[kof_ch$date==as_date("2021-04-07")], ymin = kof_ch$ch.kof.stringency.ch.stringency_plus[kof_ch$date==as_date("2021-04-07")], ymax = kof_ch$ch.kof.stringency.ch.stringency_plus[kof_ch$date==as_date("2021-04-07")]+27,colour = col_9[9], size = 0.8)+
  annotate("text",hjust =0.5,vjust =0, x=as_date("2021-04-07"), y=kof_ch$ch.kof.stringency.ch.stringency_plus[kof_ch$date==as_date("2021-04-07")]+30, label= "free self-tests",size=1.8,colour = "black") + 
  
  annotate("pointrange", x = as_date("2021-04-19"), y = kof_ch$ch.kof.stringency.ch.stringency_plus[kof_ch$date==as_date("2021-04-19")], ymin = kof_ch$ch.kof.stringency.ch.stringency_plus[kof_ch$date==as_date("2021-04-19")], ymax = kof_ch$ch.kof.stringency.ch.stringency_plus[kof_ch$date==as_date("2021-04-19")]-32,colour = col_9[9], size = 0.8)+
  annotate("text",hjust =0.5,vjust =1, x=as_date("2021-04-19"), y=kof_ch$ch.kof.stringency.ch.stringency_plus[kof_ch$date==as_date("2021-04-19")]-35, label= "further re-openings, \ne.g., restaurants outdoors \nbut with capacity limits",size=1.8,colour = "black") + 
  
  annotate("pointrange", x = as_date("2021-05-31"), y = kof_ch$ch.kof.stringency.ch.stringency_plus[kof_ch$date==as_date("2021-05-31")], ymin = kof_ch$ch.kof.stringency.ch.stringency_plus[kof_ch$date==as_date("2021-05-31")], ymax = kof_ch$ch.kof.stringency.ch.stringency_plus[kof_ch$date==as_date("2021-05-31")]+37,colour = col_9[9], size = 0.8)+
  annotate("text",hjust =0.5,vjust =0, x=as_date("2021-05-31"), y=kof_ch$ch.kof.stringency.ch.stringency_plus[kof_ch$date==as_date("2021-05-31")]+40, label= "further re-openings, \ne.g., restaurants opened indoors",size=1.8,colour = "black") + 
  
  annotate("pointrange", x = as_date("2021-06-26"), y = kof_ch$ch.kof.stringency.ch.stringency_plus[kof_ch$date==as_date("2021-06-26")], ymin = kof_ch$ch.kof.stringency.ch.stringency_plus[kof_ch$date==as_date("2021-06-26")], ymax = kof_ch$ch.kof.stringency.ch.stringency_plus[kof_ch$date==as_date("2021-06-26")]-22,colour = col_9[9], size = 0.8)+
  annotate("text",hjust =0.5,vjust =1, x=as_date("2021-06-26"), y=kof_ch$ch.kof.stringency.ch.stringency_plus[kof_ch$date==as_date("2021-06-26")]-25, label= "wide re-opening, \nbut with COVID certificate",size=1.8,colour = "black") + 
  
  scale_y_continuous(limits = c(0,100))+
  labs(tag="",x = "", y =bquote("KOF Stringency-\nPlus Index"))


re_plot <- ggplot()+
  geom_hline(aes(yintercept=1),
             linetype="dashed", size=1, colour="black")+
  theme_minimal()+
  geom_ribbon(data= cov_ch, aes(x = date, y = median_R_mean, ymin = median_R_lowHPD,ymax = median_R_highHPD),fill= col_9[9],alpha=0.4)+
  scale_x_date(date_breaks = "1 month", 
               date_labels = "%b", limits = c(time_window_plots[1],time_window_plots[2]))+
  labs(tag="",x = "", y =bquote(atop("Estimated"~italic("R"["e"]))))#," \n (by ETH Zurich)"


proportion_sequenced_plot <- ggplot()+
  geom_bar(data = cov_ch, aes(x = date, y = prop_seq_cases), fill = col_9[9], color = col_9[9],stat = 'identity', position = 'stack',width=1)+
  scale_x_date(date_breaks = "1 month", 
               date_labels = "%b",
               limits = as_date(c(time_window_plots[1],time_window_plots[2])))+
  theme_minimal()+
  #coord_cartesian(expand = FALSE)+
  labs(tag=bquote(.("")),subtitle = "", x = "", y =bquote("Proportion sequenced \nof reported SARS-CoV-2"))
ggsave(proportion_sequenced_plot, filename = paste0("./data/figures/sequences_prob_",format(Sys.time(), "%Y-%m-%d"), ".png"), height = 2, width = 5,  bg = "transparent")



seq_variants_plot <- ggplot()+
  geom_bar(data = variants_melt, aes(x = date, y = num_seq, color = who_alpha_delta,fill = who_alpha_delta),stat = 'identity', position = 'stack')+
  scale_x_date(date_breaks = "1 month", 
               date_labels = "%b",
               limits = as_date(c(time_window_plots[1],time_window_plots[2])))+
  theme_minimal()+
theme(legend.position=c(.35,.8),legend.direction="horizontal")+
  scale_y_continuous(limits = c(0,max(variants_melt$num_seq)))+
  scale_fill_manual(values= col_variants,name="") +#SARS-CoV-2 variants
  scale_color_manual(values= col_variants,name="") +#SARS-CoV-2 variants
  #coord_cartesian(expand = FALSE)+
  labs(tag=bquote(.("")),subtitle = "", x = "", y =bquote("Number of sequenced \nSARS-CoV-2"))
ggsave(seq_variants_plot, filename = paste0("./data/figures/sequences_num_col_",format(Sys.time(), "%Y-%m-%d"), ".png"), height = 2, width = 5,  bg = "transparent")

# 
for(i in variants_melt$date){
  variants_melt$frac_seq[variants_melt$date %in% i] <- variants_melt$num_seq[variants_melt$date %in% i]/cov_ch$seq_num[cov_ch$date %in% i]

}



seq_variants_frac_plot <- ggplot()+
  geom_bar(data = variants_melt, aes(x = date, y = frac_seq, color = who_alpha_delta,fill = who_alpha_delta),stat = 'identity', position = 'stack')+
  scale_x_date(date_breaks = "1 month", 
               date_labels = "%b",
               limits = as_date(c(time_window_plots[1],time_window_plots[2])))+
  theme_minimal()+
  #theme(legend.position=c(0.35,1.2),legend.direction="horizontal")+
  scale_y_continuous(limits = c(0,1))+
  scale_fill_manual(values= col_variants,name="") +#SARS-CoV-2 variants
  scale_color_manual(values= col_variants,name="") +#SARS-CoV-2 variants
  #coord_cartesian(expand = FALSE)+
  labs(tag=bquote(.("")),subtitle = "", x = "", y =bquote("Fraction of sequenced \n SARS-CoV-2 variants"))
ggsave(seq_variants_frac_plot, filename = paste0("./data/figures/sequences_fraction_",format(Sys.time(), "%Y-%m-%d"), ".png"), height = 3, width = 5,  bg = "transparent")



### binomial confidence intervals, will be showed by week
variant_data$date <- as_date(variant_data$date)
variant_data$date_num <- as.numeric(variant_data$date)
variant_data$week <- week(variant_data$date)
variant_data$year <- year(variant_data$date)

variant_week_ch <- variant_data %>% group_by(week, year, who_variants)  %>% summarise(variants= length(who_variants))
seq_week_ch <- variant_week_ch %>% group_by(week, year)  %>% summarise(total= sum(variants))
variant_week <- merge(variant_week_ch, seq_week_ch,by=c("week", "year"))
variant_week$week <-  ifelse(variant_week$week ==53 &variant_week$year ==2020, 1,variant_week$week)
variant_week$year <-ifelse(variant_week$week ==1 &variant_week$year ==2020, 2021,variant_week$year)
variant_week$region <- "CH"
remove(seq_week_ch)
remove(variant_week_ch)

variant_week$total <- as.numeric(variant_week$total)
variant_week$variants <- as.numeric(variant_week$variants)
variant_week <- variant_week[variant_week$total!=0,]
variant_week$who_variants <- factor(variant_week$who_variants, levels = lev)

lower <- numeric()
upper <- numeric()
conf <- numeric()
sampled_ci <- function(data){
  for(i in 1:length(data$who_variants)) {
    int <- binom.test(data$variants[i], data$total[i])#int <- binom.test(data$value[i], data$sequences[i])
    lower[i] <- int$conf.int[1]
    upper[i] <- int$conf.int[2]
    conf[i] <- int$estimate
  }
  data <- cbind(data, conf=conf, lower = lower, upper = upper)
  return(data)
}
variant_week <- sampled_ci(variant_week)
variant_week <- as.data.frame(variant_week)
#time_window_plots <- as_date(paste0(gsub('.{2}$', '', time_window),"01"))
for(i in 1:length(variant_week$year)){
  variant_week$week_day[i] <- as_date(as.Date(paste(variant_week$year[i], variant_week$week[i], 1, sep="-"), "%Y-%U-%u"))
}
variant_week$week_day <- as_date(variant_week$week_day)
remove(conf)
remove(upper)
remove(lower)


### prepare data for multinomial function and plot
## reorganize dataset
variants_ch <- variant_data
## overall Swiss data
variants_ch$who_variants <- factor(variants_ch$who_variants, levels = lev)
### CH variants multinomial models:
mnom_date_spline <- multinom(who_variants ~ ns(date_num, df=2), data = variants_ch)
## predict for mnom_week model
predict_date <- data.frame(date_num = seq((min(variants_ch$date_num)),as.numeric(time_window_plots[2]),1))
predict.eff_date <- Effect("date_num", mnom_date_spline,level=0.95,
                           xlevels=list(date_num=seq(min(variants_ch$date_num),as.numeric(time_window_plots[2]),1)))
predict.eff_date <- cbind(predict_date, melt(predict.eff_date$prob), melt(predict.eff_date$lower.prob),melt(predict.eff_date$upper.prob))
predict.eff_date <- predict.eff_date[,-c(2,5,6,8,9)]
colnames(predict.eff_date) <- c("date_num","who_variants", "prob","lower", "upper")
predict.eff_date$who_variants <- gsub("prob.", "", predict.eff_date$who_variants)

predict.eff_date$who_variants <- factor(predict.eff_date$who_variants, levels = lev)
predict.eff_date$date <- as_date(predict.eff_date$date_num)
remove(predict_date)

## overall Swiss multi-nominal figure
variants_plot_model <- ggplot() + 
  geom_line(data=predict.eff_date[predict.eff_date$who_variants %in% c("Alpha","Delta","others"),], aes(x = date, y = prob, color = who_variants))+
  geom_ribbon(data=predict.eff_date[predict.eff_date$who_variants %in% c("Alpha","Delta","others"),], aes(x = date, y = prob, ymin = lower,ymax = upper,fill=who_variants),alpha=0.4)+
  geom_errorbar(data= variant_week[na.omit(variant_week$region=="CH") & variant_week$who_variants %in% c("Alpha","Delta","others"),], aes(x = week_day, ymin=lower, ymax=upper, color = who_variants), width=.1) +
  geom_point(data= variant_week[na.omit(variant_week$region=="CH") & variant_week$who_variants %in% c("Alpha","Delta","others"),],aes(x = week_day, y=conf, color = who_variants))+
  geom_rect(aes(xmin = as_date(max(variant_data$date)), ymin = 0, xmax = time_window_plots[2], ymax = 1), fill= "#e8e8e8", colour= "transparent", alpha=0.4)+
  scale_x_date(date_breaks = "1 month", 
               date_labels = "%b",
               limits = as_date(c(time_window_plots[1],time_window_plots[2])))+
  scale_color_manual(values= col_variants,name="SARS-CoV-2 variants") +
  scale_fill_manual(values= col_variants,name="SARS-CoV-2 variants") +
  theme_minimal()+
  theme(legend.position = "none")+
  theme(plot.subtitle = element_text(hjust = 0.5),
        axis.title.y = element_text(size = 10),
        title = element_text(size = 10))+
  scale_y_continuous(limits = c(0,1))+
  #coord_cartesian(expand = FALSE)+
  labs(tag=bquote(.("")),subtitle = "", x = "", y =bquote("Proportion of \nSARS-CoV-2 variants\n"))
ggsave(variants_plot_model, filename = paste0("./data/figures/variants_multinomial_model",format(Sys.time(), "%b"), ".png"), height = 2, width = 5,  bg = "transparent")


variants_plot_model_log <- variants_plot_model+
  coord_cartesian(ylim = c(10^-3, 10^0))+
  theme(legend.position = "none")+
  scale_y_continuous(trans='log10',labels=yscaling,limits = c(min(predict.eff_date$lower),1))+
  labs(y =bquote("Proportion of SARS-CoV-2 variants \n(log10-scale)"))


plot_overview_ch <- NULL
plot_overview_ch <- ggarrange(cov_ch_plot, re_plot,kof_plot, proportion_sequenced_plot,seq_variants_plot,
                              variants_plot_model,
                              labels = c("A", "B", "C","D", "E", "F"),
                              ncol = 2, nrow = 3)
ggsave(plot_overview_ch, filename = paste0("./data/figures/Figure1_overview_",format(Sys.time(), "%Y-%m-%d"), ".pdf"), height = 8, width = 8,  bg = "transparent")
ggsave(plot_overview_ch, filename = paste0("./data/figures/Figure1_overview_",format(Sys.time(), "%Y-%m-%d"), ".png"), height = 8, width = 8,  bg = "transparent")





sum(cov_ch$cases_num)
length(seq_ch$date)
sum(seq_ch$who_variants %in% "Alpha")

sum(seq_ch$who_variants%in% "Delta")

