#' ---
#' title: Imports and the Swiss SARS-CoV-2 Epidemic
#' subtitle: data loading and decriptive part
#' author: "Martina Reichmuth"
#' date: "'01/12/2021"
#' ---




## Swiss surveillance of SARS-CoV-2 meta-data (reported by Federal Office of Public Health):
url <- readLines("https://www.covid19.admin.ch/api/data/context", warn=FALSE)
url <- gsub("^(.*)https:", "https:", url[20])
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

unlink("temp_data", recursive = TRUE)
unlink(paste0("covid19_bag_", Sys.Date(),".zip"), recursive = TRUE)
remove(swiss_cov_metadata)
remove(swiss_cov_tests)
remove(swiss_cov_Re)
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
#url <- GET("https://cov-spectrum.ethz.ch/gisaid/api/v1/sample/aggregated?country=Switzerland&fields=date,division,pangoLineage&accessKey=9Cb3CqmrFnVjO3XCxQLO6gUnKPd")#Used since 17 Nov 2021
url <- GET("https://lapis.cov-spectrum.org/open/v1/sample/aggregated?country=Switzerland&fields=date,division,pangoLineage")#Used since Oct 2022
jsonRespParsed<- content(url,as="parsed", encoding="UTF-8") 
seq_ch <- jsonRespParsed%>%bind_rows#%>%select(date,division,pangolinLineage)# %>%subset(.,country %in% "Switzerland") #%>%
seq_ch <- seq_ch[,c("date","division","pangoLineage","count")]
seq_ch <- seq_ch[seq_ch$count>0&!is.na(seq_ch$count),]
seq_ch$country <- "CH"
seq_ch <- seq_ch[rep(row.names(seq_ch), seq_ch$count), c(1,2,3,5)]

seq_ch$who_variants <- "others"
seq_ch$who_variants[seq_ch$pangoLineage %in% c("Alpha","alpha","B.1.1.7","Q.1","Q.2","Q.3","Q.4","Q.6")] <- "Alpha"
seq_ch$who_variants[grepl(c("Delta|delta|B.1.617.2|AY.1.1|AY.2|AY.3|AY.3.1|AY.4|AY.5|AY.5.1|AY.5.2|AY.6|AY.7|AY.7.1|AY.7.2|AY.8|AY.9|AY.10|AY.11|AY.12|AY.13|AY.14|AY.15|AY.16|AY.17|AY.18|AY.19|AY.20|AY.21|AY.22|AY.23|AY.24|AY.25|AY.26|AY.27|AY.28|AY.29|AY.30|AY.31|AY.32|AY.33|AY.34|AY.35|AY.36|AY.37"),seq_ch$pangoLineage)] <- "Delta"

lev <- c("Alpha", "Delta","others")
seq_ch$who_variants <- factor(seq_ch$who_variants, levels = lev)

seq_ch$date <- as_date(seq_ch$date)
seq_ch <- seq_ch[!is.na(seq_ch$date),]
seq_ch <- seq_ch[order(seq_ch$date),]
seq_ch <- seq_ch[!(grepl("Delta", seq_ch$who_variants)&grepl(as_date("2021-02-13"), seq_ch$date)),]#misclassified Delta

seq_ch$value <- 1
seq_ch_date  <- aggregate(seq_ch["value"], by=seq_ch["date"], sum)
colnames(seq_ch_date)[2] <- "num_seq"
swiss_cov <- merge(swiss_cov, seq_ch_date, by="date", all.x = TRUE)
for(i in 1:length(swiss_cov$date)){
  swiss_cov$prop_seq_cases[i] <- ifelse(!is.na(swiss_cov$num_seq[i]),swiss_cov$num_seq[i]/swiss_cov$cases_num[i],NA)
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


#clean some data (data labeling)
seq_ch$date[seq_ch$who_variants=="Alpha" &seq_ch$date== "2020-01-10"] <- as_date("2020-10-01")
seq_ch$date[seq_ch$who_variants=="Alpha" &seq_ch$date== "2020-05-11"] <- as_date("2020-11-05")
seq_ch$date[seq_ch$who_variants=="Alpha" &seq_ch$date== "2020-03-09"] <- NA
head(seq_ch$date[seq_ch$who_variants=="Delta"]) 
seq_ch$date[seq_ch$who_variants=="Delta" &seq_ch$date== "2020-01-02"] <- NA
seq_ch$date[seq_ch$who_variants=="Delta" &seq_ch$date== "2020-01-02"] <- NA
seq_ch$date[seq_ch$who_variants=="Delta" &seq_ch$date== "2020-01-08"] <- NA

for (n in c("Alpha","Delta")) {
  if(n =="Alpha"){
    variant_ch <- subset(seq_ch, as_date(date) %in% seq(as_date("2020-10-01"),as_date("2021-05-01"),1))
    cov_ch <- subset(swiss_cov, as_date(date) %in% seq(min(variant_ch$date),max(variant_ch$date),1))
    kof_ch <- subset(KOF, as_date(date) %in% seq(min(variant_ch$date),max(variant_ch$date),1))
    time_window <- c(min(variant_ch$date),max(variant_ch$date))
    period <- seq(min(variant_ch$date),max(variant_ch$date),1)
    
    variant_data <- table(variant_ch$date, variant_ch$who_variants=="Alpha")
  }
  else if(n =="Delta"){
    variant_ch <- subset(seq_ch, as_date(date) %in% seq(as_date("2021-02-01"),as_date("2021-09-01"),1))
    cov_ch <- subset(swiss_cov, as_date(date) %in% seq(min(variant_ch$date),max(variant_ch$date),1))
    kof_ch <- subset(KOF, as_date(date) %in% seq(min(variant_ch$date),max(variant_ch$date),1))
    time_window <- c(min(variant_ch$date),max(variant_ch$date))
    period <- seq(min(variant_ch$date),max(variant_ch$date),1)
    
    variant_data <- table(variant_ch$date, variant_ch$who_variants=="Delta")
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
    alpha_ch <- variant_data
  }
  else if(n =="Delta"){
    delta_ch <- variant_data
  }
}
write.csv(swiss_cov, paste0("swiss_cov",Sys.Date(),".csv"))

cov_ch <- subset(swiss_cov, as_date(date) %in% seq(min(alpha_ch$date),max(delta_ch$date),1))
variant_data <- subset(seq_ch, as_date(date) %in% seq(min(alpha_ch$date),max(delta_ch$date),1))

variants_melt  <- aggregate(seq_ch["value"], by=seq_ch[c("date", "who_variants")], sum)
colnames(variants_melt)[3] <- "num_seq_variant"
variants_melt$who_variants <- factor(variants_melt$who_variants, levels = lev)

kof_ch <- subset(KOF, as_date(date) %in% seq(min(variant_data$date),max(variant_data$date),1))
time_window <- c(min(variant_data$date),max(variant_data$date))
time_window_plots <- as_date(paste0(gsub('.{2}$', '', time_window),"01"))
period <- seq(min(variant_data$date),max(variant_data$date),1)


# Figure 3 (and Sup Fig 1) (overview of the Swiss epidemic)
col_9 <- (brewer.pal(9,"Set1"))
col_variants <- c("#fc8d62", "#8da0cb", "#5e5e5d")
col_alpha <- c("#fcbda4", "#fc8d62", "#94533a")
col_delta <- c("#a9b5cf", "#8da0cb", "#5f6b87")

# level data
variants_melt$who_variants <- factor(variants_melt$who_variants, levels = lev)
seq_ch$who_variants <- factor(seq_ch$who_variants, levels = lev)

# color cases due to sequencing (proportion of variants)
variants_extra <- aggregate(seq_ch["value"], by=seq_ch[c("date", "who_variants")], sum)
variants_extra <- merge(swiss_cov[,c("weigthed_cases", "date")], variants_extra, by="date",all=TRUE)
variants_extra <- merge(seq_ch_date, variants_extra, by="date",y.all=TRUE)
variants_extra$who_variants <- factor(variants_extra$who_variants, levels = lev)

variants_extra$num_variant_extrapolated <- variants_extra$weigthed_cases*variants_extra$num_seq_variant/variants_extra$num_seq
cov_ch_plot <- ggplot(data = variants_extra, aes(x = date, y = num_variant_extrapolated, color = who_variants,fill = who_variants))+ 
  geom_bar(stat='identity')+
  scale_y_continuous(limits = c(0,1e4))+
  scale_x_date(date_breaks = "1 month", 
               date_labels = "%b", limits = c(time_window_plots[1],time_window_plots[2]))+
  scale_fill_manual(values= col_variants,name="") +#SARS-CoV-2 variants
  scale_color_manual(values= col_variants,name="") +#SARS-CoV-2 variants
  theme_minimal()+
  theme(legend.position=c(.4,.9),legend.direction="horizontal")+
  labs(tag=bquote(.("")),subtitle = "", x = "", y =bquote("Number of reported \nSARS-CoV-2 cases"))
ggsave(cov_ch_plot, filename = paste0("./data/figures/cases_num_",format(Sys.time(), "%Y-%m-%d"), ".png"), height = 2, width = 5,  bg = "transparent")


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
  labs(tag=bquote(.("")),subtitle = "", x = "", y =bquote("Proportion sequenced \nof reported SARS-CoV-2"))

seq_variants_plot <- ggplot()+
  geom_bar(data = variants_extra, aes(x = date, y = num_seq, color = who_variants,fill = who_variants),stat = 'identity', position = 'stack')+
  scale_x_date(date_breaks = "1 month", 
               date_labels = "%b",
               limits = as_date(c(time_window_plots[1],time_window_plots[2])))+
  theme_minimal()+
  theme(legend.position=c(.35,.8),legend.direction="horizontal")+
  scale_y_continuous(limits = c(0,max(variants_melt$num_seq)))+
  scale_fill_manual(values= col_variants,name="") +#SARS-CoV-2 variants
  scale_color_manual(values= col_variants,name="") +#SARS-CoV-2 variants
  labs(tag=bquote(.("")),subtitle = "", x = "", y =bquote("Number of sequenced \nSARS-CoV-2"))
ggsave(seq_variants_plot, filename = paste0("./data/figures/sequences_num_col_",format(Sys.time(), "%Y-%m-%d"), ".png"), height = 2, width = 5,  bg = "transparent")


plot_overview_ch <- NULL
plot_overview_ch <- ggarrange(kof_plot,proportion_sequenced_plot,seq_variants_plot,
                              labels = c("A", "B", "C"),
                              ncol = 1, nrow = 3)
ggsave(plot_overview_ch, filename = paste0("./data/figures/SF1_",format(Sys.time(), "%Y-%m-%d"), ".png"), height =8, width = 5,  bg = "transparent")

plot_overview_ch <- ggarrange(cov_ch_plot, re_plot,
                              labels = c("A", "B"),
                              ncol = 1, nrow = 2)
ggsave(plot_overview_ch, filename = paste0("./data/figures/Figure3_",format(Sys.time(), "%Y-%m-%d"), ".pdf"), height =6, width = 5,  bg = "transparent")
ggsave(plot_overview_ch, filename = paste0("./data/figures/Figure3_",format(Sys.time(), "%Y-%m-%d"), ".png"), height =6, width = 5,  bg = "transparent")
ggsave(plot_overview_ch, filename = paste0("./data/figures/Figure3_",format(Sys.time(), "%Y-%m-%d")), height =6, width = 5,  bg = "transparent",device='tiff', dpi=700)



