




variation_scnearios <- list.files( "./data/105", pattern="csv", full.names=TRUE, recursive=FALSE)
variance_outputs <- c()
for (l in 1:length(variation_scnearios)) {
  new <- read.csv(variation_scnearios[l])
  new$variant <- ifelse(grepl("delta", import_scnearios[l]), "Delta","Alpha")
  variance_outputs <- rbind(variance_outputs,new)
}
import_scnearios <- list.files( "./data/105", pattern=".rds", full.names=TRUE, recursive=FALSE)
import_scnearios_output<- c()
for (l in 1:length(import_scnearios)) {
  new <- readRDS(import_scnearios[l])
  #cumulative incidence limit of model was 1e6 
  new <- new[new$cum_cases<=(5*1e5),]
  new$variant <- ifelse(grepl("delta", import_scnearios[l]), "delta","alpha")
  import_scnearios_output <- rbind(import_scnearios_output,new)
}

import_scnearios_output %>% group_by(seeds) %>% 
  summarise(cum_incidence = quantile(cum_cases, probs=c(probs = c(.025,.5,.975))),
            mean_incidence = mean(cum_cases),
            sd_incidence = sd(cum_cases))

for (n in 1:2) {
  if(n %in% c(1)){
    variant <- "Alpha"
    col_voc <- col_alpha

    }
  if(n %in% c(2)){
     variant <- "Delta"
    col_voc <- col_delta
  }
col_voc <- col_alpha

max_xaxis<- as.numeric(round(quantile(import_scnearios_output$cum_cases[import_scnearios_output$seeds==1], probs = 0.99)))
max_yaxis<- max(variance_outputs$sd)
var_inc_plot <- ggplot()+
  theme_minimal()+
  geom_line(data=variance_outputs, aes( x=(value), y=(sd),group=seeds,color= as.character(seeds)))+
  #facet_wrap(vars(seeds), nrow = 4)+
  scale_color_manual(name="Estimated imports",
                     values = col_voc,
                     labels = c("1", "1:10", "1:100"))+
  coord_cartesian(ylim = c(1, max_yaxis),xlim = c(500,max_xaxis)) +
  labs( x = "Cumulative incidence", y =bquote("Variation \n(standard deviation)"), tag=variant)


if(n %in% c(1)){
  alpha_var_inc_plot <- var_inc_plot
  }

}
variation_incidence_plot <- ggarrange(alpha_var_inc_plot, var_inc_plot,
                               ncol = 2, nrow = 1)

ggsave(variation_incidence_plot, filename = paste0("./data/figures/Figure_variation_incidence_",format(Sys.time(), "%Y-%m-%d"), ".png"), height = 3, width = 10,  bg = "transparent")
ggsave(variation_incidence_plot, filename = paste0("./data/figures/Figure_variation_incidence_",format(Sys.time(), "%Y-%m-%d"), ".pdf"), height = 3, width = 10,  bg = "transparent")



#old

#
import_scnearios <- list.files( "./data/105", pattern=".rds", full.names=TRUE, recursive=FALSE)
import_scnearios_output<- c()
for (l in 1:length(import_scnearios)) {
  new <- readRDS(import_scnearios[l])
  #cumulative incidence limit of model was 1e6 
  new <- new[new$cum_cases<=(5*1e5),]
  new$variant <- ifelse(grepl("delta", import_scnearios[l]), "delta","alpha")
  import_scnearios_output <- rbind(import_scnearios_output,new)
}
import_scnearios_output <- as.data.frame(import_scnearios_output)
colnames(import_scnearios_output)[is.na(colnames(import_scnearios_output))] <- c(1:180)

models_output_num <- as.data.frame(t(import_scnearios_output[,c(8:(length(import_scnearios_output)-1))]))#8
#models_output_info <- as.data.frame((cbind(rownames(import_scnearios_output), import_scnearios_output[,c(1,length(import_scnearios_output))])))#8
#models_output_info$scenario_num <- models_output_info$`rownames(import_scnearios_output)`
  
models_output_num$date <- c(1:180)
models_output_num<- melt(models_output_num, id.vars=c("date"))
#models_output_num$scenario_num <- rep(rownames(import_scnearios_output),each=180)

models_output_num$variant <- c(rep("Alpha",180*sum(grepl("lpha",import_scnearios_output$variant))),rep("Delta",180*sum(grepl("elta",import_scnearios_output$variant))))
models_output_num$seeds <- c(rep(1,180*sum(import_scnearios_output$seeds %in%1)),rep(10,180*sum(import_scnearios_output$seeds %in%10)),rep(100,180*sum(import_scnearios_output$seeds %in%100)))


for (n in 1:2) {
  if(n %in% c(1)){
    # Alpha imports from Emma Hodcroft
    import <- imports_alpha
    import$approach <- ifelse(import$HeadOfClade!="True" & is.na(import$Clade)| import$HeadOfClade=="True" & !is.na(import$Clade)  & import$Clade>-1,"conservative", "liberal" )
    variant_data <- alpha_ch
    variant <- "Alpha"
    col_voc <- col_alpha
    models_output <- models_output_num[grepl("lpha",import_scnearios_output$variant),]
  }
  if(n %in% c(2)){
    # Delta imports from Emma Hodcroft
    import <- imports_delta
    import$approach <- ifelse(import$HeadOfClade!="True" & is.na(import$Clade)| import$HeadOfClade=="True" & !is.na(import$Clade)  & import$Clade>-1,"conservative", "liberal" )
    variant_data <- delta_ch
    variant <- "Delta"
    col_voc <- col_delta
    models_output <- models_output_num[grepl("elta",import_scnearios_output$variant),]
    
  }
first_import_date <- as_date(min(import$MinDateSwissChildren))-7
last_import_date <- max(as_date(import$MinDateSwissChildren))-7
length_imports <- length(first_import_date:last_import_date)

models_output$time <- models_output$date
models_output$date <-  models_output$date+first_import_date

max_incidence_perseed<- models_output %>% group_by(seeds) %>% summarise(max_incidence_perseed = max(value))
variance_outputs <- data.frame(array(NA, dim = c(sum(max_incidence_perseed$max_incidence_perseed)+3, 4)))
colnames(variance_outputs) <- c("incidence","variance","sd", "seeds")

for(k in 1:3){
   i <- c(1,10,100)[k]
  subset_variance_outputs <- models_output[models_output$seeds==i,]
  l <- sum(max_incidence_perseed$max_incidence_perseed[0:(k -1)])
  for(j in 1:max(subset_variance_outputs$value)){
    variance_outputs$variance[l+j+1] <- var(subset_variance_outputs$time[subset_variance_outputs$value==j])
    variance_outputs$sd[l+j+1] <- sd(subset_variance_outputs$time[subset_variance_outputs$value==j])
    variance_outputs$incidence[l+j+1] <- j
    variance_outputs$seeds[l+j+1] <- i
  }
}
variance_outputs<- variance_outputs[variance_outputs$seeds!=0,]
variance_outputs<- variance_outputs[!is.na(variance_outputs$seeds),]
var_inc_plot <- ggplot()+
  theme_minimal()+
  geom_line(data=variance_outputs, aes( x=(incidence), y=(variance),group=seeds,color= as.character(seeds)))+
  scale_color_manual(name="Estimated imports",
                     values = col_voc,
                     labels = c("1", "1:10", "1:100"))+
  coord_cartesian(ylim = c(0, 250),xlim = c(0,1500)) +
  labs(tag="", x = "Incidence", y =bquote("Variance"))


inc_time_plot <- ggplot()+
  theme_minimal()+
  geom_line(data=models_output, aes( x=(date), y=(value),group=scenario_num,color= as.character(seeds)))+
  scale_color_manual(name="Estimated imports",
                     values = col_voc,
                     labels = c("1", "1:10", "1:100"))+
  coord_cartesian(ylim = c(0, 1000),xlim = c(0,1500)) +
  labs(tag="", x = "", y =paste0("Number of",variant," cases"))


if(n %in% c(1)){
  alpha_var_inc_plot <- var_inc_plot
  alpha_inc_time_plot<- inc_time_plot
  alphavariance_outputs <- variance_outputs
}
  
}
stochastic_plot <- grid.arrange(alpha_inc_time_plot, alpha_var_inc_plot,
                                inc_time_plot, var_inc_plot, nrow = 2)


ggsave(stochastic_plot, filename = paste0("./data/figures/stochastic_figure",format(Sys.time(), "%b"), ".pdf"), height = 6, width = 8,  bg = "transparent")
ggsave(stochastic_plot, filename = paste0("./data/figures/stochastic_figure",format(Sys.time(), "%b"), ".png"), height = 6, width = 8,  bg = "transparent")


models_output[,c("variable","cum_inc")]<- models_output %>% group_by(variable) %>% summarise(cum_inc=cumsum(value))


max_cumincidence_perseed<- models_output %>% group_by(seeds) %>% summarise(max_cumincidence_perseed = max(cum_inc))

variance_outputs <- data.frame(array(NA, dim = c(sum(max_cumincidence_perseed$max_cumincidence_perseed)+3, 5)))
colnames(variance_outputs) <- c("value","variance","sd","range", "seeds")

for(k in 1:3){
 # k <- 3
  i <- c(1,10,100)[k]
  subset_variance_outputs <- models_output[models_output$seeds==i,]
  l <- sum(max_cumincidence_perseed$max_cumincidence_perseed[0:(k -1)])
  #cum_values <-  1:max(max_cumincidence_perseed$max_cumincidence_perseed[k])[1:max(max_cumincidence_perseed$max_cumincidence_perseed[k]) %in% unique(subset_variance_outputs$cum_inc)]
  for(j in 500:2000){
    #s<- sort(as.numeric(unique(subset_variance_outputs$cum_inc)))[j]
    time_j<- subset_variance_outputs[subset_variance_outputs$cum_inc%in%j,]
    time_j <- time_j$time[!duplicated(time_j$variable)]
    time_j<- time_j - min(time_j)
    variance_outputs$variance[l+j+1] <- var(time_j)
    variance_outputs$sd[l+j+1] <- sd(time_j)
    variance_outputs$range[l+j+1] <- max(time_j)
    variance_outputs$value[l+j+1] <- j
    variance_outputs$seeds[l+j+1] <- i
  }
}

variance_outputs<- variance_outputs[variance_outputs$seeds!=0,]
variance_outputs<- variance_outputs[!is.na(variance_outputs$seeds),]
var_inc_plot <- ggplot()+
  theme_minimal()+
  geom_line(data=variance_outputs, aes( x=(value), y=(sd),group=seeds,color= as.character(seeds)))+
  #facet_wrap(vars(seeds), nrow = 4)+
  scale_color_manual(name="Estimated imports",
                     values = col_voc,
                     labels = c("1", "1:10", "1:100"))+
  #coord_cartesian(ylim = c(1, 250),xlim = c(1,1500)) +
  labs(tag="", x = "Incidence", y =bquote("Variation \n(standard deviation of time)"))


#models_output <- models_output[models_output$time <=length_imports,]

# Fit cases to generalized negative binomial model
library(MASS)
r_all <- c()
fit <- glm.nb(T2 ~ time , data = sim_alpha[sim_alpha$approach=="liberal"&sim_alpha$date>first_import_date&sim_alpha$date<last_import_date,]) 
r_all <- c(coef(summary(fit))[2, 1]- 1.96 * coef(summary(fit))[2, 2],coef(summary(fit))[2, 1],coef(summary(fit))[2, 1]+ 1.96 * coef(summary(fit))[2, 2])
model_se <- summary(fit)$SE.theta

# within expectation of number of VoC cases (that were simulated)
probs = c(.025,.5,.975)
n_sim <- 1e4
cases <- numeric(n_sim)
# potential final incidence
cum <- sim_alpha$T2[sim_alpha$approach=="liberal"&sim_alpha$proportion<=0.5]#sim_alpha$proportion<=0.5]
sim_alpha$time[sim_alpha$approach=="liberal"&sim_alpha$proportion<=0.5]
final <- tail(cum, 7)
round(mean(final))
sum(cum)
cum_inc<- final_inc <- numeric(n_sim)
final_inc <- numeric(n_sim)
for(j in 1:n_sim) {
  cum_inc[j] <- sum(rnbinom(length(cum), size =model_se, mu = cum))
  final_inc[j] <- sum(rnbinom(length(final), size = model_se, mu = final))
}
cum_final_expected <- round(rbind(quantile(cum_inc, probs), quantile(final_inc, probs)/7))
cum_expected_seq <- c(cum_final_expected[1,1]:cum_final_expected[1,3])
final_expected_seq <- c(cum_final_expected[2,1]:cum_final_expected[2,3])

models_output$accepted_min[as.numeric(models_output$cum_cases) >= cum_final_expected[1,1] & as.numeric(models_output$final_incidence) >= cum_final_expected[2,1]] <-1
simulation_accepted_min <- import_scnearios_output[import_scnearios_output$accepted_min%in%1,]
table(simulation_accepted_min$seeds)/1e5*100

test <- models_output %>% group_by(variable) %>% summarise(accepted=cumsum(as.numeric(models_output$value)) %in%cum_expected_seq)




models_output$simulation_accepted[as.numeric(models_output$cum_cases) %in% cum_expected_seq & as.numeric(models_output$final_incidence) %in% final_expected_seq] <-1
simulation_accept <- import_scnearios_output[import_scnearios_output$simulation_accepted%in%1,]
table(simulation_accept$seeds)/1e5*100

simulation_accept_num <- as.data.frame(t(simulation_accept[,c(8:(length(simulation_accept)-2))]))#8
simulation_accept_info <- as.data.frame((cbind(rownames(simulation_accept), simulation_accept[,c(1,length(simulation_accept)-1)])))#8
simulation_accept_info$scenario_num <- rownames(simulation_accept)

simulation_accept_num$time <- c(1:180)
simulation_accept_num<- melt(simulation_accept_num, id.vars=c("time"))
simulation_accept_num$scenario_num <- rep(rownames(simulation_accept),each=180)
simulation_accept_num <- merge(simulation_accept_info[,c("scenario_num","seeds", "variant")],simulation_accept_num, by="scenario_num", y.all=TRUE)
simulation_accept_num$date <- as_date(min(imports_alpha$MinDateSwissChildren))-7+simulation_accept_num$time


# also plot VoC curve from simulation and prediction?
ggplot()+
  theme_minimal()+
  geom_line(data=models_output[models_output$variant=="alpha",], aes( x=(date), y=(value),group=scenario_num,color= as.character(seeds)), alpha=0.8)+
scale_color_manual(name="Estimated imports",
                  values = col_alpha,
                  labels = c("1", "1:10", "1:100"))+
  labs(tag="", x = "", y =bquote("Number of Alpha cases"))

  


# time until final and cumulative size are reached
# add imports to counts!
import_scnearios_output$accepted_min[as.numeric(import_scnearios_output$cum_cases) >= cum_final_expected[1,1] & as.numeric(import_scnearios_output$final_incidence) >= cum_final_expected[2,1]] <-1
simulation_accepted_min <- import_scnearios_output[import_scnearios_output$accepted_min%in%1,]
table(simulation_accepted_min$seeds)/1e5*100

simulation_accept_min_num <- as.data.frame(t(simulation_accepted_min[,c(8:(length(simulation_accepted_min)-3))]))#8
simulation_accept_min_info <- as.data.frame((cbind(rownames(simulation_accepted_min), simulation_accepted_min[,c(1,length(simulation_accepted_min)-2)])))#8
simulation_accept_min_info$scenario_num <- rownames(simulation_accepted_min)

simulation_accept_min_num$time <- c(1:180)
simulation_accept_min_num<- melt(simulation_accept_min_num, id.vars=c("time"))
simulation_accept_min_num$scenario_num <- rep(rownames(simulation_accepted_min),each=180)
simulation_accept_min_num <- merge(simulation_accept_min_info[,c("scenario_num","seeds", "variant")],simulation_accept_min_num, by="scenario_num", y.all=TRUE)
simulation_accept_min_num$date <- as_date(min(imports_alpha$MinDateSwissChildren))-7+simulation_accept_min_num$time
time_start_import<- as_date(min(imports_alpha$MinDateSwissChildren))-7

simulation_accept_min_num1 <- simulation_accept_min_num[simulation_accept_min_num$value >= cum_final_expected[2,2],]
simulation_accept_min_num1<- simulation_accept_min_num1[!duplicated(simulation_accept_min_num1$scenario_num),]
time_start_import+quantile(simulation_accept_min_num1$time[simulation_accept_min_num1$variant=="alpha"], probs)
time_start_import+quantile(simulation_accept_min_num1$time[simulation_accept_min_num1$variant=="alpha"&simulation_accept_min_num1$seeds==1], probs)
time_start_import+quantile(simulation_accept_min_num1$time[simulation_accept_min_num1$variant=="alpha"&simulation_accept_min_num1$seeds==10], probs)
time_start_import+quantile(simulation_accept_min_num1$time[simulation_accept_min_num1$variant=="alpha"&simulation_accept_min_num1$seeds==100], probs)

paste(min(alpha_ch$date[alpha_ch$prediction>=0.5]), "(",min(alpha_ch$date[alpha_ch$prediction.lower>=0.5]), min(alpha_ch$date[alpha_ch$prediction.upper>=0.5]),")")

# time graph with prediction, simulation and time by different scenarios

