
#' ---
#' title: Imports and the Swiss SARS-CoV-2 Epidemic
#' subtitle: model spread of new variants in Switzerland
#' author: "Martina Reichmuth"
#' date: "'01/12/2021"
#' ---


imports_003_SEIR <- function(imports, variant){ # variant == "alpha" or "delta"
  popsize <- unique(cov_ch$pop[cov_ch$geoRegion=="CH"])
  if(variant=="alpha"){shift <- 3}
  if(variant=="delta"){shift <- 2}
  wt_cases <- round(mean(na.omit(swiss_cov$weigthed_cases[as_date(swiss_cov$date) %in%  ((period[1:14])-shift)])))
  times <- 1:times_length
  R <- swiss_cov$median_R_mean[swiss_cov$date %in% period] #might be taken randomly from interval
  #R_e = (1 - p)*R_w + p*R_v
  #R_v = (1 + kappa)*R_w
  #R_w = R_e/(1 + p*kappa)
  # Settings
  parameters <- c(#beta = R_init/2.6/state[["S"]], # R = beta*S*D = beta*S/gamma -> beta = R*gamma/S
    alpha = as.numeric(parameters_output$median[parameters_output$parameter=="kappa"&parameters_output$variant=="Alpha"]),
    delta = as.numeric(parameters_output$median[parameters_output$parameter=="kappa"&parameters_output$variant=="Delta"]),
    sigma = 1/2.6, # might be adapted with sd and normal distribution
    gamma = 1/2.6,
    zeta = 1/2)
  ascertainment <- 0.5
  introduction <- 0
  names(parameters)[names(parameters) == variant] <- 'variant_advantage'
  #init_cases <- wt_cases/ascertainment*(1/as.numeric(unname(parameters["gamma"]))+1/as.numeric(unname(parameters["zeta"]))) #cases/ascertainment  *incubation time+ time to test
  init_tested <- round(wt_cases/as.numeric(unname(parameters["zeta"])))
  init_cases <- round(wt_cases/(ascertainment*as.numeric(unname(parameters["gamma"])))) #cases/ascertainment  *incubation time+ time to test
  
  recovered <- popsize*seroprevalence_ch$p[seroprevalence_ch$date %in% period][1]#seroprevalence_ge$k/seroprevalence_ge$n
  
  states <- c(S = popsize - 2*init_cases- init_tested - recovered - introduction,#
              E1 = init_cases,
              E2 = 0,
              I1 = init_cases,
              I2 = introduction,
              C1 = init_tested,
              C2= 0,
              T1 = wt_cases,
              T2= 0,
              R1 = recovered,
              R2 = 0,
              Imports=0)
  # Run simulation
  SEIR_variants = function(state, parameter,time){
    require(deSolve) # (Soetaert et al, 2010)
    SEIR <- function(t, state, parameters) {
      with(as.list(c(state, parameters)), {
        p <-  (as.numeric(E2)/(as.numeric(E2)+as.numeric(E1)))
        omega <- imports[pmax(1,ceiling(t))]
        beta1 <- as.numeric(R[pmax(1,ceiling(t))])/(1 + p*variant_advantage)*gamma/S 
        beta2 <- as.numeric(R[pmax(1,ceiling(t))])*(1 + variant_advantage)/(1 + p*variant_advantage)*gamma/S#R[round(t)]
        dS <- - beta1*S*I1 - beta2*S*I2
        dE1 <- beta1*S*I1 - sigma*E1
        dE2 <- beta2*S*I2 - sigma*E2 + omega
        dI1 <- sigma*E1 - gamma*I1
        dI2 <- sigma*E2 - gamma*I2
        dC1 <- gamma*I1 *ascertainment - zeta*C1# got tested after being infectious, e.g., zeta = 1/2, ascertainment 50%
        dC2 <-  gamma*I2 *ascertainment - zeta*C2# got tested after being infectious, e.g., zeta = 1/2, ascertainment 50%
        dT1 <- zeta*C1
        dT2 <- zeta*C2
        dR1 <- gamma*I1*(1-ascertainment)+zeta*C1 - omega
        dR2 <- gamma*I2*(1-ascertainment)+zeta*C2
        dImports <- omega
        der <- c(dS, dE1, dE2, dI1, dI2, dC1, dC2,dT1, dT2,dR1, dR2,dImports)
        list(der)
      })
    }
    out= ode(state, time, SEIR, parameter)
    out[-1,c("T1", "T2","Imports")] <- round(diff(out[,c("T1", "T2","Imports")]))
    return(as.data.frame(out))
  }
  sim <- SEIR_variants(states,parameters,times)
  
  sim$date <- sim$time+min(variant_data$date)-1
  sim$T_total <- sim$T1+sim$T2
  sim$proportion <- sim$T2/sim$T_total
  sim$R_total <- sim$R1+sim$R2
  return(as.data.frame(sim))
}
