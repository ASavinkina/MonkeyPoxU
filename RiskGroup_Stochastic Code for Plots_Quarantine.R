# Stochastic university monkeypox model
# Used for graphs for presentations and paper

# 11/11/22

# Import libraries

library(adaptivetau)
library(ggplot2)
library(tidyverse)
library(scales)
library(gridExtra)


# Set parameters of interest

# Number of initial infections at day 0 of 100

initial_inf= 1 

# Average number of exogenous infections over 100 days

exogshock =1

# Proportion of cases detected and isolated

diagPerc = 0.5

# Number of student contacts of each detected and isolated case
#   That are quarantined

studentcontacts = 2

#Name scenarios for later saving of plots

isolation = ifelse(diagPerc==0, "NoIsolation",
              ifelse(diagPerc==0.2, "20pIso",
                ifelse(diagPerc==0.5, "50pIso",
                    ifelse(diagPerc==0.8, "80pIso",
                        ifelse(diagPerc==0.9, "90pIso",NA)))))

quarantine = ifelse(studentcontacts==0, "NoQuarantine",
           ifelse(studentcontacts==2, "2studentcontacts",
                  ifelse(studentcontacts==6, "6studentcontacts",
                         ifelse(studentcontacts==20, "20studentcontacts",NA))))

# Additional variable definitions

Yale_undegrad_pop = 6500 #input$Yale_undegrad_pop # total population
exograte = exogshock/100  # rate of exogenous shocks
recoveryrate = 1/21 # recovery rate
diagrate = (diagPerc*recoveryrate)/(1-diagPerc) # daily rate of diagnosis of infectious cases
isorate = 1/28 #duration of isolation
quarrate = 1/14  # duration of quarantined 
R0_h= 2.4 #R0
R0_l = 0.8 #R0

HR_prop = 0.1
LR_prop = 1-HR_prop



# Set initial conditions for model

init.values = c(
  S_h = Yale_undegrad_pop*HR_prop-initial_inf, P_h = 0, I_h = initial_inf, Dx0_h=0, Dx_h=0, 
  Qs_h=0, Qi_h=0, R_h=0,
  S_l = Yale_undegrad_pop*LR_prop, P_l = 0, I_l = 0, Dx0_l=0, Dx_l=0, 
  Qs_l=0, Qi_l=0, R_l=0
)

# Specify all transitions

transitions = list(
  c(S_h = -1, P_h = +1), # movement from susceptible to presymptomatic, endogenous infection
  c(S_h = -1, P_h = +1), # movement from susceptible to presymptomatic, exogenous infection 
  c(S_h= -1, Qs_h= +1), #movement from susceptible to quarantined susceptible
  c(P_h= -1, Qi_h= +1), #movement from presymptomatuc to quarantined infected
  c(Qs_h= -1, S_h= +1), #movement from quarantined susceptible back to susceptible
  c(Qi_h= -1, I_h= +1), #movement from quarantined infected to infected
  c(P_h = -1, I_h = +1), #movement from presymptomatic to infected 
  c(I_h = -1, Dx0_h = +1),  #movement from infected to newly diagnosed
  c(Dx0_h = -1, Dx_h = +1), #movement from newly diagnosed to diagnosed
  c(Qi_h = -1, Dx_h = +1), #movement from quarantined infected to diagnosed
  c(I_h = -1, R_h = +1), #movement from infected to recovered
  c(Dx_h= -1, R_h= +1), #movement from diagnosed to recovered
  
  c(S_l = -1, P_l = +1), # movement from susceptible to presymptomatic, endogenous infection
  c(S_l= -1, Qs_l= +1), #movement from susceptible to quarantined susceptible
  c(P_l= -1, Qi_l= +1), #movement from presymptomatuc to quarantined infected
  c(Qs_l= -1, S_l= +1), #movement from quarantined susceptible back to susceptible
  c(Qi_l= -1, I_l= +1), #movement from quarantined infected to infected
  c(P_l = -1, I_l = +1), #movement from presymptomatic to infected 
  c(I_l = -1, Dx0_l = +1),  #movement from infected to newly diagnosed
  c(Dx0_l = -1, Dx_l = +1), #movement from newly diagnosed to diagnosed
  c(Qi_l = -1, Dx_l = +1), #movement from quarantined infected to diagnosed
  c(I_l = -1, R_l = +1), #movement from infected to recovered
  c(Dx_l= -1, R_l= +1) #movement from diagnosed to recovered
  
)

# Rates for all transitions
# (In same order as "transitions")
RateF <- function(x, pars, times) {
  return(c(
    pars$beta_l*x["S_h"]*x["I_l"]/(x["S_l"] + x["P_l"] + x["I_l"]  + x["R_l"]+ x["S_h"] + x["P_h"] + x["I_h"] + x["R_h"])+
      pars$beta_l*x["S_h"]*x["I_h"]/(x["S_l"] + x["P_l"] + x["I_l"] + x["R_l"]+ x["S_h"] + x["P_h"] + x["I_h"] + x["R_h"])+
      pars$beta_h*x["S_h"]*x["I_h"]/(x["S_h"] + x["P_h"] + x["I_h"] + x["R_h"]),                                                                                
    pars$theta, # movement from susceptible to presymptomatic, exogenous infection
    ifelse(x["P_h"]>1,pars$iota*x["Dx0_h"]*(1-pars$attackrate),pars$iota*x["Dx0_h"]),
    ifelse(x["P_h"]>1, pars$iota*x["Dx0_h"]*pars$attackrate, 0), #movement from susceptible to quarantined infected (based on number newly diagnosed)
    pars$omega*x["Qs_h"], #movement from quarantined susceptible back to susceptible
    pars$omega*x["Qi_h"], #movement from quarantined infected to infected
    pars$gamma*x["P_h"], #movement from presymptomatic to infected (duration of incubation)
    pars$delta*x["I_h"], #movement from infected to diagnosed
    pars$tau*x["Dx0_h"], #movement from newly diagnosed to diagnosed
    pars$mu*x["Qi_h"], #movement from quarantined infected to diagnosed
    pars$rho*x["I_h"], #movement from infected to recovered
    pars$omicron*x["Dx_h"], #movement from diagnosed to recovered 
    
    pars$beta_l*x["S_l"]*x["I_l"]/(x["S_l"] + x["P_l"] + x["I_l"]  + x["R_l"]+ x["S_h"] + x["P_h"] + x["I_h"] + x["R_h"])+
      pars$beta_l*x["S_l"]*x["I_h"]/(x["S_l"] + x["P_l"] + x["I_l"] + x["R_l"]+ x["S_h"] + x["P_h"] + x["I_h"] + x["R_h"]),                                                                         
    ifelse(x["P_l"]>1,pars$iota*x["Dx0_l"]*(1-pars$attackrate),pars$iota*x["Dx0_l"]),
    ifelse(x["P_l"]>1, pars$iota*x["Dx0_l"]*pars$attackrate, 0), #movement from susceptible to quarantined infected (based on number newly diagnosed)
    pars$omega*x["Qs_l"], #movement from quarantined susceptible back to susceptible
    pars$omega*x["Qi_l"], #movement from quarantined infected to infected
    pars$gamma*x["P_l"], #movement from presymptomatic to infected (duration of incubation)
    pars$delta*x["I_l"], #movement from infected to diagnosed
    pars$tau*x["Dx0_l"], #movement from newly diagnosed to diagnosed
    pars$mu*x["Qi_l"], #movement from quarantined infected to diagnosed
    pars$rho*x["I_l"], #movement from infected to recovered
    pars$omicron*x["Dx_l"] #movement from diagnosed to recovered 
  ))
}


# Setting parameters
pars = list(
  beta_h= R0_h * (1/21), #(6*0.1)*(1/14), # beta
  beta_l = R0_l * (1/21),
  gamma= 1/7.6, #incubation period
  delta= diagrate,  #diagnosis rate
  rho= recoveryrate, #recovery rate
  tau = 1, #rate of moving from newly diagnosed to diagnosed
  iota= studentcontacts, # students vaccinated per diagnosed case
  mu = 1/7.6, #incubation period for those in vaccination
  omega = quarrate, #length of quarantine for susceptible
  attackrate=0.2, 
  theta= exograte, #rate of exogenous shocks
  omicron = isorate
)


#Running stochastic model 1,000 times

#  Create dataset



results_all_presymptomatic <- data.frame(matrix(0, nrow=0, ncol=4))
colnames(results_all_presymptomatic) <- c("time", "P_h","P_l", "run")
results_all_infected <- data.frame(matrix(0, nrow=0, ncol=4))
colnames(results_all_infected) <- c("time", "I_h","I_l", "run")
results_all_diagnosed <- data.frame(matrix(0, nrow=0, ncol=4))
colnames(results_all_diagnosed) <- c("time", "Dx_h","Dx_l", "run")
results_all_recovered <- data.frame(matrix(0, nrow=0, ncol=4))
colnames(results_all_recovered) <- c("time", "R_h","R_l", "run")
results_all_newlydiagnosed <- data.frame(matrix(0, nrow=0, ncol=4))
colnames(results_all_newlydiagnosed) <- c("time", "Dx0_h","Dx0_l", "run")
results_all_quarantined<- data.frame(matrix(0, nrow=0, ncol=4))
colnames(results_all_quarantined) <- c("time", "Qs_h", "Qs_l","run")

runs=1000

set.seed(1)

for (i in 1:runs) {
  
  
  results = as.data.frame(ssa.adaptivetau(init.values, 
                                          transitions, 
                                          RateF, 
                                          pars, 
                                          tf=100))
  
  results_presymptomatic_i <- results[,c(1,3,11)]
  results_presymptomatic_i$run <- paste0(i)
  
  results_all_presymptomatic <- rbind(results_all_presymptomatic, results_presymptomatic_i)
  
  results_infected_i <- results[,c(1,4,12)]
  results_infected_i$run <- paste0(i)
  
  results_all_infected <- rbind(results_all_infected, results_infected_i)
  
  results_newlydiagnosed_i <- results[,c(1,5,13)]
  results_newlydiagnosed_i$run <- paste0(i)
  
  results_all_newlydiagnosed <- rbind(results_all_newlydiagnosed, results_newlydiagnosed_i)
  
  results_diagnosed_i <- results[,c(1,6,14)]
  results_diagnosed_i$run <- paste0(i)
  
  results_all_diagnosed <- rbind(results_all_diagnosed, results_diagnosed_i)
  
  results_recovered_i <- results[,c(1,9,17)] 
  results_recovered_i$run <- paste0(i)
  
  results_all_recovered <- rbind(results_all_recovered, results_recovered_i)
  
  results_quarantined_i <- results[,c(1,7,15)] 
  results_quarantined_i[,2] <- results_quarantined_i[,2] + results[,8] 
  results_quarantined_i[,3] <- results_quarantined_i[,3] + results[,16]
  results_quarantined_i$run <- paste0(i)
  
  results_all_quarantined <- rbind(results_all_quarantined, results_quarantined_i)
  
  
}

# Trajectory graphs

# Infectious students plot

# Average infections students

results_all_infected2 <- results_all_infected

results_all_infected2$I <- results_all_infected$I_h + results_all_infected$I_l

results_all_infected2$time2 <- round(results_all_infected2$time)

results_all_infected3 <- results_all_infected2 %>%
  group_by(run,time2) %>%
  summarise_at(vars(I), list(maxI = max))

results_all_infected4 <- results_all_infected3 %>%
  group_by(time2) %>%
  summarise_at(vars(maxI), list(nmin=min, Q1=~quantile(., probs = 0.25), Q95l=~quantile(., probs = 0.05),
                                median=median, Q3=~quantile(., probs = 0.75),Q95u=~quantile(., probs = 0.95),
                                max=max))
# Total infections

total_infections <- results_all_recovered

total_infections$allinfections_hr <- total_infections$R_h + results_all_diagnosed$Dx_h + 
  results_all_newlydiagnosed$Dx0_h +
  results_all_infected$I_h +
  results_all_presymptomatic$P_h 

total_infections$allinfections_lr <- total_infections$R_l + 
  results_all_diagnosed$Dx_l + 
  results_all_newlydiagnosed$Dx0_l +
  results_all_infected$I_l+ 
  results_all_presymptomatic$P_l


total_infections2 <- total_infections[which(total_infections$time==100),]

#   Calculate number of expected index cases

expected_cases <- initial_inf+exogshock

#   Calculate proportion of times there are more cases than the expected number

additional_cases_hr <- length(total_infections2$allinfections_hr[total_infections2$allinfections_hr>expected_cases])
additional_cases_lr <- length(total_infections2$allinfections_lr[total_infections2$allinfections_lr>0])


#   Calculate mean number of additional cases given secondary case

additional_cases_mean_hr <- mean(total_infections2$allinfections_hr[total_infections2$allinfections_hr>expected_cases])
additional_cases_mean_lr <- mean(total_infections2$allinfections_lr[total_infections2$allinfections_lr>expected_cases])



# Average quarantine beds


results_all_quarantined2 <- results_all_quarantined

results_all_quarantined2$Q <- results_all_quarantined2$Qs_h + results_all_quarantined2$Qs_l

results_all_quarantined2$time2 <- round(results_all_quarantined2$time)

results_all_quarantined3 <- results_all_quarantined2 %>%
  group_by(run,time2) %>%
  summarise_at(vars(Q), list(maxQ = max))


results_all_quarantined_likelihood <- results_all_quarantined3 %>%
  group_by(run) %>%
  summarise_at(vars(maxQ), list(maxQQ = max))


results_all_quarantined4 <- results_all_quarantined3 %>%
  group_by(time2) %>%
  summarise_at(vars(maxQ), list(nmin=min, Q1=~quantile(., probs = 0.25), Q95l=~quantile(., probs = 0.05),
                                median=median, Q3=~quantile(., probs = 0.75),Q95u=~quantile(., probs = 0.95),
                                max=max))

# Calculated likelihood of quarantined students beyond the index case, calculate average number of quarantined students,
# and calculate max number of quarantined students

quarantined_number <- sum(results_all_quarantined_likelihood$maxQQ>studentcontacts)
quarantined_average <- mean(results_all_quarantined_likelihood$maxQQ)
quarantined_max <- max(results_all_quarantined_likelihood$maxQQ)


# Average isolated students 


#results_all_diagnosed2 <- results_all_diagnosed 

results_all_diagnosed2 <- cbind(results_all_diagnosed, results_all_newlydiagnosed$Dx0_h, results_all_newlydiagnosed$Dx0_l)

results_all_diagnosed2$total <- results_all_diagnosed2$Dx_h + results_all_diagnosed2$Dx_l +
  results_all_diagnosed2$`results_all_newlydiagnosed$Dx0_h`+ 
  results_all_diagnosed2$`results_all_newlydiagnosed$Dx0_l`

results_all_diagnosed2$time2 <- round(results_all_diagnosed2$time)

results_all_diagnosed3 <- results_all_diagnosed2 %>%
  group_by(run,time2) %>%
  summarise_at(vars(total), list(maxD = max))

results_all_diagnosed_likelihood <- results_all_diagnosed3 %>%
  group_by(run) %>%
  summarise_at(vars(maxD), list(maxDD = max))

# Calculated likelihood of isolated students beyond the index case, calculate average number of isolated students,
# and calculate max number of isolated students

isolated_number <- sum(results_all_diagnosed_likelihood$maxDD>1)
#isolated_average <- mean(results_all_diagnosed_likelihood$maxDD)
isolated_average <- mean(results_all_diagnosed_likelihood$maxDD[results_all_diagnosed_likelihood$maxDD>1])
isolated_max <- max(results_all_diagnosed_likelihood$maxDD)

#Table of Values for number of infections, isolations, and quarantine

Table <- matrix(c("% more than initial cases. high-risk","% cases, low-risk","Mean total cases if more than initial, high risk",
                  "Mean total cases, low risk","% more than 1 isolated case",
                  "mean number of isolated cases","isolated max", "mean number of quarantined people","quarantined max",
                  additional_cases_hr,additional_cases_lr,additional_cases_mean_hr, additional_cases_mean_lr, isolated_number,isolated_average,isolated_max,
                  quarantined_average, quarantined_max), ncol=2)




#Plots

avg_infectious_plot_i <- ggplot(data=results_all_infected4, aes(x=time2, y=median)) + geom_line() +
  theme_classic() + theme(legend.position = "none", text = element_text(size = 20)) + ylim(0,80) +
  geom_ribbon(aes(ymin = Q95l, ymax = Q95u), alpha = 0.1) + xlab("Days") + ylab("Average number of \ninfectious students") +
  ggtitle("Number of infectious \nstudents by day, and 95% interval")

avg_isolated_plot_i <- ggplot(data=results_all_diagnosed4, aes(x=time2, y=median)) + geom_line() +
  theme_classic() + theme(legend.position = "none", text = element_text(size = 20)) + ylim(0,25) +
  geom_ribbon(aes(ymin = Q95l, ymax = Q95u), alpha = 0.1) + xlab("Days") + ylab("Average number of \nisolated students") +
  ggtitle("Number of isolated \nstudents by day, and 95% interval")

avg_quarantine_plot_i <- ggplot(data=results_all_quarantined4, aes(x=time2, y=median)) + geom_line() +
  theme_classic() + theme(legend.position = "none", text = element_text(size = 20)) + ylim(0,250) +
  geom_ribbon(aes(ymin = Q95l, ymax = Q95u), alpha = 0.1) + xlab("Days") + ylab("Average number of \nquarantined students") +
  ggtitle("Number of quarantined \nstudents by day, and 95% interval")

maxquarplot <- ggplot(results_all_quarantined_likelihood, aes(x=maxQQ)) +
    geom_histogram(fill='grey', color='black') +
    theme_classic() + theme(text = element_text(size = 20))+ scale_x_continuous(breaks=c(0,50,100,150,200,250), limits = c(0,250)) +
    xlab('Maximum number of quarantined students') + ylab('Likelihood')

maxisoplot <- ggplot(results_all_diagnosed_likelihood, aes(x=maxDD)) +
  geom_histogram(fill='grey', color='black') +
  theme_classic() + scale_x_continuous(breaks=c(0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100), limits = c(0,100)) +
  xlab('Maximum number of isolated students') + ylab('Likelihood')

infections_plot <- ggplot(total_infections2, aes(x=allinfections)) +
  geom_histogram(fill='grey', color='black', bins=15) +
  theme_classic() + theme(text = element_text(size = 20))+ scale_x_continuous(breaks=c(0,25,50,75,100,125), limits = c(0,125)) +
  xlab('Cumulative number of infections over 100 days') + ylab('Likelihood')


destination = paste0('Plots/','Quar_Plots_',isolation,'_',quarantine,'.pdf')
destination2 = paste0('Plots/','Quar_Table_',isolation,'_',quarantine,'_',initial_inf,'.pdf')


#open PDF
pdf(file=destination)
#specify to save plots in 2x2 grid


#save plots to PDF
plot(avg_infectious_plot_i)
plot(avg_isolated_plot_i)
plot(avg_quarantine_plot_i)
plot(infections_plot)
plot(maxisoplot)
plot(maxquarplot)

#turn off PDF plotting
dev.off() 

pdf(file=destination2)
grid.table(Table)
dev.off() 





