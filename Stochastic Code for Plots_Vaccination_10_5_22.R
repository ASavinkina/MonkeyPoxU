# Stochastic university monkeypox model
# Used for graphs for presentations and paper

# 10/6/22

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

exogshock =0 

  # Proportion of cases detected and isolated

diagPerc = 0.5

  # Number of student contacts of each detected and isolated case
  #   That are vaccinated

studentcontacts = 0 

#Name scenarios for later saving of plots

isolation = ifelse(diagPerc==0, "NoIsolation",
              ifelse(diagPerc==0.2, "20pIso",
                ifelse(diagPerc==0.5, "50pIso",
                    ifelse(diagPerc==0.8, "80pIso",
                        ifelse(diagPerc==0.9, "90pIso",NA)))))

vaccination = ifelse(studentcontacts==0, "NoVaccination",
                     ifelse(studentcontacts==1, "1studentcontacts",      
           ifelse(studentcontacts==2, "2studentcontacts",
                  ifelse(studentcontacts==6, "6studentcontacts",
                         ifelse(studentcontacts==20, "20studentcontacts",NA)))))


# Additional variable definitions

Vaccine_Efficacy = 0.8 # Vaccine effectiveness for the non-infected
Vaccine_Efficacy_PEP = Vaccine_Efficacy * 0.5 # Vaccine effectiveness for those who are infected

Yale_undegrad_pop = 6500 #input$Yale_undegrad_pop # total population
exograte = exogshock/100  # rate of exogenous shocks
recoveryrate = 1/21 # recovery rate
diagrate = (diagPerc*recoveryrate)/(1-diagPerc) # daily rate of diagnosis of infectious cases
isorate = 1/28 #duration of isolation
quarrate = 1/14  # duration of vaccination 
R0_h= 1.4 #R0


# Set initial conditions for model

init.values = c(
  S_h = Yale_undegrad_pop-initial_inf, P_h = 0, I_h = initial_inf, Dx0_h=0, Dx_h=0, 
  Vs_h=0, Vi_h=0, R_h=0
)

# Specify all transitions

transitions = list(
  c(S_h = -1, P_h = +1), # movement from susceptible to presymptomatic, endogenous infection
  c(S_h = -1, P_h = +1), # movement from susceptible to presymptomatic, exogenous infection 
  c(S_h= -1, Vs_h= +1), #movement from susceptible to vaccinated susceptible
  c(P_h= -1, Vi_h= +1), #movement from presymptomatuc to vaccinated infected
  #c(Qs_h= -1, S_h= +1), #movement from vaccinated susceptible back to susceptible
  #c(Qi_h= -1, I_h= +1), #movement from vaccinated infected to infected
  c(P_h = -1, I_h = +1), #movement from presymptomatic to infected 
  c(I_h = -1, Dx0_h = +1),  #movement from infected to newly diagnosed
  c(Dx0_h = -1, Dx_h = +1), #movement from newly diagnosed to diagnosed
  #c(Qi_h = -1, Dx_h = +1), #movement from vaccinated infected to diagnosed
  c(I_h = -1, R_h = +1), #movement from infected to recovered
  c(Dx_h= -1, R_h= +1) #movement from diagnosed to recovered
  
)

# Rates for all transitions
# (In same order as "transitions")
RateF <- function(x, pars, times) {
  return(c(
    pars$beta_h*x["S_h"]*x["I_h"]/(x["S_h"] + x["P_h"] + x["I_h"] + x["R_h"]), # movement from susceptible to presymptomatic, endogenous infection
    pars$theta, # movement from susceptible to presymptomatic, exogenous infection
    ifelse(x["P_h"]>1,pars$iota*x["Dx0_h"]*(1-pars$attackrate)*pars$Vaccine_Efficacy,pars$iota*x["Dx0_h"]*pars$Vaccine_Efficacy),
    #pars$iota*x["Dx0_h"]*(1-pars$attackrate), #movement from susceptible to vaccinated susceptible (based on number newly diagnosed)
    ifelse(x["P_h"]>1, pars$iota*x["Dx0_h"]*pars$attackrate*pars$Vaccine_Efficacy_PEP, 0), #movement from susceptible to vaccinated infected (based on number newly diagnosed)
    #pars$iota*x["Dx0_h"]*pars$attackrate
    #pars$omega*x["Qs_h"], #movement from vaccinated susceptible back to susceptible
    #pars$omega*x["Qi_h"], #movement from vaccinated infected to infected
    pars$gamma*x["P_h"], #movement from presymptomatic to infected (duration of incubation)
    pars$delta*x["I_h"], #movement from infected to diagnosed
    pars$tau*x["Dx0_h"], #movement from newly diagnosed to diagnosed
    #pars$mu*x["Qi_h"], #movement from vaccinated infected to diagnosed
    pars$rho*x["I_h"], #movement from infected to recovered
    pars$omicron*x["Dx_h"] #movement from diagnosed to recovered 
  ))
}


# Setting parameters
pars = list(
  beta_h= R0_h * (1/21), #(6*0.1)*(1/14), # beta
  gamma= 1/7.6, #incubation period
  delta= diagrate,  #diagnosis rate
  rho= recoveryrate, #recovery rate
  tau = 1, #rate of moving from newly diagnosed to diagnosed
  iota= studentcontacts, # students vaccinated per diagnosed case
  mu = 1/7.6, #incubation period for those in vaccination
  omega = quarrate, #length of vaccination for susceptible
  attackrate=0.2, 
  theta= exograte, #rate of exogenous shocks
  omicron = isorate,
  Vaccine_Efficacy=Vaccine_Efficacy,
  Vaccine_Efficacy_PEP = Vaccine_Efficacy_PEP
)


#Running stochastic model 1,000 times

#  Create dataset



results_all_presymptomatic <- data.frame(matrix(0, nrow=0, ncol=3))
colnames(results_all_presymptomatic) <- c("time", "P_h", "run")
results_all_infected <- data.frame(matrix(0, nrow=0, ncol=3))
colnames(results_all_infected) <- c("time", "I_h", "run")
results_all_diagnosed <- data.frame(matrix(0, nrow=0, ncol=3))
colnames(results_all_diagnosed) <- c("time", "Dx_h", "run")
results_all_recovered <- data.frame(matrix(0, nrow=0, ncol=3))
colnames(results_all_recovered) <- c("time", "R_h", "run")
results_all_newlydiagnosed <- data.frame(matrix(0, nrow=0, ncol=3))
colnames(results_all_newlydiagnosed) <- c("time", "Dx0_h", "run")
results_all_vaccinated<- data.frame(matrix(0, nrow=0, ncol=3))
colnames(results_all_vaccinated) <- c("time", "Dx0_h", "run")

runs=1000

set.seed(1)

for (i in 1:runs) {
  
  
  results = as.data.frame(ssa.adaptivetau(init.values, 
                                          transitions, 
                                          RateF, 
                                          pars, 
                                          tf=100))
  
  results_presymptomatic_i <- results[,c(1,3)]
  results_presymptomatic_i$run <- paste0(i)
  
  results_all_presymptomatic <- rbind(results_all_presymptomatic, results_presymptomatic_i)
  
  results_infected_i <- results[,c(1,4)]
  results_infected_i$run <- paste0(i)
  
  results_all_infected <- rbind(results_all_infected, results_infected_i)
  
  results_newlydiagnosed_i <- results[,c(1,5)]
  results_newlydiagnosed_i$run <- paste0(i)
  
  results_all_newlydiagnosed <- rbind(results_all_newlydiagnosed, results_newlydiagnosed_i)
  
  results_diagnosed_i <- results[,c(1,6)]
  results_diagnosed_i$run <- paste0(i)
  
  results_all_diagnosed <- rbind(results_all_diagnosed, results_diagnosed_i)
  
  results_recovered_i <- results[,c(1,9)] 
  results_recovered_i$run <- paste0(i)
  
  results_all_recovered <- rbind(results_all_recovered, results_recovered_i)
  
  results_vaccinated_i <- results[,c(1,7)]/Vaccine_Efficacy 
  results_vaccinated_i[,2] <- results_vaccinated_i[,2] + results[,8]/Vaccine_Efficacy_PEP 
  results_vaccinated_i$run <- paste0(i)
  
  results_all_vaccinated <- rbind(results_all_vaccinated, results_vaccinated_i)
  
  
}

# Trajectory graphs

# Infectious students plot

# Average infections students

results_all_infected2 <- results_all_infected

results_all_infected2$time2 <- round(results_all_infected2$time)

results_all_infected3 <- results_all_infected2 %>%
  group_by(run,time2) %>%
  summarise_at(vars(I_h), list(maxI = max))

results_all_infected4 <- results_all_infected3 %>%
  group_by(time2) %>%
  summarise_at(vars(maxI), list(nmin=min, Q1=~quantile(., probs = 0.25), Q95l=~quantile(., probs = 0.05),
                                median=median, Q3=~quantile(., probs = 0.75),Q95u=~quantile(., probs = 0.95),
                                max=max))

# Total infections

total_infections <- results_all_recovered

total_infections$allinfections <- total_infections$R_h + results_all_diagnosed$Dx_h + results_all_newlydiagnosed$Dx0_h +
  results_all_infected$I_h + results_all_presymptomatic$P_h

#   Calculate total number of infectious students 

total_infections2 <- total_infections[which(total_infections$time==100),]

#   Calculate number of expected index cases

expected_cases <- initial_inf+exogshock

#   Calculate proportion of times there are more cases than the expected number

additional_cases <- length(total_infections2$allinfections[total_infections2$allinfections>expected_cases])

#   Calculate mean number of additional cases given secondary case

additional_cases_mean <- mean(total_infections2$allinfections[total_infections2$allinfections>expected_cases])

  
# Average vaccination beds

results_all_vaccinated2 <- results_all_vaccinated

results_all_vaccinated2$time2 <- round(results_all_vaccinated2$time)

results_all_vaccinated3 <- results_all_vaccinated2 %>%
  group_by(run,time2) %>%
  summarise_at(vars(Vs_h), list(maxQ = max))


results_all_vaccinated_likelihood <- results_all_vaccinated3 %>%
  group_by(run) %>%
  summarise_at(vars(maxQ), list(maxQQ = max))


results_all_vaccinated4 <- results_all_vaccinated3 %>%
  group_by(time2) %>%
  summarise_at(vars(maxQ), list(nmin=min, Q1=~quantile(., probs = 0.25), Q95l=~quantile(., probs = 0.05),
                                median=median, Q3=~quantile(., probs = 0.75),Q95u=~quantile(., probs = 0.95),
                                max=max))

# Calculated likelihood of vaccinated students beyond the index case, calculate average number of vaccinated students,
# and calculate max number of vaccinated students

vaccinated_number <- sum(results_all_vaccinated_likelihood$maxQQ>studentcontacts)
vaccinated_average <- mean(results_all_vaccinated_likelihood$maxQQ)
vaccinated_max <- max(results_all_vaccinated_likelihood$maxQQ)


# Average isolated students 


results_all_diagnosed2 <- results_all_diagnosed 

results_all_diagnosed2 <- cbind(results_all_diagnosed, results_all_newlydiagnosed$Dx0_h)

results_all_diagnosed2$total <- results_all_diagnosed2$Dx_h + results_all_diagnosed2$`results_all_newlydiagnosed$Dx0_h`

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


results_all_diagnosed4 <- results_all_diagnosed3 %>%
  group_by(time2) %>%
  summarise_at(vars(maxD), list(nmin=min, Q1=~quantile(., probs = 0.25), Q95l=~quantile(., probs = 0.05),
                                median=median, Q3=~quantile(., probs = 0.75),Q95u=~quantile(., probs = 0.95),
                                max=max))


#Table of Values for number of infections, isolations, and vaccinatins

Table <- matrix(c("% more than 1 case","Mean total cases if more than 1","% more than 1 isolated case",
                  "mean number of isolated cases","isolated max", "mean number of vaccinated people","vaccinated max",
  additional_cases,additional_cases_mean,isolated_number,isolated_average,isolated_max,vaccinated_average, vaccinated_max), ncol=2)

#Plots

# Infections plot for no vaccination
# avg_infectious_plot_i <- ggplot(data=results_all_infected4, aes(x=time2, y=median)) + geom_line() +
#   theme_classic() + theme(legend.position = "none",text = element_text(size = 20)) + ylim(0,80) +
#   geom_ribbon(aes(ymin = Q95l, ymax = Q95u), alpha = 0.1) + xlab("Days") + ylab("Average number of \ninfectious students") 
#   #+ ggtitle("Number of infectious \nstudents by day, and 95% interval")

# # Infections plot for vaccination
# avg_infectious_plot_i <- ggplot(data=results_all_infected4, aes(x=time2, y=median)) + geom_line() +
#   theme_classic() + theme(legend.position = "none",text = element_text(size = 20)) + ylim(0,25) +
#   geom_ribbon(aes(ymin = Q95l, ymax = Q95u), alpha = 0.1) + xlab("Days") + ylab("Average number of \ninfectious students") 
# #+ ggtitle("Number of infectious \nstudents by day, and 95% interval")

# Infections plot for vaccination
avg_infectious_plot_i <- ggplot(data=results_all_infected3, aes(x=time2, y=maxI, group=run)) + geom_line(alpha=0.3) +
  theme_classic() + theme(legend.position = "none",text = element_text(size = 20)) + ylim(0,25) + xlab("Days") + ylab("Number of \ninfectious students") 
 #+ ggtitle("Number of infectious \nstudents by day")


avg_isolated_plot_i <- ggplot(data=results_all_diagnosed4, aes(x=time2, y=median)) + geom_line() +
  theme_classic() + theme(legend.position = "none",text = element_text(size = 20)) + ylim(0,10) +
  geom_ribbon(aes(ymin = Q95l, ymax = Q95u), alpha = 0.1) + xlab("Days") + ylab("Average number of \nisolated students") 
#+ggtitle("Number of isolated \nstudents by day, and 95% interval")

avg_vaccination_plot_i <- ggplot(data=results_all_vaccinated4, aes(x=time2, y=median)) + geom_line() +
  theme_classic() + theme(legend.position = "none",text = element_text(size = 20)) + ylim(0,250) +
  geom_ribbon(aes(ymin = Q95l, ymax = Q95u), alpha = 0.1) + xlab("Days") + ylab("Average number of \nvaccinated students") +
  ggtitle("Number of vaccinated \nstudents by day, and 95% interval")

maxvaxplot <- ggplot(results_all_vaccinated_likelihood, aes(x=maxQQ)) +
    geom_histogram(fill='grey', color='black') +
    theme_classic() +theme(text = element_text(size = 20))+ scale_x_continuous(breaks=c(0,50,100,150,200,250), limits = c(0,250)) +
    xlab('Maximum number of vaccinated students') + ylab('Likelihood')

maxisoplot <- ggplot(results_all_diagnosed_likelihood, aes(x=maxDD)) +
  geom_histogram(fill='grey', color='black') +
  theme_classic() +theme(text = element_text(size = 20)) + scale_x_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10), limits = c(0,10)) +
  xlab('Maximum number of isolated students') + ylab('Likelihood')

# Infections plot for no vaccination 

# infections_plot <- ggplot(total_infections2, aes(x=allinfections)) +
#   geom_histogram(fill='grey', color='black', bins=15) +
#   theme_classic() + theme(text = element_text(size = 20))+ scale_x_continuous(breaks=c(0,50,100,150,200,250,300,350), limits = c(0,350)) +
#   xlab('Cumulative number of infections over 100 days') + ylab('Likelihood')

# Infections plot for vaccination scenario

infections_plot <- ggplot(total_infections2, aes(x=allinfections)) +
  geom_histogram(fill='grey', color='black', bins=30) +
  theme_classic() + theme(text = element_text(size = 20))+ scale_x_continuous(breaks=c(0,5,10,15,20,25), limits = c(0,25)) +
  ylim(0,1000) +
  xlab('Cumulative number of infections over 100 days') + ylab('Likelihood') 


destination = paste0('Plots/','Vax_Plots_',isolation,'_',vaccination,'_',exogshock,'_',initial_inf,'_',runs,'.pdf')
destination2 = paste0('Plots/','Vax_Table_',isolation,'_',vaccination,'_',exogshock,'_',initial_inf,'_',runs,'.pdf')

#open PDF
pdf(file=destination)
#specify to save plots in 2x2 grid


#save plots to PDF
plot(avg_infectious_plot_i)
plot(avg_isolated_plot_i)
plot(avg_vaccination_plot_i)
plot(infections_plot)
plot(maxisoplot)
plot(maxvaxplot)

#turn off PDF plotting
dev.off() 

pdf(file=destination2)
grid.table(Table)
dev.off() 

#citation("adaptivetau")

