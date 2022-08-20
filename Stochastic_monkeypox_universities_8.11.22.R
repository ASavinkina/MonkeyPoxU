library(adaptivetau)
library(ggplot2)
library(tidyverse)
library(scales)


# Intiial conditions

Yale_undegrad_pop = 6500 # total population
initial_inf = 10 # initial infections
exogshock =0 #size of exogenous shock. 0 for no exogenous shocks, >1 for superspreader events
exograte = 1/30 # rate of exogenous shocks
diagrate = 1/10 # daily rate of diagnosis of infectious cases
quarduration = 1/14 # duration of quarantine for non-infected
studentcontacts = 5 # students quarantined per diagnosed case
R0_h = 1.5

init.values = c(
  S_h = Yale_undegrad_pop-initial_inf, P_h = 0, I_h = initial_inf, Dx0_h=0, Dx_h=0, 
  Qs_h=0, Qi_h=0, R_h=0
)

# Specify all transitions
transitions = list(
  c(S_h = -1, P_h = +1), # movement from susceptible to presymptomatic, endogenous infection
  c(S_h = -exogshock, P_h = +exogshock), # movement from susceptible to presymptomatic, exogenous infection 
  c(S_h= -1, Qs_h= +1), #movement from susceptible to quarantined susceptible
  c(P_h= -1, Qi_h= +1), #movement from susceptible to quarantined infected
  c(Qs_h= -1, S_h= +1), #movement from quarantined susceptible back to susceptible
  c(P_h = -1, I_h = +1), #movement from presymptomatic to infected 
  c(I_h = -1, Dx0_h = +1),  #movement from infected to newly diagnosed
  c(Dx0_h = -1, Dx_h = +1), #movement from newly diagnosed to diagnosed
  c(Qi_h = -1, Dx_h = +1), #movement from quarantined infected to diagnosed
  c(I_h = -1, R_h = +1), #movement from infected to recovered
  c(Dx_h= -1, R_h= +1) #movement from diagnosed to recovered
  
)

# Rates for all transitions
# (In same order as "transitions")
RateF <- function(x, pars, times) {
  return(c(
    pars$beta_h*x["S_h"]*x["I_h"]/(x["S_h"] + x["P_h"] + x["I_h"] + x["R_h"]), # movement from susceptible to presymptomatic, endogenous infection
    pars$theta, # movement from susceptible to presymptomatic, exogenous infection
    ifelse(x["P_h"]>1,pars$iota*x["Dx0_h"]*(1-pars$attackrate),pars$iota*x["Dx0_h"]),
    #pars$iota*x["Dx0_h"]*(1-pars$attackrate), #movement from susceptible to quarantined susceptible (based on number newly diagnosed)
    ifelse(x["P_h"]>1, pars$iota*x["Dx0_h"]*pars$attackrate, 0), #movement from susceptible to quarantined infected (based on number newly diagnosed)
    #pars$iota*x["Dx0_h"]*pars$attackrate
    pars$omega*x["Qs_h"], #movement from quarantined susceptible back to susceptible
    pars$gamma*x["P_h"], #movement from presymptomatic to infected (duration of incubation)
    pars$delta*x["I_h"], #movement from infected to diagnosed
    pars$tau*x["Dx0_h"], #movement from newly diagnosed to diagnosed
    pars$mu*x["Qi_h"], #movement from quarantined infected to diagnosed
    pars$rho*x["I_h"], #movement from infected to recovered
    pars$rho*x["Dx_h"] #movement from diagnosed to recovered 
  ))
}

# Setting parameters
pars = list(
  beta_h= R0_h * (1/14), #(6*0.1)*(1/14), # beta
  gamma= 1/14, #incubation period
  delta= diagrate,  #diagnosis rate
  rho= 1/14, #recovery rate
  tau = 1, #rate of moving from newly diagnosed to diagnosed
  iota= studentcontacts, # students quarantined per diagnosed case
  mu = 1/14, #incubation period for those in quarantine
  omega = quarduration, #length of quarantine for susceptible
  attackrate=0.1, 
  theta= exograte #rate of exogenous shocks
)


# Running stochastic model

results = as.data.frame(ssa.adaptivetau(init.values, 
                                        transitions, 
                                        RateF, 
                                        pars, 
                                        tf=100))

#Running stochastic model 1,000 times

#  Create dataset

runs=1000


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
results_all_quarantined<- data.frame(matrix(0, nrow=0, ncol=3))
colnames(results_all_quarantined) <- c("time", "Dx0_h", "run")

runs=100

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
  
  results_quarantined_i <- results[,c(1,7)] 
  results_quarantined_i[,2] <- results_quarantined_i[,2] + results[,8] 
  results_quarantined_i$run <- paste0(i)
  
  results_all_quarantined <- rbind(results_all_quarantined, results_quarantined_i)
  
  
}


# # Quarantine beds version 1: likelihood of needed beds by time 
# 
# results_all_quarantined$summary <- ifelse(results_all_quarantined$Qs_h==0, "Zero",
#                                           ifelse(results_all_quarantined$Qs_h>50, "Greater than 50",
#                                                  ifelse(results_all_quarantined$Qs_h>30, "30 to 40",
#                                                         ifelse(results_all_quarantined$Qs_h>20, "20 to 30",
#                                                                ifelse(results_all_quarantined$Qs_h>10, "10 to 20",
#                                                                       ifelse(results_all_quarantined$Qs_h>5, "5 to 10","Less than 5"))))))
# 
# quarantine_summary <- data.frame(matrix(nrow=7, ncol=2))
# 
# quarantine_summary[,1] <- c("0","Less than 5", "Greater than 5", "Greater than 10", "Greater than 20", "Greater than 30",
#                             "Greater than 50")
# 
# 
# quarantine_summary[1,2] <- length(which(results_all_quarantined$Qs_h==0))/length(results_all_quarantined$Qs_h)
# quarantine_summary[2,2] <- length(which(results_all_quarantined$Qs_h<=5))/length(results_all_quarantined$Qs_h)
# quarantine_summary[3,2] <- length(which(results_all_quarantined$Qs_h>5))/length(results_all_quarantined$Qs_h)
# quarantine_summary[4,2] <- length(which(results_all_quarantined$Qs_h>10))/length(results_all_quarantined$Qs_h)
# quarantine_summary[5,2] <- length(which(results_all_quarantined$Qs_h>20))/length(results_all_quarantined$Qs_h)
# quarantine_summary[6,2] <- length(which(results_all_quarantined$Qs_h>30))/length(results_all_quarantined$Qs_h)
# quarantine_summary[7,2] <- length(which(results_all_quarantined$Qs_h>50))/length(results_all_quarantined$Qs_h)
# 
# 
# quarantine_summary$X1 <- factor(quarantine_summary$X1,levels = c("0","Less than 5", "Greater than 5", "Greater than 10", "Greater than 20", "Greater than 30",
#                                                                  "Greater than 50"))
# 
# quarantinebeds <- ggplot(quarantine_summary, aes(x=X2, y=X1, fill=X1)) + geom_bar(stat="identity") + 
#   theme_classic() + theme(legend.position="none") + xlab("") + ylab("") + 
#   scale_x_continuous(labels = scales::percent) + ggtitle("Number of quarantine beds needed, % likelihood")


# Average quarantine beds

results_all_quarantined2 <- results_all_quarantined

results_all_quarantined2$time2 <- round(results_all_quarantined2$time)

results_all_quarantined3 <- results_all_quarantined2 %>%
  group_by(run,time2) %>%
  summarise_at(vars(Qs_h), list(maxQ = max))

# Time exceeding quarantine capacity

time_quar_past_cap <- percent(length(which(results_all_quarantined3$maxQ>quarantine_capacity_count))/length(results_all_quarantined3$maxQ))

# Likelihood exceeding quarantine capacity

results_all_quarantined_likelihood <- results_all_quarantined3 %>%
  group_by(run) %>%
  summarise_at(vars(maxQ), list(maxQQ = max))

likelihood_quar_past_cap <- percent(length(which(results_all_quarantined_likelihood$maxQQ>quarantine_capacity_count))/length(results_all_quarantined_likelihood$maxQQ))

results_all_quarantined4 <- results_all_quarantined3 %>%
  group_by(time2) %>%
  summarise_at(vars(maxQ), list(nmin=min, Q1=~quantile(., probs = 0.25), Q95l=~quantile(., probs = 0.05),
                                median=median, Q3=~quantile(., probs = 0.75),Q95u=~quantile(., probs = 0.95),
                                max=max))
ggplot(data=results_all_quarantined4, aes(x=time2, y=median)) + geom_line() +
  theme_classic() + theme(legend.position = "none") + 
  geom_ribbon(aes(ymin = Q95l, ymax = Q95u), alpha = 0.1) + xlab("Days") + ylab("Average number of quarantined students") +
  ggtitle("Average number of quarantined \nstudents by day, and 95% interval")

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
ggplot(data=results_all_infected4, aes(x=time2, y=median)) + geom_line() +
  theme_classic() + theme(legend.position = "none") + 
  geom_ribbon(aes(ymin = Q95l, ymax = Q95u), alpha = 0.1) + xlab("Days") + ylab("Average number of infectious students") +
  ggtitle("Average number of infectious \nstudents by day, and 95% interval")

# Average isolated students 

isolation_capacity_count= 5

results_all_diagnosed2 <- results_all_diagnosed 

results_all_diagnosed2 <- cbind(results_all_diagnosed, results_all_newlydiagnosed$Dx0_h)

results_all_diagnosed2$total <- results_all_diagnosed2$Dx_h + results_all_diagnosed2$`results_all_newlydiagnosed$Dx0_h`

results_all_diagnosed2$time2 <- round(results_all_diagnosed2$time)

results_all_diagnosed3 <- results_all_diagnosed2 %>%
  group_by(run,time2) %>%
  summarise_at(vars(total), list(maxD = max))

time_iso_past_cap <- percent(length(which(results_all_diagnosed3$maxD>isolation_capacity_count))/length(results_all_diagnosed3$maxD))

results_all_diagnosed_likelihood <- results_all_diagnosed3 %>%
  group_by(run) %>%
  summarise_at(vars(maxD), list(maxDD = max))

likelihood_iso_past_cap <- percent(length(which(results_all_diagnosed_likelihood$maxDD>isolation_capacity_count))/length(results_all_diagnosed_likelihood$maxDD))

results_all_diagnosed4 <- results_all_diagnosed3 %>%
  group_by(time2) %>%
  summarise_at(vars(maxD), list(nmin=min, Q1=~quantile(., probs = 0.25), Q95l=~quantile(., probs = 0.05),
                                median=median, Q3=~quantile(., probs = 0.75),Q95u=~quantile(., probs = 0.95),
                                max=max))

ggplot(data=results_all_diagnosed4, aes(x=time2, y=median)) + geom_line() +
  theme_classic() + theme(legend.position = "none") + 
  geom_ribbon(aes(ymin = Q95l, ymax = Q95u), alpha = 0.1) + xlab("Days") + ylab("Average number of isolated students") +
  ggtitle("Average number of isolated \nstudents by day, and 95% interval")

# Likelihood of an outbreak

# Likelihood of exceeding isolation capacity

isolation_capacity_count <- 5

isolation_capacity <- cbind(results_all_diagnosed, results_all_newlydiagnosed$Dx0_h)

isolation_capacity$total <- isolation_capacity$Dx_h + isolation_capacity$`results_all_newlydiagnosed$Dx0_h`

isolation_capacity_likelihood <- isolation_capacity %>%
  group_by(run) %>%
  summarise_at(vars(total), list(maxL = max))

likelihood_iso_past_cap <- percent(length(which(isolation_capacity_likelihood$maxL>isolation_capacity_count))/length(isolation_capacity_likelihood$maxL))

# Time exceeding isolation capacity

time_iso_past_cap <- percent(length(which(isolation_capacity$total>isolation_capacity_count))/length(isolation_capacity$total))


# Total infections

total_infections <- results_all_recovered
  
total_infections$allinfections <- total_infections$R_h + results_all_diagnosed$Dx_h + results_all_newlydiagnosed$Dx0_h +
                    results_all_infected$I_h + results_all_presymptomatic$P_h

total_infections2 <- total_infections[which(total_infections$time==100),]

median_infections <- median(total_infections2$allinfections)
l95_infections <- quantile(total_infections2$allinfections, probs = 0.05)
u95_infections <- quantile(total_infections2$allinfections, probs = 0.95)


# Max number of quarantined students  

results_all_quarantined2 <- results_all_quarantined %>% group_by(run) %>% summarise(Qs_h = max(Qs_h))
results_all_quarantined2$summary <- ifelse(results_all_quarantined2$Qs_h==0, "Zero",
                                           ifelse(results_all_quarantined2$Qs_h>50, "Greater than 50",
                                                  ifelse(results_all_quarantined2$Qs_h>30, "30 to 50",
                                                         ifelse(results_all_quarantined2$Qs_h>20, "20 to 30",
                                                                ifelse(results_all_quarantined2$Qs_h>10, "10 to 20",
                                                                       ifelse(results_all_quarantined2$Qs_h>5, "5 to 10","Less than 5"))))))


results_all_quarantined2$summary <- factor(results_all_quarantined2$summary,levels = c("Zero","Less than 5", "5 to 10", "10 to 20", "20 to 30", "30 to 50",
                                                                 "Greater than 50"))
ggplot(results_all_quarantined2, aes(y=summary,fill=summary)) + geom_bar()+ 
  theme_classic() + theme(legend.position="none") + xlab("") + ylab("") + 
  ggtitle("Max number of quarantine beds needed, % likelihood")


# Max number of isolated students  

results_all_diagnosed <- cbind(results_all_diagnosed, results_all_newlydiagnosed$Dx0_h)


results_all_diagnosed$total <- results_all_diagnosed$Dx_h + results_all_diagnosed$`results_all_newlydiagnosed$Dx0_h`


results_all_diagnosed2 <- results_all_diagnosed %>% group_by(run) %>% summarise(total = max(total))
results_all_diagnosed2$summary <- ifelse(results_all_diagnosed2$total==0, "Zero",
                                           ifelse(results_all_diagnosed2$total>50, "Greater than 50",
                                                  ifelse(results_all_diagnosed2$total>30, "30 to 50",
                                                         ifelse(results_all_diagnosed2$total>20, "20 to 30",
                                                                ifelse(results_all_diagnosed2$total>10, "10 to 20",
                                                                       ifelse(results_all_diagnosed2$total>5, "5 to 10","Less than 5"))))))


results_all_diagnosed2$summary <- factor(results_all_diagnosed2$summary,levels = c("Zero","Less than 5", "5 to 10", "10 to 20", "20 to 30", "30 to 50",
                                                                                       "Greater than 50"))
ggplot(results_all_diagnosed2, aes(y=summary,fill=summary)) + geom_bar()+ 
  theme_classic() + theme(legend.position="none") + xlab("") + ylab("") + 
  ggtitle("Max number of quarantine beds needed, % likelihood")


# graph of presymptomatic over time

ggplot(data=results_all_presymptomatic, aes(x=time, y=P_h, color=run)) + geom_line() +
  theme_classic() + theme(legend.position = "none") 

# graph of infected over time

ggplot(data=results_all_infected, aes(x=time, y=I_h, color=run)) + geom_line() +
  theme_classic() + theme(legend.position = "none") 

# graph of new diagnoses

ggplot(data=results_all_newlydiagnosed, aes(x=time, y=Dx0_h, color=run)) + geom_line()+
  theme_classic() + theme(legend.position = "none") 

# graph of diagnosed (and isolated) over time

ggplot(data=results_all_diagnosed, aes(x=time, y=Dx_h, color=run)) + geom_line()+
  theme_classic() + theme(legend.position = "none") 

#graph of recovered over time

ggplot(data=results_all_recovered, aes(x=time, y=R_h, color=run)) + geom_line()+
  theme_classic() + theme(legend.position = "none") 

# graph of all isolated over time

ggplot(data=results_all_quarantined, aes(x=time, y=Qs_h, color=run)) + geom_line()+
  theme_classic() + theme(legend.position = "none") 

mean(results_all_quarantined$Qs_h)
min(results_all_quarantined$Qs_h)
max(results_all_quarantined$Qs_h)





