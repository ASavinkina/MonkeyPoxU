# Stochastic university monkeypox model
# Used for graphs for presentations and paper


# 11/1/22

# Import libraries

library(adaptivetau)
library(ggplot2)
library(tidyverse)
library(scales)
library(gridExtra)
library(RColorBrewer)


Model_runs = 1000

# Set parameters of interest

  # Number of initial infections at day 0 of 100 (high-risk)

initial_inf= 1 

  # Average number of exogenous infections over 100 days

exogshock =1

  # Proportion of cases detected and isolated

diagPerc = 0.2

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

R0_l = 0.8 #R0

#HR_prop = 0.1
#LR_prop = 1-HR_prop


# Set initial conditions for model

Proportion_Outbreak <-data.frame(matrix(NA, nrow = 100, ncol = 4))
colnames(Proportion_Outbreak) <- c("R0_h","PropHR","PropOutbreaks","MeanInf")

Proportion_Outbreak$R0_h <- rep(seq(0.8,3.5, by=0.3),10)
Proportion_Outbreak$PropHR <-  rep(seq(0.05,0.50, by=0.05),each=10)

for (k in 1:100) {


HR_prop <- Proportion_Outbreak[k,2]  
LR_prop <- 1- HR_prop

R0_h <- Proportion_Outbreak[k,1]  
#R0_l <- R0_h/2

init.values = c(
  S_h = round(Yale_undegrad_pop*HR_prop-initial_inf,0), P_h = 0, I_h = initial_inf, Dx0_h=0, Dx_h=0, 
  Vs_h=0, Vi_h=0, R_h=0,
  S_l = round(Yale_undegrad_pop*LR_prop,0), P_l = 0, I_l = 0, Dx0_l=0, Dx_l=0, 
  Vs_l=0, Vi_l=0, R_l=0
)

# Specify all transitions

transitions = list(
  # high risk population
  c(S_h = -1, P_h = +1), # movement from susceptible to presymptomatic, endogenous infection
  c(S_h = -1, P_h = +1), # movement from susceptible to presymptomatic, exogenous infection 
  c(S_h= -1, Vs_h= +1), #movement from susceptible to vaccinated susceptible
  c(P_h= -1, Vi_h= +1), #movement from presymptomatuc to vaccinated infected
  c(P_h = -1, I_h = +1), #movement from presymptomatic to infected 
  c(I_h = -1, Dx0_h = +1),  #movement from infected to newly diagnosed
  c(Dx0_h = -1, Dx_h = +1), #movement from newly diagnosed to diagnosed
  c(I_h = -1, R_h = +1), #movement from infected to recovered
  c(Dx_h= -1, R_h= +1), #movement from diagnosed to recovered
  
  #low risk population
  c(S_l = -1, P_l = +1), # movement from susceptible to presymptomatic, endogenous infection
  c(S_l= -1, Vs_l= +1), #movement from susceptible to vaccinated susceptible
  c(P_l= -1, Vi_l= +1), #movement from presymptomatuc to vaccinated infected
  c(P_l = -1, I_l = +1), #movement from presymptomatic to infected 
  c(I_l = -1, Dx0_l = +1),  #movement from infected to newly diagnosed
  c(Dx0_l = -1, Dx_l = +1), #movement from newly diagnosed to diagnosed
  c(I_l = -1, R_l = +1), #movement from infected to recovered
  c(Dx_l= -1, R_l= +1) #movement from diagnosed to recovered
  
)

# Rates for all transitions
# (In same order as "transitions")
RateF <- function(x, pars, times) {
  return(c(
    ifelse(x["S_h"]>1, pars$beta_l*x["S_h"]*x["I_l"]/(x["S_l"] + x["P_l"] + x["I_l"] + x["Vi_l"]+  x["Vs_l"] + x["R_l"]+ x["S_h"] + x["P_h"] + x["I_h"] + x["Vi_h"]+  x["Vs_h"]+ x["R_h"])+
      pars$beta_l*x["S_h"]*x["I_h"]/(x["S_l"] + x["P_l"] + x["I_l"] + x["Vi_l"]+  x["Vs_l"] + x["R_l"]+ x["S_h"] + x["P_h"] + x["I_h"] + x["Vi_h"]+  x["Vs_h"]+ x["R_h"])+
      pars$beta_h*x["S_h"]*x["I_h"]/(x["S_h"] + x["P_h"] + x["I_h"] + x["Vi_h"]+  x["Vs_h"]+ x["R_h"]),0),                                                                            
    pars$theta, # movement from susceptible to presymptomatic, exogenous infection
    ifelse(x["P_h"]>1,pars$iota*x["Dx0_h"]*(1-pars$attackrate)*pars$Vaccine_Efficacy,pars$iota*x["Dx0_h"]*pars$Vaccine_Efficacy),
    ifelse(x["P_h"]>1, pars$iota*x["Dx0_h"]*pars$attackrate*pars$Vaccine_Efficacy_PEP, 0), #movement from susceptible to vaccinated infected (based on number newly diagnosed)
    pars$gamma*x["P_h"], #movement from presymptomatic to infected (duration of incubation)
    pars$delta*x["I_h"], #movement from infected to diagnosed
    pars$tau*x["Dx0_h"], #movement from newly diagnosed to diagnosed
    pars$rho*x["I_h"], #movement from infected to recovered
    pars$omicron*x["Dx_h"], #movement from diagnosed to recovered 
    
    ifelse(x["S_l"]>1, pars$beta_l*x["S_l"]*x["I_l"]/(x["S_l"] + x["P_l"] + x["I_l"] + x["Vi_l"]+  x["Vs_l"] + x["R_l"]+ x["S_h"] + x["P_h"] + x["I_h"] + x["Vi_h"]+  x["Vs_h"]+ x["R_h"])+
    pars$beta_l*x["S_l"]*x["I_h"]/(x["S_l"] + x["P_l"] + x["I_l"] + x["Vi_l"]+  x["Vs_l"] + x["R_l"]+ x["S_h"] + x["P_h"] + x["I_h"] + x["Vi_h"]+  x["Vs_h"]+ x["R_h"]),0),
    ifelse(x["P_l"]>1,pars$iota*x["Dx0_l"]*(1-pars$attackrate)*pars$Vaccine_Efficacy,pars$iota*x["Dx0_l"]*pars$Vaccine_Efficacy),
    ifelse(x["P_l"]>1, pars$iota*x["Dx0_l"]*pars$attackrate*pars$Vaccine_Efficacy_PEP, 0), #movement from susceptible to vaccinated infected (based on number newly diagnosed)
    pars$gamma*x["P_l"], #movement from presymptomatic to infected (duration of incubation)
    pars$delta*x["I_l"], #movement from infected to diagnosed
    pars$tau*x["Dx0_l"], #movement from newly diagnosed to diagnosed
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
  omega = quarrate, #length of vaccination for susceptible
  attackrate=0.2, 
  theta= exograte, #rate of exogenous shocks
  omicron = isorate,
  Vaccine_Efficacy=Vaccine_Efficacy,
  Vaccine_Efficacy_PEP = Vaccine_Efficacy_PEP
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

runs=Model_runs

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
  
  # results_vaccinated_i <- results[,c(1,7,15)]
  # results_vaccinated_i[,2] <- results_vaccinated_i[,2]/Vaccine_Efficacy + results[,8]/Vaccine_Efficacy_PEP 
  # results_vaccinated_i[,3] <- results_vaccinated_i[,3]/Vaccine_Efficacy + results[,16]/Vaccine_Efficacy_PEP 
  # results_vaccinated_i$run <- paste0(i)
  # 
  # results_all_vaccinated <- rbind(results_all_vaccinated, results_vaccinated_i)
  
  
}


total_infections <- results_all_recovered

total_infections$R <- results_all_recovered$R_h + results_all_recovered$R_l

total_infections$allinfections <- total_infections$R + results_all_diagnosed$Dx_h + 
  results_all_diagnosed$Dx_l + 
  results_all_newlydiagnosed$Dx0_h +
  results_all_newlydiagnosed$Dx0_l +
  results_all_infected$I_h +
  results_all_infected$I_l+ 
  results_all_presymptomatic$P_h +
  results_all_presymptomatic$P_l

#   Calculate total number of infectious students 

total_infections2 <- total_infections[which(total_infections$time==100),]

#   Calculate number of expected index cases

expected_cases <- initial_inf+exogshock

#   Calculate proportion of times there are more cases than the expected number

additional_cases <- length(total_infections2$allinfections[total_infections2$allinfections>expected_cases])/Model_runs

Proportion_Outbreak[k,3]  <- additional_cases

#   Calculate mean number of additional cases given secondary case

additional_cases_mean <- mean(total_infections2$allinfections[total_infections2$allinfections>expected_cases])

Proportion_Outbreak[k,4]  <- additional_cases_mean

}


write.csv(Proportion_Outbreak, file="TwoWaySA_HR_R0_11.12.22.csv")



Prop_Outbreak_Graph <-ggplot(data=Proportion_Outbreak, aes(x=Proportion_Outbreak$R0_h, y=Proportion_Outbreak$PropHR, fill=Proportion_Outbreak$PropOutbreaks)) + 
  geom_tile() + theme_classic() + xlab("R0, high risk group") + ylab("Proportion of population that is high-risk") 

Prop_Outbreak_Graph + theme(legend.title = element_blank(),text=element_text(size = 15))

Blues <- brewer.pal(n = 7, name = 'Blues')


Proportion_Outbreak$Mean2 <- round(Proportion_Outbreak$MeanInf,0)
Proportion_Outbreak$Mean3 <- ifelse(Proportion_Outbreak$Mean2<20,"<20",
                                    ifelse(Proportion_Outbreak$Mean2>=20 & Proportion_Outbreak$Mean2<30, "20-30",
                                           ifelse(Proportion_Outbreak$Mean2>=30 & Proportion_Outbreak$Mean2<50, "30-50",
                                                  ifelse(Proportion_Outbreak$Mean2>=50 & Proportion_Outbreak$Mean2<100, "50-100",
                                                         ifelse(Proportion_Outbreak$Mean2>=100 & Proportion_Outbreak$Mean2<150, "100-150",
                                                                ifelse(Proportion_Outbreak$Mean2>=150 & Proportion_Outbreak$Mean2<300, "150-300",
                                                                       ifelse(Proportion_Outbreak$Mean2>=300, "300+",NA)))))))


Proportion_Outbreak$Mean4 <- factor(Proportion_Outbreak$Mean3, levels = c("<20", "20-30", "30-50","50-100","100-150","150-300","300+"))


Cases_Outbreak <- ggplot(data=Proportion_Outbreak, aes(x=Proportion_Outbreak$R0_h, y=Proportion_Outbreak$PropHR, fill=Proportion_Outbreak$Mean4)) +
  geom_tile() + scale_fill_manual(values=Blues) + theme_classic() + xlab("R0, high risk group") + ylab("Proportion of population that is high-risk") 

Cases_Outbreak + theme(legend.title = element_blank(),text = element_text(size = 15)) 

















