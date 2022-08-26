#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(adaptivetau)
library(ggplot2)
library(tidyverse)
library(shinydashboard)
library(scales)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
    
    output$I_plot <- renderPlot({
        
        ############# Multiple plot function ###################
        ## http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/ ##
        
        multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
            library(grid)
            
            # Make a list from the ... arguments and plotlist
            plots <- c(list(...), plotlist)
            
            numPlots = length(plots)
            
            # If layout is NULL, then use 'cols' to determine layout
            if (is.null(layout)) {
                # Make the panel
                # ncol: Number of columns of plots
                # nrow: Number of rows needed, calculated from # of cols
                layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                                 ncol = cols, nrow = ceiling(numPlots/cols))
            }
            
            if (numPlots==1) {
                print(plots[[1]])
                
            } else {
                # Set up the page
                grid.newpage()
                pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
                
                # Make each plot, in the correct location
                for (i in 1:numPlots) {
                    # Get the i,j matrix positions of the regions that contain this subplot
                    matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
                    
                    print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                                    layout.pos.col = matchidx$col))
                }
            } 
        }
    
    
    # Intiial conditions
    
    Yale_undegrad_pop = input$Yale_undegrad_pop #input$Yale_undegrad_pop # total population
    initial_inf= input$initial_inf #input$initial_inf # initial infections

    exogshock =1 #size of exogenous shock. 0 for no exogenous shocks, >1 for superspreader events
    exograte = ifelse(input$exograte==0, 0, input$exograte/30)  # rate of exogenous shocks
    recoveryrate = 1/21
    diagrate = ifelse(input$diagrate>0.95, 1, (input$diagrate*recoveryrate)/(1-input$diagrate)) # daily rate of diagnosis of infectious cases
    isorate = 1/input$isoduration #duration of isolation
    quarrate = 1/input$quarduration  # duration of quarantine 
    studentcontacts = input$studentcontacts # students quarantined per diagnosed case
    R0_h= input$R0_h
        
        # Yale_undegrad_pop = 6500 # total population
        # initial_inf = 10 # initial infections
        # exogshock =0 #size of exogenous shock. 0 for no exogenous shocks, >1 for superspreader events
        # exograte = 1/30 # rate of exogenous shocks
        # diagrate = 1/10 # daily rate of diagnosis of infectious cases
        # quarduration = 1/14 # duration of quarantine for non-infected
        # studentcontacts = 5 # students quarantined per diagnosed case
        # R0_h = 1.5
    
    
        init.values = c(
            S_h = Yale_undegrad_pop-initial_inf, P_h = 0, I_h = initial_inf, Dx0_h=0, Dx_h=0, 
            Qs_h=0, Qi_h=0, R_h=0
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
            pars$omega*x["Qi_h"], #movement from quarantined infected to infected
            pars$gamma*x["P_h"], #movement from presymptomatic to infected (duration of incubation)
            pars$delta*x["I_h"], #movement from infected to diagnosed
            pars$tau*x["Dx0_h"], #movement from newly diagnosed to diagnosed
            pars$mu*x["Qi_h"], #movement from quarantined infected to diagnosed
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
        iota= studentcontacts, # students quarantined per diagnosed case
        mu = 1/7.6, #incubation period for those in quarantine
        omega = quarrate, #length of quarantine for susceptible
        attackrate=0.2, 
        theta= exograte, #rate of exogenous shocks
        omicron = isorate
    )
    
    
    # Running stochastic model
    
    
    # results = as.data.frame(ssa.adaptivetau(init.values, 
    #                                         transitions, 
    #                                         RateF, 
    #                                         pars, 
    #                                         tf=100))
    
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
    results_all_quarantined<- data.frame(matrix(0, nrow=0, ncol=3))
    colnames(results_all_quarantined) <- c("time", "Dx0_h", "run")
    
    runs=100
    
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
        
        results_quarantined_i <- results[,c(1,7)] 
        results_quarantined_i[,2] <- results_quarantined_i[,2] + results[,8] 
        results_quarantined_i$run <- paste0(i)
        
        results_all_quarantined <- rbind(results_all_quarantined, results_quarantined_i)
        
        
    }
    

    # Total infections
    
    total_infections <- results_all_recovered
    
    total_infections$allinfections <- total_infections$R_h + results_all_diagnosed$Dx_h + results_all_newlydiagnosed$Dx0_h +
        results_all_infected$I_h + results_all_presymptomatic$P_h
    
    total_infections2 <- total_infections[which(total_infections$time==100),]
    
    output$maxinfectionsplot <- renderPlot(hist(total_infections2$allinfections,
                                             main="",xlab="Cumulative number of cases at 100 days",
                                             ylab="Percent likelihood"),
                                            width=300, height=300)
    
    # maxinfectionsplot <- ggplot(total_infections2, aes(x=allinfections)) + 
    #     geom_histogram(fill='grey', color='black') +
    #     theme_classic() + 
    #     # scale_x_continuous(breaks=c(0,500), limits = c(0,500)) +
    #     xlab('Cumulative infections in 100 days') + ylab('Likelihood')

    median_infections <- median(total_infections2$allinfections)
    nonewinfections <- percent(length(which(total_infections2$allinfections==input$initial_inf))/length(total_infections2$allinfections))
    
    min_infections <- min(total_infections2$allinfections)
    max_infections <- max(total_infections2$allinfections)
    
    # Average quarantine beds
    
    quarantine_capacity_count <- input$quarcapacity
    
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
    
    # maxquarplot <- ggplot(results_all_quarantined_likelihood, aes(x=maxQQ)) + 
    #     geom_histogram(fill='grey', color='black') +
    #     theme_classic() + 
    #     # scale_x_continuous(breaks=c(0,100), limits = c(0,100)) +
    #     xlab('Maximum number of quarantined students') + ylab('Likelihood')
    
    output$maxquarplot <- renderPlot(hist(results_all_quarantined_likelihood$maxQQ,
                                         main="",xlab="Max number of students in quarantine",
                                         ylab="Percent likelihood"), height=300, width=300)
    
             
    
    likelihood_quar_past_cap <- percent(length(which(results_all_quarantined_likelihood$maxQQ>quarantine_capacity_count))/length(results_all_quarantined_likelihood$maxQQ))
    
    results_all_quarantined4 <- results_all_quarantined3 %>%
        group_by(time2) %>%
        summarise_at(vars(maxQ), list(nmin=min, Q1=~quantile(., probs = 0.25), Q95l=~quantile(., probs = 0.05),
                                      median=median, Q3=~quantile(., probs = 0.75),Q95u=~quantile(., probs = 0.95),
                                      max=max))
    
    avg_quarantine_plot <- ggplot(data=results_all_quarantined4, aes(x=time2, y=median)) + geom_line() +
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
    avg_infectious_plot <- ggplot(data=results_all_infected4, aes(x=time2, y=median)) + geom_line() +
        theme_classic() + theme(legend.position = "none") + 
        geom_ribbon(aes(ymin = Q95l, ymax = Q95u), alpha = 0.1) + xlab("Days") + ylab("Average number of infectious students") +
        ggtitle("Average number of infectious \nstudents by day, and 95% interval")
    
    # Average isolated students 
    
    isolation_capacity_count <- input$isocapacity
    
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
    
    output$maxisoplot <- renderPlot(hist(results_all_diagnosed_likelihood$maxDD,
                                         main="",xlab="Max number of students in isolation",
                                         ylab="Percent likelihood"), height=300, width=300)


    # maxisoplot <- ggplot(results_all_diagnosed_likelihood, aes(x=maxDD)) + 
    #     geom_histogram(fill='grey', color='black') +
    #     theme_classic() + 
    #     # scale_x_continuous(breaks=c(0,50), limits = c(0,50)) +
    #     xlab('Maximum number of isolated students') + ylab('Likelihood')
    
    likelihood_iso_past_cap <- percent(length(which(results_all_diagnosed_likelihood$maxDD>isolation_capacity_count))/length(results_all_diagnosed_likelihood$maxDD))
    
    results_all_diagnosed4 <- results_all_diagnosed3 %>%
        group_by(time2) %>%
        summarise_at(vars(maxD), list(nmin=min, Q1=~quantile(., probs = 0.25), Q95l=~quantile(., probs = 0.05),
                                      median=median, Q3=~quantile(., probs = 0.75),Q95u=~quantile(., probs = 0.95),
                                      max=max))
    
    avg_isolated_plot <- ggplot(data=results_all_diagnosed4, aes(x=time2, y=median)) + geom_line() +
        theme_classic() + theme(legend.position = "none") + 
        geom_ribbon(aes(ymin = Q95l, ymax = Q95u), alpha = 0.1) + xlab("Days") + ylab("Average number of isolated students") +
        ggtitle("Average number of isolated \nstudents by day, and 95% interval")
    
    
    # graph of presymptomatic over time
    
    presymp_plot <- ggplot(data=results_all_presymptomatic, aes(x=time, y=P_h, color=run)) + geom_line() +
        theme_classic() + theme(legend.position = "none") + ylab("Number infected presymptomatic") +
        xlab("Days")
    
    
    # graph of infected over time
    
    infected_plot <- ggplot(data=results_all_infected, aes(x=time, y=I_h, color=run)) + geom_line() +
        theme_classic() + theme(legend.position = "none") + ylab("Number infected symptomatic") +
        xlab("Days")
    
    
    # graph of new diagnoses
    
    ggplot(data=results_all_newlydiagnosed, aes(x=time, y=Dx0_h, color=run)) + geom_line()+
        theme_classic() + theme(legend.position = "none") + ylab("Number newly diagnosed") +
        xlab("Days")
    
    
    # graph of diagnosed (and isolated) over time
    
    diagnosed_plot <- ggplot(data=results_all_diagnosed, aes(x=time, y=Dx_h, color=run)) + geom_line()+
        theme_classic() + theme(legend.position = "none") + ylab("Number in isolation (diagnosed+)") +
        xlab("Days")
    
    
    #graph of recovered over time
    
    recovered_plot <- ggplot(data=results_all_recovered, aes(x=time, y=R_h, color=run)) + geom_line()+
        theme_classic() + theme(legend.position = "none") + ylab("Number recovered (cumulative)") +
        xlab("Days")
    
    
    # graph of all isolated over time
    
    quarantined_plot <- ggplot(data=results_all_quarantined, aes(x=time, y=Qs_h, color=run)) + geom_line()+
        theme_classic() + theme(legend.position = "none") + ylab("Number in quarantine (+ and -)") +
        xlab("Days")
    
    mean(results_all_quarantined$Qs_h)
    min(results_all_quarantined$Qs_h)
    max(results_all_quarantined$Qs_h)
    
    
    output$DPlot <- renderPlot(avg_isolated_plot, height=300, width=300)
    output$QPlot <- renderPlot(avg_quarantine_plot, height=300, width=300)
    output$IPlot <- renderPlot(avg_infectious_plot, height=300, width=300)
    output$D1Plot <- renderPlot(diagnosed_plot, height=300, width=300)
    output$Q1Plot <- renderPlot(quarantined_plot, height=300, width=300)
    output$I1Plot <- renderPlot(infected_plot, height=300, width=300)
    output$R1Plot <- renderPlot(recovered_plot, height=300, width=300)
    # output$maxquarplot <- renderPlot(maxquarplot, height=300, width=300)
    # output$maxinfectionsplot <- renderPlot(maxinfectionsplot, height=300, width=300)
    # output$maxisoplot <- renderPlot(maxisoplot, height=300, width=300)
    #output$total_inf_hist <- renderPlot(total_inf_hist, height=200, width=200)
    #output$total_iso_hist
    
    #multiplot1(avg_isolated_plot, avg_quarantine_plot, avg_infectious_plot,quarantined_plot,diagnosed_plot,infected_plot,recovered_plot, cols=2)
    
    #multiplot2(avg_isolated_plot, avg_quarantine_plot, avg_infectious_plot, cols=2)
    
    output$isocaplikelihood <- renderValueBox({
        valueBox(likelihood_iso_past_cap, "Likelihood of exceeding isolation capacity")
    })
    
    output$isocaptime <- renderValueBox({
        valueBox(time_iso_past_cap, "Percent of time isolation capacity exceeded")
    })
    
    output$quarcaplikelihood <- renderValueBox({
        valueBox(likelihood_quar_past_cap, "Likelihood of exceeding quarantine capacity")
    })
    
    output$quarcaptime <- renderValueBox({
        valueBox(time_quar_past_cap, "Percent of time quarantine capacity exceeded")
    })
    
    output$medianinfections <- renderValueBox({
        valueBox(paste(median_infections,"[",min_infections,",",max_infections,"]"), "Median [min,max] number infections in 100 days")
    })
    

    output$nonewinfections <- renderValueBox({
        valueBox(nonewinfections, "Likelihood of no infections post-incoming infections at start of semester")
    })

    
    
})
    
     })
