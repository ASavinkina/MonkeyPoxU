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
    exograte = ifelse(input$exograte==0, 0, input$exograte/7)  # rate of exogenous shocks
    diagrate = input$diagrate # daily rate of diagnosis of infectious cases
    quarduration = 1/input$quarduration  # duration of quarantine for non-infected
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
    
    runs=100
    
    
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
    
    
    # graph of presymptomatic over time
    
    presymp_plot <- ggplot(data=results_all_presymptomatic, aes(x=time, y=P_h, color=run)) + geom_line() +
        theme_classic() + theme(legend.position = "none") + ylab("Number infected presymptomatic")
    
    
    # graph of infected over time
    
    infected_plot <- ggplot(data=results_all_infected, aes(x=time, y=I_h, color=run)) + geom_line() +
        theme_classic() + theme(legend.position = "none") + ylab("Number infected symptomatic")
    
    
    # graph of new diagnoses
    
    ggplot(data=results_all_newlydiagnosed, aes(x=time, y=Dx0_h, color=run)) + geom_line()+
        theme_classic() + theme(legend.position = "none") + ylab("Number newly diagnosed")
    
    
    # graph of diagnosed (and isolated) over time
    
    diagnosed_plot <- ggplot(data=results_all_diagnosed, aes(x=time, y=Dx_h, color=run)) + geom_line()+
        theme_classic() + theme(legend.position = "none") + ylab("Number diagnosed (cumulative)")
    
    
    #graph of recovered over time
    
    recovered_plot <- ggplot(data=results_all_recovered, aes(x=time, y=R_h, color=run)) + geom_line()+
        theme_classic() + theme(legend.position = "none") + ylab("Number recovered (cumulative)")
    
    
    # graph of all isolated over time
    
    quarantined_plot <- ggplot(data=results_all_quarantined, aes(x=time, y=Qs_h, color=run)) + geom_line()+
        theme_classic() + theme(legend.position = "none") + ylab("Number in quarantine (infected and not infected)")
    
    mean(results_all_quarantined$Qs_h)
    min(results_all_quarantined$Qs_h)
    max(results_all_quarantined$Qs_h)
    
    # output$QPlot <- renderPlot(quarantined_plot)
    # output$DPlot <- renderPlot(diagnosed_plot)
    # output$IPlot <- renderPlot(infected_plot)
    # output$RPlot <- renderPlot(recovered_plot)
    
    multiplot(quarantined_plot,diagnosed_plot,infected_plot,recovered_plot, cols=2)


})
    
})
