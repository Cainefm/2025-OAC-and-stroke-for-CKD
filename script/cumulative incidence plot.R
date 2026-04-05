K=60
pool.logistic.plot <-speedglm(out ~ expo + I(time) + I(timesqr) + Sex + age_indx + time_af_index + 
                                chadsvas_score + dx.cbd + dx.copd + dx.dementia + dx.stroke_embo + 
                                dx.ra + dx.dm_com0 + dx.htn + antiplatelet + I(time*expo) + I(timesqr*expo),
                              family=binomial(link="logit"),data=person_trial)

# plot --------------------------------------------------------------------
{
# Create dataset with all time points for each individual under each treatment level
treat0 <- expandRows(person_trial[which(person_trial$time==0),], count=K, count.is.col=F) 
treat0$time <- rep(seq(0, 59), nrow(person_trial[which(person_trial$time==0),]))
treat0$timesqr <- treat0$time^2

# Under flu 
treat0$expo <- as.numeric(as.character(treat0$expo))
treat0$expo <- 0
# Under covid 
treat1 <- treat0
treat1$expo <- 1


# Extract predicted values from pooled logistic regression model for each person-time row
# Predicted values correspond to discrete-time hazards
treat0$p.event0 <- predict(pool.logistic.plot, treat0, type="response")
treat1$p.event1 <- predict(pool.logistic.plot, treat1, type="response")

# Obtain predicted survival probabilities from discrete-time hazards
treat0.surv <- treat0 %>% group_by(id) %>% mutate(surv0 = cumprod(1 - p.event0)) %>% ungroup()
treat1.surv <- treat1 %>% group_by(id) %>% mutate(surv1 = cumprod(1 - p.event1)) %>% ungroup()

# Estimate risks from survival probabilities
# Risk = 1 - S(t)
treat0.surv$risk0 <- 1 - treat0.surv$surv0
treat1.surv$risk1 <- 1 - treat1.surv$surv1

# Get the mean in each treatment group at each time point from 0 to 29 (30 time points in total)
treat0.surv$expo <- as.numeric(as.character(treat0.surv$expo))
risk0 <- aggregate(treat0.surv[c("expo", "time", "risk0")], by=list(treat0.surv$time), FUN=mean)[c("expo", "time", "risk0")] 
treat1.surv$expo <- as.numeric(as.character(treat1.surv$expo))
risk1 <- aggregate(treat1.surv[c("expo", "time", "risk1")], by=list(treat1.surv$time), FUN=mean)[c("expo", "time", "risk1")] 

# Prepare data
graph.pred <- merge(risk0, risk1, by=c("time"))
# Edit data frame to reflect that risks are estimated at the END of each interval
graph.pred$time_0 <- graph.pred$time + 1
zero <- data.frame(cbind(0,0,0,1,0,0))
zero <- setNames(zero,names(graph.pred))
graph <- rbind(zero, graph.pred)

### Use pooled logistic regression estimates to compute causal estimates ###

# 30-days risk in flu
risk0.plr <- graph$risk0[which(graph$time_0==K-1)]
risk0.plr

# 30-days risk in COVID
risk1.plr <- graph$risk1[which(graph$time_0==K-1)]
risk1.plr

# 30-days risk difference
rd.plr <- risk1.plr - risk0.plr
rd.plr

# 30-days risk ratio
rr.plr <- risk1.plr / risk0.plr
rr.plr

ggplot(graph, 
       aes(x=time_0, y=risk)) + # set x and y axes
  geom_line(aes(y = risk1, # create line for COVID group
                color = "OAC"),
            size = 1.5) + 
  geom_line(aes(y = risk0, # create line for no FLU group
                color = "NonUsers"),
            size = 1.5) +
  xlab("Years") + # label x axis
  scale_x_continuous(limits = c(0, 60), # format x axis
                     breaks=c(0, 12, 24, 36, 48, 60),
                     labels =c(0:5) ) +
  ylab("Cumulative Incidence (%)") + # label y axis
  scale_y_continuous(limits=c(0, 0.6), # format y axis
                     breaks=c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6),
                     labels=c("0%","10%",  "20%", "30%", "40%","50%","60%")) +
  theme_minimal()+ # set plot theme elements
scale_color_manual(name="OAC users",
                   values=c("#386cb0","#f87f01"), # set colors
                   breaks=c('NonUsers', 'OAC'),
                   labels=c("None users","users")) +
theme(plot.title = element_text(face = "bold",
                                size = rel(1.2), hjust = 0.5),
      text = element_text(),
      panel.background = element_rect(colour = NA),
      plot.background = element_rect(colour = NA),
      # panel.border = element_rect(colour = NA),
      axis.title = element_text(face = "bold",size = rel(1)),
      axis.title.y = element_text(angle=90,vjust =2),
      axis.title.x = element_text(vjust = -0.2),
      axis.text = element_text(), 
      axis.line = element_line(colour="black"),
      axis.ticks = element_line(),
      panel.grid.major = element_line(colour="#f0f0f0"),
      panel.grid.minor = element_blank(),
      legend.key = element_rect(colour = NA),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.key.size= unit(0.2, "cm"),
      legend.title = element_text(face="italic"),
      plot.margin=unit(c(10,5,5,5),"mm"),
      strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
      strip.text = element_text(face="bold")
) -> plot_main
}


# male
{
  data4plot <- person_trial[Sex=="M"]
  pool.logistic.plot <-speedglm(out ~ expo + I(time) + I(timesqr) + age_indx + time_af_index + 
                                  chadsvas_score + dx.cbd + dx.copd + dx.dementia + dx.stroke_embo + 
                                  dx.ra + dx.dm_com0 + dx.htn + antiplatelet + I(time*expo) + I(timesqr*expo),
                                family=binomial(link="logit"),data=data4plot, maxit = 100 )
  # Create dataset with all time points for each individual under each treatment level
  treat0 <- expandRows(data4plot[which(data4plot$time==0),], count=K, count.is.col=F) 
  treat0$time <- rep(seq(0, 59), nrow(data4plot[which(data4plot$time==0),])) 
  treat0$timesqr <- treat0$time^2
  
  # Under flu 
  treat0$expo <- as.numeric(as.character(treat0$expo))
  treat0$expo <- 0
  # Under covid 
  treat1 <- treat0
  treat1$expo <- 1
  
  
  # Extract predicted values from pooled logistic regression model for each person-time row
  # Predicted values correspond to discrete-time hazards
  treat0$p.event0 <- predict(pool.logistic.plot, treat0, type="response")
  treat1$p.event1 <- predict(pool.logistic.plot, treat1, type="response")
  
  # Obtain predicted survival probabilities from discrete-time hazards
  treat0.surv <- treat0 %>% group_by(id) %>% mutate(surv0 = cumprod(1 - p.event0)) %>% ungroup()
  treat1.surv <- treat1 %>% group_by(id) %>% mutate(surv1 = cumprod(1 - p.event1)) %>% ungroup()
  
  # Estimate risks from survival probabilities
  # Risk = 1 - S(t)
  treat0.surv$risk0 <- 1 - treat0.surv$surv0
  treat1.surv$risk1 <- 1 - treat1.surv$surv1
  
  # Get the mean in each treatment group at each time point from 0 to 29 (30 time points in total)
  treat0.surv$expo <- as.numeric(as.character(treat0.surv$expo))
  risk0 <- aggregate(treat0.surv[c("expo", "time", "risk0")], by=list(treat0.surv$time), FUN=mean)[c("expo", "time", "risk0")] 
  treat1.surv$expo <- as.numeric(as.character(treat1.surv$expo))
  risk1 <- aggregate(treat1.surv[c("expo", "time", "risk1")], by=list(treat1.surv$time), FUN=mean)[c("expo", "time", "risk1")] 
  
  # Prepare data
  graph.pred <- merge(risk0, risk1, by=c("time"))
  # Edit data frame to reflect that risks are estimated at the END of each interval
  graph.pred$time_0 <- graph.pred$time + 1
  zero <- data.frame(cbind(0,0,0,1,0,0))
  zero <- setNames(zero,names(graph.pred))
  graph <- rbind(zero, graph.pred)
  
  ### Use pooled logistic regression estimates to compute causal estimates ###
  
  # 30-days risk in flu
  risk0.plr <- graph$risk0[which(graph$time_0==K-1)]
  risk0.plr
  
  # 30-days risk in COVID
  risk1.plr <- graph$risk1[which(graph$time_0==K-1)]
  risk1.plr
  
  # 30-days risk difference
  rd.plr <- risk1.plr - risk0.plr
  rd.plr
  
  # 30-days risk ratio
  rr.plr <- risk1.plr / risk0.plr
  rr.plr
  
  ggplot(graph, 
         aes(x=time_0, y=risk)) + # set x and y axes
    geom_line(aes(y = risk1, # create line for COVID group
                  color = "OAC"),
              size = 1.5) + 
    geom_line(aes(y = risk0, # create line for no FLU group
                  color = "NonUsers"),
              size = 1.5) +
    xlab("Years") + # label x axis
    scale_x_continuous(limits = c(0, 60), # format x axis
                       breaks=c(0, 12, 24, 36, 48, 60),
                       labels =c(0:5) ) +
    ylab("Cumulative Incidence (%)") + # label y axis
    scale_y_continuous(limits=c(0, 0.6), # format y axis
                       breaks=c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6),
                       labels=c("0%","10%",  "20%", "30%", "40%","50%","60%")) +
    theme_minimal()+ # set plot theme elements
    scale_color_manual(name="OAC users",
                       values=c("#386cb0","#f87f01"), # set colors
                       breaks=c('NonUsers', 'OAC'),
                       labels=c("None users","users")) +
    theme(plot.title = element_text(face = "bold",
                                    size = rel(1.2), hjust = 0.5),
          text = element_text(),
          panel.background = element_rect(colour = NA),
          plot.background = element_rect(colour = NA),
          # panel.border = element_rect(colour = NA),
          axis.title = element_text(face = "bold",size = rel(1)),
          axis.title.y = element_text(angle=90,vjust =2),
          axis.title.x = element_text(vjust = -0.2),
          axis.text = element_text(), 
          axis.line = element_line(colour="black"),
          axis.ticks = element_line(),
          panel.grid.major = element_line(colour="#f0f0f0"),
          panel.grid.minor = element_blank(),
          legend.key = element_rect(colour = NA),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.key.size= unit(0.2, "cm"),
          legend.title = element_text(face="italic"),
          plot.margin=unit(c(10,5,5,5),"mm"),
          strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
          strip.text = element_text(face="bold")
    ) -> plot_male
}

# male
{
  data4plot <- person_trial[Sex=="F"]
  pool.logistic.plot <-speedglm(out ~ expo + I(time) + I(timesqr) + age_indx + time_af_index + 
                                  chadsvas_score + dx.cbd + dx.copd + dx.dementia + dx.stroke_embo + 
                                  dx.ra + dx.dm_com0 + dx.htn + antiplatelet + I(time*expo) + I(timesqr*expo),
                                family=binomial(link="logit"),data=data4plot, maxit = 100 )
  # Create dataset with all time points for each individual under each treatment level
  treat0 <- expandRows(data4plot[which(data4plot$time==0),], count=K, count.is.col=F) 
  treat0$time <- rep(seq(0, 59), nrow(data4plot[which(data4plot$time==0),])) 
  treat0$timesqr <- treat0$time^2
  
  # Under flu 
  treat0$expo <- as.numeric(as.character(treat0$expo))
  treat0$expo <- 0
  # Under covid 
  treat1 <- treat0
  treat1$expo <- 1
  
  
  # Extract predicted values from pooled logistic regression model for each person-time row
  # Predicted values correspond to discrete-time hazards
  treat0$p.event0 <- predict(pool.logistic.plot, treat0, type="response")
  treat1$p.event1 <- predict(pool.logistic.plot, treat1, type="response")
  
  # Obtain predicted survival probabilities from discrete-time hazards
  treat0.surv <- treat0 %>% group_by(id) %>% mutate(surv0 = cumprod(1 - p.event0)) %>% ungroup()
  treat1.surv <- treat1 %>% group_by(id) %>% mutate(surv1 = cumprod(1 - p.event1)) %>% ungroup()
  
  # Estimate risks from survival probabilities
  # Risk = 1 - S(t)
  treat0.surv$risk0 <- 1 - treat0.surv$surv0
  treat1.surv$risk1 <- 1 - treat1.surv$surv1
  
  # Get the mean in each treatment group at each time point from 0 to 29 (30 time points in total)
  treat0.surv$expo <- as.numeric(as.character(treat0.surv$expo))
  risk0 <- aggregate(treat0.surv[c("expo", "time", "risk0")], by=list(treat0.surv$time), FUN=mean)[c("expo", "time", "risk0")] 
  treat1.surv$expo <- as.numeric(as.character(treat1.surv$expo))
  risk1 <- aggregate(treat1.surv[c("expo", "time", "risk1")], by=list(treat1.surv$time), FUN=mean)[c("expo", "time", "risk1")] 
  
  # Prepare data
  graph.pred <- merge(risk0, risk1, by=c("time"))
  # Edit data frame to reflect that risks are estimated at the END of each interval
  graph.pred$time_0 <- graph.pred$time + 1
  zero <- data.frame(cbind(0,0,0,1,0,0))
  zero <- setNames(zero,names(graph.pred))
  graph <- rbind(zero, graph.pred)
  
  ### Use pooled logistic regression estimates to compute causal estimates ###
  
  # 30-days risk in flu
  risk0.plr <- graph$risk0[which(graph$time_0==K-1)]
  risk0.plr
  
  # 30-days risk in COVID
  risk1.plr <- graph$risk1[which(graph$time_0==K-1)]
  risk1.plr
  
  # 30-days risk difference
  rd.plr <- risk1.plr - risk0.plr
  rd.plr
  
  # 30-days risk ratio
  rr.plr <- risk1.plr / risk0.plr
  rr.plr
  
  ggplot(graph, 
         aes(x=time_0, y=risk)) + # set x and y axes
    geom_line(aes(y = risk1, # create line for COVID group
                  color = "OAC"),
              size = 1.5) + 
    geom_line(aes(y = risk0, # create line for no FLU group
                  color = "NonUsers"),
              size = 1.5) +
    xlab("Years") + # label x axis
    scale_x_continuous(limits = c(0, 60), # format x axis
                       breaks=c(0, 12, 24, 36, 48, 60),
                       labels =c(0:5) ) +
    ylab("Cumulative Incidence (%)") + # label y axis
    scale_y_continuous(limits=c(0, 0.6), # format y axis
                       breaks=c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6),
                       labels=c("0%","10%",  "20%", "30%", "40%","50%","60%")) +
    theme_minimal()+ # set plot theme elements
    scale_color_manual(name="OAC users",
                       values=c("#386cb0","#f87f01"), # set colors
                       breaks=c('NonUsers', 'OAC'),
                       labels=c("None users","users")) +
    theme(plot.title = element_text(face = "bold",
                                    size = rel(1.2), hjust = 0.5),
          text = element_text(),
          panel.background = element_rect(colour = NA),
          plot.background = element_rect(colour = NA),
          # panel.border = element_rect(colour = NA),
          axis.title = element_text(face = "bold",size = rel(1)),
          axis.title.y = element_text(angle=90,vjust =2),
          axis.title.x = element_text(vjust = -0.2),
          axis.text = element_text(), 
          axis.line = element_line(colour="black"),
          axis.ticks = element_line(),
          panel.grid.major = element_line(colour="#f0f0f0"),
          panel.grid.minor = element_blank(),
          legend.key = element_rect(colour = NA),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.key.size= unit(0.2, "cm"),
          legend.title = element_text(face="italic"),
          plot.margin=unit(c(10,5,5,5),"mm"),
          strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
          strip.text = element_text(face="bold")
    ) -> plot_female
}

person_trial[,age_indx:=as.numeric(age_indx)]
# recurrent
# {
#   data4plot <- person_trial[dx.stroke_embo==1]
#   pool.logistic.plot <-speedglm(out ~ expo + I(time) + I(timesqr) + Sex + age_indx + time_af_index + 
#                                   chadsvas_score + dx.cbd + dx.copd + dx.dementia + 
#                                   dx.ra + dx.dm_com0 + dx.htn + antiplatelet + I(time*expo) + I(timesqr*expo),
#                                 family=binomial(link="logit"),data=data4plot, maxit = 100 )
#   # Create dataset with all time points for each individual under each treatment level
#   treat0 <- expandRows(data4plot[which(data4plot$time==0),], count=K, count.is.col=F) 
#   treat0$time <- rep(seq(0, 59), nrow(data4plot[which(data4plot$time==0),])) 
#   treat0$timesqr <- treat0$time^2
#   
#   # Under flu 
#   treat0$expo <- as.numeric(as.character(treat0$expo))
#   treat0$expo <- 0
#   # Under covid 
#   treat1 <- treat0
#   treat1$expo <- 1
#   
#   
#   # Extract predicted values from pooled logistic regression model for each person-time row
#   # Predicted values correspond to discrete-time hazards
#   treat0$p.event0 <- predict(pool.logistic.plot, treat0, type="response")
#   treat1$p.event1 <- predict(pool.logistic.plot, treat1, type="response")
#   
#   # Obtain predicted survival probabilities from discrete-time hazards
#   treat0.surv <- treat0 %>% group_by(id) %>% mutate(surv0 = cumprod(1 - p.event0)) %>% ungroup()
#   treat1.surv <- treat1 %>% group_by(id) %>% mutate(surv1 = cumprod(1 - p.event1)) %>% ungroup()
#   
#   # Estimate risks from survival probabilities
#   # Risk = 1 - S(t)
#   treat0.surv$risk0 <- 1 - treat0.surv$surv0
#   treat1.surv$risk1 <- 1 - treat1.surv$surv1
#   
#   # Get the mean in each treatment group at each time point from 0 to 29 (30 time points in total)
#   treat0.surv$expo <- as.numeric(as.character(treat0.surv$expo))
#   risk0 <- aggregate(treat0.surv[c("expo", "time", "risk0")], by=list(treat0.surv$time), FUN=mean)[c("expo", "time", "risk0")] 
#   treat1.surv$expo <- as.numeric(as.character(treat1.surv$expo))
#   risk1 <- aggregate(treat1.surv[c("expo", "time", "risk1")], by=list(treat1.surv$time), FUN=mean)[c("expo", "time", "risk1")] 
#   
#   # Prepare data
#   graph.pred <- merge(risk0, risk1, by=c("time"))
#   # Edit data frame to reflect that risks are estimated at the END of each interval
#   graph.pred$time_0 <- graph.pred$time + 1
#   zero <- data.frame(cbind(0,0,0,1,0,0))
#   zero <- setNames(zero,names(graph.pred))
#   graph <- rbind(zero, graph.pred)
#   
#   ### Use pooled logistic regression estimates to compute causal estimates ###
#   
#   # 30-days risk in flu
#   risk0.plr <- graph$risk0[which(graph$time_0==K-1)]
#   risk0.plr
#   
#   # 30-days risk in COVID
#   risk1.plr <- graph$risk1[which(graph$time_0==K-1)]
#   risk1.plr
#   
#   # 30-days risk difference
#   rd.plr <- risk1.plr - risk0.plr
#   rd.plr
#   
#   # 30-days risk ratio
#   rr.plr <- risk1.plr / risk0.plr
#   rr.plr
#   
#   ggplot(graph, 
#          aes(x=time_0, y=risk)) + # set x and y axes
#     geom_line(aes(y = risk1, # create line for COVID group
#                   color = "OAC"),
#               size = 1.5) + 
#     geom_line(aes(y = risk0, # create line for no FLU group
#                   color = "NonUsers"),
#               size = 1.5) +
#     xlab("Years") + # label x axis
#     scale_x_continuous(limits = c(0, 60), # format x axis
#                        breaks=c(0, 12, 24, 36, 48, 60),
#                        labels =c(0:5) ) +
#     ylab("Cumulative Incidence (%)") + # label y axis
#     scale_y_continuous(limits=c(0, 0.6), # format y axis
#                        breaks=c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6),
#                        labels=c("0%","10%",  "20%", "30%", "40%","50%","60%")) +
#     theme_minimal()+ # set plot theme elements
#     scale_color_manual(name="OAC users",
#                        values=c("#386cb0","#f87f01"), # set colors
#                        breaks=c('NonUsers', 'OAC'),
#                        labels=c("None users","users")) +
#     theme(plot.title = element_text(face = "bold",
#                                     size = rel(1.2), hjust = 0.5),
#           text = element_text(),
#           panel.background = element_rect(colour = NA),
#           plot.background = element_rect(colour = NA),
#           # panel.border = element_rect(colour = NA),
#           axis.title = element_text(face = "bold",size = rel(1)),
#           axis.title.y = element_text(angle=90,vjust =2),
#           axis.title.x = element_text(vjust = -0.2),
#           axis.text = element_text(), 
#           axis.line = element_line(colour="black"),
#           axis.ticks = element_line(),
#           panel.grid.major = element_line(colour="#f0f0f0"),
#           panel.grid.minor = element_blank(),
#           legend.key = element_rect(colour = NA),
#           legend.position = "bottom",
#           legend.direction = "horizontal",
#           legend.key.size= unit(0.2, "cm"),
#           legend.title = element_text(face="italic"),
#           plot.margin=unit(c(10,5,5,5),"mm"),
#           strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
#           strip.text = element_text(face="bold")
#     ) -> plot_recurrent
# }

# incidence

{
  data4plot <- person_trial[dx.stroke_embo==0]
  pool.logistic.plot <-speedglm(out ~ expo + I(time) + I(timesqr) + Sex + age_indx + time_af_index + 
                                  chadsvas_score + dx.cbd + dx.copd + dx.dementia  + 
                                  dx.ra + dx.dm_com0 + dx.htn + antiplatelet + I(time*expo) + I(timesqr*expo),
                                family=binomial(link="logit"),data=data4plot, maxit = 100 )
  # Create dataset with all time points for each individual under each treatment level
  treat0 <- expandRows(data4plot[which(data4plot$time==0),], count=K, count.is.col=F) 
  treat0$time <- rep(seq(0, 59), nrow(data4plot[which(data4plot$time==0),])) 
  treat0$timesqr <- treat0$time^2
  
  # Under flu 
  treat0$expo <- as.numeric(as.character(treat0$expo))
  treat0$expo <- 0
  # Under covid 
  treat1 <- treat0
  treat1$expo <- 1
  
  
  # Extract predicted values from pooled logistic regression model for each person-time row
  # Predicted values correspond to discrete-time hazards
  treat0$p.event0 <- predict(pool.logistic.plot, treat0, type="response")
  treat1$p.event1 <- predict(pool.logistic.plot, treat1, type="response")
  
  # Obtain predicted survival probabilities from discrete-time hazards
  treat0.surv <- treat0 %>% group_by(id) %>% mutate(surv0 = cumprod(1 - p.event0)) %>% ungroup()
  treat1.surv <- treat1 %>% group_by(id) %>% mutate(surv1 = cumprod(1 - p.event1)) %>% ungroup()
  
  # Estimate risks from survival probabilities
  # Risk = 1 - S(t)
  treat0.surv$risk0 <- 1 - treat0.surv$surv0
  treat1.surv$risk1 <- 1 - treat1.surv$surv1
  
  # Get the mean in each treatment group at each time point from 0 to 29 (30 time points in total)
  treat0.surv$expo <- as.numeric(as.character(treat0.surv$expo))
  risk0 <- aggregate(treat0.surv[c("expo", "time", "risk0")], by=list(treat0.surv$time), FUN=mean)[c("expo", "time", "risk0")] 
  treat1.surv$expo <- as.numeric(as.character(treat1.surv$expo))
  risk1 <- aggregate(treat1.surv[c("expo", "time", "risk1")], by=list(treat1.surv$time), FUN=mean)[c("expo", "time", "risk1")] 
  
  # Prepare data
  graph.pred <- merge(risk0, risk1, by=c("time"))
  # Edit data frame to reflect that risks are estimated at the END of each interval
  graph.pred$time_0 <- graph.pred$time + 1
  zero <- data.frame(cbind(0,0,0,1,0,0))
  zero <- setNames(zero,names(graph.pred))
  graph <- rbind(zero, graph.pred)
  
  ### Use pooled logistic regression estimates to compute causal estimates ###
  
  # 30-days risk in flu
  risk0.plr <- graph$risk0[which(graph$time_0==K-1)]
  risk0.plr
  
  # 30-days risk in COVID
  risk1.plr <- graph$risk1[which(graph$time_0==K-1)]
  risk1.plr
  
  # 30-days risk difference
  rd.plr <- risk1.plr - risk0.plr
  rd.plr
  
  # 30-days risk ratio
  rr.plr <- risk1.plr / risk0.plr
  rr.plr
  
  ggplot(graph, 
         aes(x=time_0, y=risk)) + # set x and y axes
    geom_line(aes(y = risk1, # create line for COVID group
                  color = "OAC"),
              size = 1.5) + 
    geom_line(aes(y = risk0, # create line for no FLU group
                  color = "NonUsers"),
              size = 1.5) +
    xlab("Years") + # label x axis
    scale_x_continuous(limits = c(0, 60), # format x axis
                       breaks=c(0, 12, 24, 36, 48, 60),
                       labels =c(0:5) ) +
    ylab("Cumulative Incidence (%)") + # label y axis
    scale_y_continuous(limits=c(0, 0.6), # format y axis
                       breaks=c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6),
                       labels=c("0%","10%",  "20%", "30%", "40%","50%","60%")) +
    theme_minimal()+ # set plot theme elements
    scale_color_manual(name="OAC users",
                       values=c("#386cb0","#f87f01"), # set colors
                       breaks=c('NonUsers', 'OAC'),
                       labels=c("None users","users")) +
    theme(plot.title = element_text(face = "bold",
                                    size = rel(1.2), hjust = 0.5),
          text = element_text(),
          panel.background = element_rect(colour = NA),
          plot.background = element_rect(colour = NA),
          # panel.border = element_rect(colour = NA),
          axis.title = element_text(face = "bold",size = rel(1)),
          axis.title.y = element_text(angle=90,vjust =2),
          axis.title.x = element_text(vjust = -0.2),
          axis.text = element_text(), 
          axis.line = element_line(colour="black"),
          axis.ticks = element_line(),
          panel.grid.major = element_line(colour="#f0f0f0"),
          panel.grid.minor = element_blank(),
          legend.key = element_rect(colour = NA),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.key.size= unit(0.2, "cm"),
          legend.title = element_text(face="italic"),
          plot.margin=unit(c(10,5,5,5),"mm"),
          strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
          strip.text = element_text(face="bold")
    ) -> plot_inci
}


# incidence male

{
  data4plot <- person_trial[dx.stroke_embo==0 & Sex =="M"]
  pool.logistic.plot <-speedglm(out ~ expo + I(time) + I(timesqr)  + age_indx + time_af_index + 
                                  chadsvas_score + dx.cbd + dx.copd + dx.dementia  + 
                                  dx.ra + dx.dm_com0 + dx.htn + antiplatelet + I(time*expo) + I(timesqr*expo),
                                family=binomial(link="logit"),data=data4plot, maxit = 10000000 )
  # Create dataset with all time points for each individual under each treatment level
  treat0 <- expandRows(data4plot[which(data4plot$time==0),], count=K, count.is.col=F) 
  treat0$time <- rep(seq(0, 59), nrow(data4plot[which(data4plot$time==0),])) 
  treat0$timesqr <- treat0$time^2
  
  # Under flu 
  treat0$expo <- as.numeric(as.character(treat0$expo))
  treat0$expo <- 0
  # Under covid 
  treat1 <- treat0
  treat1$expo <- 1
  
  
  # Extract predicted values from pooled logistic regression model for each person-time row
  # Predicted values correspond to discrete-time hazards
  treat0$p.event0 <- predict(pool.logistic.plot, treat0, type="response")
  treat1$p.event1 <- predict(pool.logistic.plot, treat1, type="response")
  
  # Obtain predicted survival probabilities from discrete-time hazards
  treat0.surv <- treat0 %>% group_by(id) %>% mutate(surv0 = cumprod(1 - p.event0)) %>% ungroup()
  treat1.surv <- treat1 %>% group_by(id) %>% mutate(surv1 = cumprod(1 - p.event1)) %>% ungroup()
  
  # Estimate risks from survival probabilities
  # Risk = 1 - S(t)
  treat0.surv$risk0 <- 1 - treat0.surv$surv0
  treat1.surv$risk1 <- 1 - treat1.surv$surv1
  
  # Get the mean in each treatment group at each time point from 0 to 29 (30 time points in total)
  treat0.surv$expo <- as.numeric(as.character(treat0.surv$expo))
  risk0 <- aggregate(treat0.surv[c("expo", "time", "risk0")], by=list(treat0.surv$time), FUN=mean)[c("expo", "time", "risk0")] 
  treat1.surv$expo <- as.numeric(as.character(treat1.surv$expo))
  risk1 <- aggregate(treat1.surv[c("expo", "time", "risk1")], by=list(treat1.surv$time), FUN=mean)[c("expo", "time", "risk1")] 
  
  # Prepare data
  graph.pred <- merge(risk0, risk1, by=c("time"))
  # Edit data frame to reflect that risks are estimated at the END of each interval
  graph.pred$time_0 <- graph.pred$time + 1
  zero <- data.frame(cbind(0,0,0,1,0,0))
  zero <- setNames(zero,names(graph.pred))
  graph <- rbind(zero, graph.pred)
  
  ### Use pooled logistic regression estimates to compute causal estimates ###
  
  # 30-days risk in flu
  risk0.plr <- graph$risk0[which(graph$time_0==K-1)]
  risk0.plr
  
  # 30-days risk in COVID
  risk1.plr <- graph$risk1[which(graph$time_0==K-1)]
  risk1.plr
  
  # 30-days risk difference
  rd.plr <- risk1.plr - risk0.plr
  rd.plr
  
  # 30-days risk ratio
  rr.plr <- risk1.plr / risk0.plr
  rr.plr
  
  ggplot(graph, 
         aes(x=time_0, y=risk)) + # set x and y axes
    geom_line(aes(y = risk1, # create line for COVID group
                  color = "OAC"),
              size = 1.5) + 
    geom_line(aes(y = risk0, # create line for no FLU group
                  color = "NonUsers"),
              size = 1.5) +
    xlab("Years") + # label x axis
    scale_x_continuous(limits = c(0, 60), # format x axis
                       breaks=c(0, 12, 24, 36, 48, 60),
                       labels =c(0:5) ) +
    ylab("Cumulative Incidence (%)") + # label y axis
    scale_y_continuous(limits=c(0, 0.6), # format y axis
                       breaks=c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6),
                       labels=c("0%","10%",  "20%", "30%", "40%","50%","60%")) +
    theme_minimal()+ # set plot theme elements
    scale_color_manual(name="OAC users",
                       values=c("#386cb0","#f87f01"), # set colors
                       breaks=c('NonUsers', 'OAC'),
                       labels=c("None users","users")) +
    theme(plot.title = element_text(face = "bold",
                                    size = rel(1.2), hjust = 0.5),
          text = element_text(),
          panel.background = element_rect(colour = NA),
          plot.background = element_rect(colour = NA),
          # panel.border = element_rect(colour = NA),
          axis.title = element_text(face = "bold",size = rel(1)),
          axis.title.y = element_text(angle=90,vjust =2),
          axis.title.x = element_text(vjust = -0.2),
          axis.text = element_text(), 
          axis.line = element_line(colour="black"),
          axis.ticks = element_line(),
          panel.grid.major = element_line(colour="#f0f0f0"),
          panel.grid.minor = element_blank(),
          legend.key = element_rect(colour = NA),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.key.size= unit(0.2, "cm"),
          legend.title = element_text(face="italic"),
          plot.margin=unit(c(10,5,5,5),"mm"),
          strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
          strip.text = element_text(face="bold")
    ) -> plot_inci_male
}

# incidence female

{
  data4plot <- person_trial[dx.stroke_embo==0 & Sex =="F"]
  pool.logistic.plot <-speedglm(out ~ expo + I(time) + I(timesqr) + age_indx + time_af_index + 
                                  chadsvas_score + dx.cbd + dx.copd + dx.dementia  + 
                                  dx.ra + dx.dm_com0 + dx.htn + antiplatelet + I(time*expo) + I(timesqr*expo),
                                family=binomial(link="logit"),data=data4plot, maxit = 100 )
  # Create dataset with all time points for each individual under each treatment level
  treat0 <- expandRows(data4plot[which(data4plot$time==0),], count=K, count.is.col=F) 
  treat0$time <- rep(seq(0, 59), nrow(data4plot[which(data4plot$time==0),])) 
  treat0$timesqr <- treat0$time^2
  
  # Under flu 
  treat0$expo <- as.numeric(as.character(treat0$expo))
  treat0$expo <- 0
  # Under covid 
  treat1 <- treat0
  treat1$expo <- 1
  
  
  # Extract predicted values from pooled logistic regression model for each person-time row
  # Predicted values correspond to discrete-time hazards
  treat0$p.event0 <- predict(pool.logistic.plot, treat0, type="response")
  treat1$p.event1 <- predict(pool.logistic.plot, treat1, type="response")
  
  # Obtain predicted survival probabilities from discrete-time hazards
  treat0.surv <- treat0 %>% group_by(id) %>% mutate(surv0 = cumprod(1 - p.event0)) %>% ungroup()
  treat1.surv <- treat1 %>% group_by(id) %>% mutate(surv1 = cumprod(1 - p.event1)) %>% ungroup()
  
  # Estimate risks from survival probabilities
  # Risk = 1 - S(t)
  treat0.surv$risk0 <- 1 - treat0.surv$surv0
  treat1.surv$risk1 <- 1 - treat1.surv$surv1
  
  # Get the mean in each treatment group at each time point from 0 to 29 (30 time points in total)
  treat0.surv$expo <- as.numeric(as.character(treat0.surv$expo))
  risk0 <- aggregate(treat0.surv[c("expo", "time", "risk0")], by=list(treat0.surv$time), FUN=mean)[c("expo", "time", "risk0")] 
  treat1.surv$expo <- as.numeric(as.character(treat1.surv$expo))
  risk1 <- aggregate(treat1.surv[c("expo", "time", "risk1")], by=list(treat1.surv$time), FUN=mean)[c("expo", "time", "risk1")] 
  
  # Prepare data
  graph.pred <- merge(risk0, risk1, by=c("time"))
  # Edit data frame to reflect that risks are estimated at the END of each interval
  graph.pred$time_0 <- graph.pred$time + 1
  zero <- data.frame(cbind(0,0,0,1,0,0))
  zero <- setNames(zero,names(graph.pred))
  graph <- rbind(zero, graph.pred)
  
  ### Use pooled logistic regression estimates to compute causal estimates ###
  
  # 30-days risk in flu
  risk0.plr <- graph$risk0[which(graph$time_0==K-1)]
  risk0.plr
  
  # 30-days risk in COVID
  risk1.plr <- graph$risk1[which(graph$time_0==K-1)]
  risk1.plr
  
  # 30-days risk difference
  rd.plr <- risk1.plr - risk0.plr
  rd.plr
  
  # 30-days risk ratio
  rr.plr <- risk1.plr / risk0.plr
  rr.plr
  
  ggplot(graph, 
         aes(x=time_0, y=risk)) + # set x and y axes
    geom_line(aes(y = risk1, # create line for COVID group
                  color = "OAC"),
              size = 1.5) + 
    geom_line(aes(y = risk0, # create line for no FLU group
                  color = "NonUsers"),
              size = 1.5) +
    xlab("Years") + # label x axis
    scale_x_continuous(limits = c(0, 60), # format x axis
                       breaks=c(0, 12, 24, 36, 48, 60),
                       labels =c(0:5) ) +
    ylab("Cumulative Incidence (%)") + # label y axis
    scale_y_continuous(limits=c(0, 0.6), # format y axis
                       breaks=c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6),
                       labels=c("0%","10%",  "20%", "30%", "40%","50%","60%")) +
    theme_minimal()+ # set plot theme elements
    scale_color_manual(name="OAC users",
                       values=c("#386cb0","#f87f01"), # set colors
                       breaks=c('NonUsers', 'OAC'),
                       labels=c("None users","users")) +
    theme(plot.title = element_text(face = "bold",
                                    size = rel(1.2), hjust = 0.5),
          text = element_text(),
          panel.background = element_rect(colour = NA),
          plot.background = element_rect(colour = NA),
          # panel.border = element_rect(colour = NA),
          axis.title = element_text(face = "bold",size = rel(1)),
          axis.title.y = element_text(angle=90,vjust =2),
          axis.title.x = element_text(vjust = -0.2),
          axis.text = element_text(), 
          axis.line = element_line(colour="black"),
          axis.ticks = element_line(),
          panel.grid.major = element_line(colour="#f0f0f0"),
          panel.grid.minor = element_blank(),
          legend.key = element_rect(colour = NA),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.key.size= unit(0.2, "cm"),
          legend.title = element_text(face="italic"),
          plot.margin=unit(c(10,5,5,5),"mm"),
          strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
          strip.text = element_text(face="bold")
    ) -> plot_inci_female
}


library(cowplot)

legend <- get_legend(
  plot_main + theme(legend.box.margin = margin(0, 0, 0, 12))
)
cowplot::get_plot_component(plot_main, 'guide-box-top', return_all = TRUE)


all_plot <- 
plot_grid(plot_grid(plot_main +labs(subtitle = "A: Primary or secondary prevention")+ theme(legend.position="none"), 
                    plot_inci +labs(subtitle = "B: Primary prevention")+ theme(legend.position="none"), 
                    plot_male +labs(subtitle = "C: Subgroup analysis-male")+ theme(legend.position="none"), 
                    plot_female+labs(subtitle = "D: Subgroup analysis-female")+ theme(legend.position="none"), 
                    plot_inci_male+labs(subtitle = "E: Subgroup analysis-primary prevention in the male")+ theme(legend.position="none"),
                    plot_inci_female+labs(subtitle = "F: Subgroup analysis-primary prevention in the female")+ theme(legend.position="none"),
                    ncol=2),
          cowplot::get_plot_component(plot_main, 'guide-box-bottom', return_all = TRUE),nrow=2,rel_heights = c(1,.1))

save_plot <-function(filename, plot, width_in, height_in, dpi =300) {
  if(tools::file_ext(filename) =="pdf") {
    ggsave(filename, plot = plot, width = width_in, height = height_in, units ="in")
  }else{
    ggsave(filename, plot = plot, width = width_in * dpi, height = height_in * dpi, units ="px", dpi = dpi)
  }
}

# 使用示例
save_plot("haemorraghic stroke cummulative incidence.pdf", all_plot, width_in = 10, height_in =12)
# save_plot("plot.png", p, width_in =8, height_in =6)
# save_plot("plot.jpeg", p, width_in =8, height_in =6)


