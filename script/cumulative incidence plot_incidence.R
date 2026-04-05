generate_risk_plot <- function(formula, data, year, 
                               colors = c("NonUsers" = "#386cb0", "OAC" = "#f87f01"),
                               labels = c("NonUsers" = "None users", "OAC" = "users")) {
  
  # 1. Model Fitting --------------------------------------------------------
  model <- speedglm(formula = formula,
                    family = binomial(link = "logit"),
                    data = data,
                    maxit = 1000000000)
  
  # 2. Counterfactual Prediction --------------------------------------------
  # Expand baseline data for all time points
  treat0 <- expandRows(data[data$time == 0, ], count = year*12 , count.is.col = FALSE)
  treat0$time <- rep(seq(0, (year*12-1)), nrow(data[data$time == 0, ]))
  treat0$timesqr <- treat0$time^2
  
  # Create treatment scenarios
  treat0$expo <- 0  # Control group
  treat1 <- copy(treat0)
  treat1$expo <- 1  # Treatment group
  
  # 3. Survival Probability Calculation -------------------------------------
  treat0$p.event0 <- predict(model, treat0, type = "response")
  treat1$p.event1 <- predict(model, treat1, type = "response")
  
  calculate_survival <- function(dt, p_col) {
    dt %>%
      group_by(id) %>%
      mutate(surv = cumprod(1 - .data[[p_col]])) %>%
      ungroup() %>%
      mutate(risk = 1 - surv)
  }
  
  # 4. Risk Aggregation ----------------------------------------------------
  treat0.surv <- calculate_survival(treat0, "p.event0")
  treat1.surv <- calculate_survival(treat1, "p.event1")
 
  aggregate_risk <- function(dt) {
    aggregate(dt[c("expo", "time", "risk")], 
              by = list(time = dt$time), 
              FUN = mean)[c("expo", "time", "risk")]
  }
  
  risk0 <- aggregate_risk(treat0.surv)
  risk1 <- aggregate_risk(treat1.surv)
  
  # 5. Data Merging --------------------------------------------------------
  graph_data <- merge(risk0, risk1, by = "time") %>%
    mutate(time_0 = time + 1)  # Shift time for interval end
  
  # Add baseline zero point
  zero_row <- data.frame(time = 0, expo.x=0, risk.x = 0, expo.y=0, risk.y = 0, time_0 = 0)
  graph_data <- rbind(zero_row, graph_data)
  
  # 6. Visualization -------------------------------------------------------
  plot_obj <- ggplot(graph_data, aes(x = time_0)) +
    geom_line(aes(y = risk.y, color = "OAC"), size = 1.5) +
    geom_line(aes(y = risk.x, color = "NonUsers"), size = 1.5) +
    scale_x_continuous(
      # name = "Years",
      name = "",
      limits = c(0, year*12),
      breaks = seq(0, year)*12,
      labels = 0:year,
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      # name = "Cumulative Incidence (%)",
      name = "",
      limits = c(0, 0.8),
      labels = scales::percent_format(accuracy = 1),
      expand = c(0, 0)
    ) +
    scale_color_manual(
      name = "OAC users",
      values = colors,
      labels = labels
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold",
                                    size = rel(1.2), hjust = 0.5),
          text = element_text(),
          panel.background = element_rect(colour = NA),
          plot.background = element_rect(colour = NA),
          # panel.border = element_rect(colour = NA),
          axis.title = element_text(face = "bold",size = rel(1)),
          axis.title.y = element_text(angle=90,vjust =2,size=16),
          axis.title.x = element_text(vjust = -0.2,size=16),
          axis.text = element_text(), 
          axis.line = element_line(colour="black"),
          axis.ticks = element_line(),
          panel.grid.major = element_line(colour="#f0f0f0"),
          panel.grid.minor = element_blank(),
          legend.key = element_rect(colour = NA),
          legend.text = element_text(size=12),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.key.size= unit(0.2, "cm"),
          legend.title = element_text(face="italic",size=14),
          plot.margin=unit(c(10,5,5,5),"mm"),
          strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
          strip.text = element_text(face="bold")
    )
  
  # 7. Risk Metrics --------------------------------------------------------
  final_risk <- tail(graph_data, 1)
  risk_metrics <- list(
    risk_ratio = final_risk$risk.y / final_risk$risk.x,
    risk_difference = final_risk$risk.y - final_risk$risk.x
  )
  
  return(list(plot = plot_obj, metrics = risk_metrics))
}

# Example Usage -----------------------------------------------------------
year <- 3
# Generate main plot
main_plot <- generate_risk_plot(
  formula = formula_rr_main_cum_plot,
  data = person_trial,
  year = year
)

# Generate sex-stratified plots
male_plot <- generate_risk_plot(
  formula = formula_rr_sex_cum_plot,
  data = person_trial[Sex == "M"],
  year = year
)

female_plot <- generate_risk_plot(
  formula = formula_rr_sex_cum_plot,
  data = person_trial[Sex == "F"], 
  year = year
)

all_stroke <- list(main_plot, male_plot, female_plot)
all_bleeding <- list(main_plot, male_plot, female_plot)
ischemic_plot <- list(main_plot, male_plot, female_plot)
haemorraghic_plot <- list(main_plot, male_plot, female_plot)
mortality <- list(main_plot, male_plot, female_plot)

# Combine plots using patchwork
library(patchwork)
(ischemic_plot[[1]]$plot + theme() | haemorraghic_plot[[1]]$plot) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold"))

all_stroke[[1]]$plot + ggtitle('Composite stroke and mortality') +
all_bleeding[[1]]$plot + ggtitle('Composite bleeding events')+
ischemic_plot[[1]]$plot + ggtitle('Ischemic stroke')+
haemorraghic_plot[[1]]$plot + ggtitle('Hemorraghic stroke') + 
# mortality[[1]]$plot + ggtitle("All-cause mortality")+
wrap_elements(panel = textGrob('All population', rot=-90))+
all_stroke[[2]]$plot + 
all_bleeding[[2]]$plot +
ischemic_plot[[2]]$plot + 
haemorraghic_plot[[2]]$plot +
# mortality[[2]]$plot+
wrap_elements(panel = textGrob('Male', rot=-90))+
all_stroke[[3]]$plot +
all_bleeding[[3]]$plot  +
ischemic_plot[[3]]$plot + 
haemorraghic_plot[[3]]$plot +
# mortality[[3]]$plot +
wrap_elements(panel = textGrob('Female', rot=-90))+
plot_layout(nrow =  2,ncol=2,widths = c(1,1,1,1,.1), guides = "collect",axis_titles = "collect") & theme(legend.position = 'bottom') 



all_stroke[[1]]$plot + ggtitle('Composite stroke and mortality') +  
  all_bleeding[[1]]$plot + ggtitle('Composite bleeding events')+
  wrap_elements(panel = textGrob('All population', rot=-90))+
  all_stroke[[2]]$plot +
  all_bleeding[[2]]$plot +
  wrap_elements(panel = textGrob('Male', rot=-90))+
  all_stroke[[3]]$plot +
  all_bleeding[[3]]$plot  +
  wrap_elements(panel = textGrob('Female', rot=-90))+
  ischemic_plot[[1]]$plot + ggtitle('Ischemic stroke')+
  haemorraghic_plot[[1]]$plot + ggtitle('Hemorraghic stroke') + 
  wrap_elements(panel = textGrob('All population', rot=-90))+
  ischemic_plot[[2]]$plot + 
  haemorraghic_plot[[2]]$plot +
  wrap_elements(panel = textGrob('Male', rot=-90))+
  ischemic_plot[[3]]$plot + 
  haemorraghic_plot[[3]]$plot +
  wrap_elements(panel = textGrob('Female', rot=-90))+
  plot_layout(nrow =  8, ncol=3, guides = "collect",axis_titles = "collect",widths = c(1,1,.1)) & theme(legend.position = 'bottom') 





# plot --------------------------------------------------------------------

{
  pool.logistic.plot <-speedglm(formula_rr_main_cum_plot,
                              family=binomial(link="logit"),data=person_trial)

  # Create dataset with all time points for each individual under each treatment level
  treat0 <- expandRows(person_trial[which(person_trial$time==0),], count=K, count.is.col=F) 
  treat0$time <- rep(seq(0, K-1), nrow(person_trial[which(person_trial$time==0),]))
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
    scale_x_continuous(limits = c(0, K), # format x axis
                       breaks=c(0:year*12),
                       labels =c(0:year) ) +
    ylab("Cumulative Incidence (%)") + # label y axis
    scale_y_continuous(limits=c(0, 0.4), # format y axis
                       breaks=c(0, 0.1, 0.2, 0.3, 0.4),
                       labels=c("0%","10%",  "20%", "30%", "40%")) +
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
  pool.logistic.plot <-speedglm(formula_rr_sex_cum_plot,
                                family=binomial(link="logit"),data=data4plot, maxit = 100000 )
  # Create dataset with all time points for each individual under each treatment level
  treat0 <- expandRows(data4plot[which(data4plot$time==0),], count=K, count.is.col=F) 
  treat0$time <- rep(seq(0, K-1), nrow(data4plot[which(data4plot$time==0),])) 
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
    scale_x_continuous(limits = c(0, K), # format x axis
                       breaks=c(0:year*12),
                       labels =c(0:year) ) +
    ylab("Cumulative Incidence (%)") + # label y axis
    scale_y_continuous(limits=c(0, 0.4), # format y axis
                       breaks=c(0, 0.1, 0.2, 0.3, 0.4),
                       labels=c("0%","10%",  "20%", "30%", "40%")) +
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

# female
{
  data4plot <- person_trial[Sex=="F"]
  pool.logistic.plot <-speedglm(formula_rr_sex_cum_plot,
                                family=binomial(link="logit"),data=data4plot, maxit = 100 )
  # Create dataset with all time points for each individual under each treatment level
  treat0 <- expandRows(data4plot[which(data4plot$time==0),], count=K, count.is.col=F) 
  treat0$time <- rep(seq(0, K-1), nrow(data4plot[which(data4plot$time==0),])) 
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
    scale_x_continuous(limits = c(0, K), # format x axis
                       breaks=c(0:year*12),
                       labels =c(0:year) ) +
    ylab("Cumulative Incidence (%)") + # label y axis
    scale_y_continuous(limits=c(0, 0.4), # format y axis
                       breaks=c(0, 0.1, 0.2, 0.3, 0.4),
                       labels=c("0%","10%",  "20%", "30%", "40%")) +
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

library(cowplot)

legend <- get_legend(
  plot_main + theme(legend.box.margin = margin(0, 0, 0, 12))
)
cowplot::get_plot_component(plot_main, 'guide-box-top', return_all = TRUE)


all_plot <- 
  plot_grid(plot_grid(plot_main +labs(subtitle = "A: Primary or secondary prevention")+ theme(legend.position="none"), 
                      plot_male +labs(subtitle = "C: Subgroup analysis-male")+ theme(legend.position="none"), 
                      plot_female+labs(subtitle = "D: Subgroup analysis-female")+ theme(legend.position="none"), 
                      ncol=1),
            cowplot::get_plot_component(plot_main, 'guide-box-bottom', return_all = TRUE),nrow=2,rel_heights = c(1,.1))

all_plot


# 使用示例
save_plot("haemorraghic stroke cummulative incidence.pdf", all_plot, width_in = 10, height_in =12)
# save_plot("plot.png", p, width_in =8, height_in =6)
# save_plot("plot.jpeg", p, width_in =8, height_in =6)


