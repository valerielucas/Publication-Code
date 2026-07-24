## Valerie Lucas
## 2026-07-23

## dataset_construction: function to take in participant list, baseline data, and depression data and return analytic dataset
## uid_list: vector of participant IDs to be included in the analysis
## adapt_visit201: dataframe of baseline data
## adapt_depression: dataframe of follow-up depression data for all visits
dataset_construction <- function(uid_list, adapt_visit201, adapt_depression){
  
  # create shell variables
  visit_list <- c(201, 301, 401, 501, 701, 801)
  visit_number_list <- c(0:3,5,7)
  
  uid_all <- rep(uid_list, times = length(visit_list))
  visit_all <- rep(visit_list, each = length(uid_list))
  visit_number_all <- rep(visit_number_list, each = length(uid_list))
  
  # create empty data frame
  adapt_empty <- data.frame(uid = uid_all, visit = visit_all, visit_number = visit_number_all)
  
  # fill in exposure and covariate data
  adapt_baseline <- merge(adapt_empty, adapt_visit201, all.x = TRUE, by = c("uid", "visit")) |>
    group_by(uid) |>
    arrange(visit) |>
    fill(arm_name, double_orphan, not_care, food_insecure, ipv_physical, sex_violence, school_violence, someprimaryorless_mother, consumption, consumption_rs1, consumption_rs2, consumption_rs3, yw_age, cdi_depressed7, .direction = "down")
  
  # merge in depression data -- drop all depression data for which there is no one in the baseline data
  adapt <- merge(adapt_baseline, adapt_depression, all.x = TRUE, by = c("uid", "visit"))
  
  # create new variable indicator for when
  adapt$depression_nonmissing <- !is.na(adapt$depression)

  
  return(adapt)
  
}


# function to obtain point estimates for depression via g-computation

## adapt_no201: analytic dataset with 201 data removed
## aces: vector of ACE names
## gcomp_model_list: list of g-computation models for each ACE

## returns: dataframe of point estimatess
point_estimates_gcomp <- function(adapt_no201, aces, gcomp_model_list){
  
  # create empty data frame for results
  adapt_results <- data.frame(ace = character(0), visit_number = numeric(0), prev_diff = numeric(0))
  
  # single effect for each ACE
  
  for(i in 1:length(aces)){
    
    # create copy of dataset
    adapt_no201_copy <- adapt_no201
    
    # create weighted g-comp model
    gcomp_model_boot <- glm(gcomp_model_list[[i]], family = binomial(link = "logit"), data = adapt_no201_copy)
    
    # create g-comp model
    gcomp_model <- glm(gcomp_model_list[[i]], family = binomial(link = "logit"), data = adapt_no201_copy)
    
    # set plan for ACE to 1
    adapt_no201_copy[, (names(adapt_no201_copy) == aces[i])] <- 1
    
    # make predictions for ACE
    gcomp_model1 <- data.frame(probability1 = predict(gcomp_model, type = "response", newdata = adapt_no201_copy))
    
    # set plan for ACE to 0
    adapt_no201_copy[, (names(adapt_no201_copy) == aces[i])] <- 0
    gcomp_model0 <- data.frame(probability0 = predict(gcomp_model, type = "response", newdata = adapt_no201_copy))
    
    # bind g-comp predictions with the dataset
    adapt_result_loading <- cbind(adapt_no201_copy, gcomp_model1, gcomp_model0)
    
    # reframe to calculate prevalence difference
    adapt_result_loading <- adapt_result_loading |>
      group_by(visit_number) |>
      reframe(ace = aces[i],
              prev1 = mean(probability1),
              prev0 = mean(probability0),
              prev_diff = prev1 - prev0,
              prev_ratio = prev1 / prev0
              )
    
    # bind the results for this ACE to the overall results
    adapt_results <- rbind(adapt_results, adapt_result_loading)
  }
  
  # joint effect of physical IPV and sexual violence
  
  # reset dataset
  adapt_no201_copy <- adapt_no201
  
  # create g-comp model
  gcomp_model <- glm(gcomp_model_list[[7]], family = binomial(link = "logit"), data = adapt_no201_copy)
  
  # set plan for ACE to 1
  adapt_no201_copy[, (names(adapt_no201_copy) == aces[4] | names(adapt_no201_copy) == aces[5])] <- 1
  
  # make predictions for ACE
  gcomp_model1 <- data.frame(probability1 = predict(gcomp_model, type = "response", newdata = adapt_no201_copy))
  
  # set plan for ACE to 0
  adapt_no201_copy[, (names(adapt_no201_copy) == aces[4] | names(adapt_no201_copy) == aces[5])] <- 0
  gcomp_model0 <- data.frame(probability0 = predict(gcomp_model, type = "response", newdata = adapt_no201_copy))
  
  # bind g-comp predictions with the dataset
  adapt_result_loading <- cbind(adapt_no201_copy, gcomp_model1, gcomp_model0)
  
  # reframe to calculate prevalence difference
  adapt_result_loading <- adapt_result_loading |>
    group_by(visit_number) |>
    reframe(ace = "ipv_sexual_violence",
            prev1 = mean(probability1),
            prev0 = mean(probability0),
            prev_diff = prev1 - prev0,
            prev_ratio = prev1 / prev0)
  
  # bind the results for this ACE to the overall results
  adapt_results <- rbind(adapt_results, adapt_result_loading)
  
  # return overall results
  return(adapt_results)
}

# point_estimates_ipw: function to obtain point estimates for depression via inverse probability weighting

## adapt_no201: analytic dataset with 201 data removed
## aces: vector of ACE names
## iptw_model_list: list of IPTW models for each ACE
## iptw_model_list_st: list of IPTW stabilizing numerator models for each ACE
## ipcw_model_list: list of IPCW models for each ACE
## ipcw_st: IPCW stabilizing numerator model

## returns: dataframe of point estimates 
point_estimates_ipw <- function(adapt_no201, aces, iptw_model_list, iptw_model_list_st, ipcw_model_list, ipcw_st){
  
  # create frame for results
  adapt_results <- data.frame(ace = character(0), visit_number = numeric(0), prev1 = numeric(0), prev0 = numeric(0), prev_diff = numeric(0), prev_ratio = numeric(0))
  
  for(i in 1:length(aces)){
    
    # reset dataset
    adapt_no201_copy <- adapt_no201
    
    # save ace to new var
    adapt_no201_copy$ace_var <- adapt_no201_copy[, (names(adapt_no201_copy) == aces[i])]
    
    # fit iptw model
    iptw_model <- glm(iptw_model_list[[i]], family = binomial(link = "logit"), data = adapt_no201_copy)
    
    # predict fromm iptw model
    adapt_no201_copy$pi_a <- predict(iptw_model, type = "response", newdata = adapt_no201_copy)
    
    # fit iptw stable model
    iptw_model_st <- glm(iptw_model_list_st[[i]], family = binomial(link = "logit"), data = adapt_no201_copy)
    
    # predict fromm iptw model
    adapt_no201_copy$pi_a_st <- predict(iptw_model_st, type = "response", newdata = adapt_no201_copy)
    
    # calculate IPTW from probability of being assigned treatment
    adapt_no201_copy$iptw_a <- adapt_no201_copy$ace_var * adapt_no201_copy$pi_a_st / adapt_no201_copy$pi_a + (1 - adapt_no201_copy$ace_var) * (1 - adapt_no201_copy$pi_a_st) / (1 - adapt_no201_copy$pi_a)
    
    # fit IPCW model
    ipcw_model <- glm(ipcw_model_list[[i]], family = binomial(link = "logit"), data = adapt_no201_copy)
    
    # predict from IPCW model
    adapt_no201_copy$pi_c <- predict(ipcw_model, type = "response", newdata = adapt_no201_copy)
    
    # fit IPCW stable model
    ipcw_model_st <- glm(ipcw_st, family = binomial(link = "logit"), data = adapt_no201_copy)
    
    # predict from IPCW model
    adapt_no201_copy$pi_c_st <- predict(ipcw_model_st, type = "response", newdata = adapt_no201_copy)
    
    # calculate IPCW from probability of being observed
    adapt_no201_copy$ipcw <- adapt_no201_copy$pi_c_st / adapt_no201_copy$pi_c
    
    # calculate overall weight
    adapt_no201_copy$ipw <- adapt_no201_copy$ipcw * adapt_no201_copy$iptw
    
    adapt_result_loading <- adapt_no201_copy |>
      group_by(visit_number) |>
      reframe(ace = aces[i],
              prev1 = sum(depression * ipw * ace_var, na.rm = TRUE) / sum(ipw * ace_var * !is.na(depression)),
              prev0 = sum(depression * ipw * (1 - ace_var), na.rm = TRUE) / sum(ipw * (1 - ace_var) * !is.na(depression)),
              prev_diff = prev1 - prev0,
              prev_ratio = prev1 / prev0)
    
    
    # bind the results for this ACE to the overall results
    adapt_results <- rbind(adapt_results, adapt_result_loading)
    
  }
  
  # joint effect of physical IPV
  
  # reset dataset
  adapt_no201_copy <- adapt_no201
  
  # fit iptw A model
  iptw_model_a <- glm(iptw_model_list[[4]], family = binomial(link = "logit"), data = adapt_no201_copy)
  
  # predict from iptw A model
  adapt_no201_copy$pi_a <- predict(iptw_model_a, type = "response", newdata = adapt_no201_copy)
  
  # fit iptw stable model
  iptw_model_st <- glm(iptw_model_list_st[[4]], family = binomial(link = "logit"), data = adapt_no201_copy)
  
  # predict fromm iptw model
  adapt_no201_copy$pi_a_st <- predict(iptw_model_st, type = "response", newdata = adapt_no201_copy)
  
  # calculate IPTW from probability of being assigned treatment
  adapt_no201_copy$iptw_a <- adapt_no201_copy$ipv_physical * adapt_no201_copy$pi_a_st / adapt_no201_copy$pi_a + (1 - adapt_no201_copy$ipv_physical) * (1 - adapt_no201_copy$pi_a_st) / (1 - adapt_no201_copy$pi_a)
  
  # fit iptw B model
  iptw_model_b <- glm(iptw_model_list[[7]], family = binomial(link = "logit"), data = adapt_no201_copy)
  
  # predict from iptw B model
  adapt_no201_copy$pi_b <- predict(iptw_model_b, type = "response", newdata = adapt_no201_copy)
  
  # fit iptw stable model
  iptw_model_st <- glm(iptw_model_list_st[[5]], family = binomial(link = "logit"), data = adapt_no201_copy)
  
  # predict fromm iptw model
  adapt_no201_copy$pi_b_st <- predict(iptw_model_st, type = "response", newdata = adapt_no201_copy)
  
  # calculate IPTW from probability of being assigned treatment
  adapt_no201_copy$iptw_b <- adapt_no201_copy$sex_violence * adapt_no201_copy$pi_b_st / adapt_no201_copy$pi_b + (1 - adapt_no201_copy$sex_violence) * (1 - adapt_no201_copy$pi_b_st) / (1 - adapt_no201_copy$pi_b)
  
  # multiply weights
  adapt_no201_copy$iptw <-  adapt_no201_copy$iptw_a *  adapt_no201_copy$iptw_b
  
  # fit IPCW model
  ipcw_model <- glm(ipcw_model_list[[7]], family = binomial(link = "logit"), data = adapt_no201_copy)
  
  # predict from IPCW model
  adapt_no201_copy$pi_c <- predict(ipcw_model, type = "response", newdata = adapt_no201_copy)
  
  # fit IPCW stable model
  ipcw_model_st <- glm(ipcw_st, family = binomial(link = "logit"), data = adapt_no201_copy)
  
  # predict from IPCW model
  adapt_no201_copy$pi_c_st <- predict(ipcw_model_st, type = "response", newdata = adapt_no201_copy)
  
  # calculate IPCW from probability of being observed
  adapt_no201_copy$ipcw <- adapt_no201_copy$pi_c_st / adapt_no201_copy$pi_c
  
  # calculate overall weight
  adapt_no201_copy$ipw <- adapt_no201_copy$ipcw * adapt_no201_copy$iptw
  
  # reframe to calculate prevalence difference
  adapt_result_loading <- adapt_no201_copy |>
    group_by(visit_number) |>
    reframe(ace = "ipv_sexual_violence",
            prev1 = sum(depression * ipw * ipv_physical * sex_violence, na.rm = TRUE) / sum(ipw * ipv_physical * sex_violence * !is.na(depression)),
            prev0 = sum(depression * ipw * (1 - ipv_physical) * (1 - sex_violence), na.rm = TRUE) / sum(ipw * (1 - ipv_physical) * (1 - sex_violence) * !is.na(depression)),
            prev_diff = prev1 - prev0,
            prev_ratio = prev1 / prev0)
  
  # bind the results for this ACE to the overall results
  adapt_results <- rbind(adapt_results, adapt_result_loading)
  
  
  return(adapt_results)
}




## function to do a cluster bootsstrap of our analysis

## reps: number of reps in the bootstrap
## adapt_visit201: 
## adapt_depression:
## aces: vector of ACEs in the analysis
## analysis: function with analytic strategy
## point_estimates: point estimates from point_estimates_x function
## gcomp_model_list: if using point_estimates_gcomp, list of g-computation regression equations
## iptw_model_list: if using point_estimates_ipw, list of IPTW models for each ACE
## iptw_model_list_st: if using point_estimates_ipw, list of IPTW stabilizing numerator models for each ACE
## ipcw_model_list: if using point_estimates_ipw, list of IPCW models for each ACE
## ipcw_st: if using point_estimates_ipw, IPCW stabilizing numerator model

## returns: dataframe of point estimates and 95% confidence intervals
bootstrap_analysis <- function(reps, adapt_visit201, adapt_depression, aces, analysis, point_estimates, gcomp_model_list = NULL, iptw_model_list = NULL, iptw_model_list_st = NULL, ipcw_model_list = NULL, ipcw_st = NULL){
  
  # set up empty results
  adapt_boot_results <- data.frame(ace = character(0), visit_number = numeric(0), prev_diff = numeric(0))
  
  # set seed
  for(i in 1:reps){
    
    # sample UIDs with replacement
    uid_list <- sample(x = adapt_visit201$uid, size = length(adapt_visit201$uid), replace = TRUE)
    
    # create analytic dataset
    adapt_boot <- dataset_construction(uid_list, adapt_visit201, adapt_depression)
    
    # remove baseline visit (201), no outcomes reported
    adapt_boot <- adapt_boot |>
      filter(visit != 201)
    
    # conduct analysis
    if(!is.null(gcomp_model_list)){
      adapt_boot_result_loading <- analysis(adapt_boot, aces, gcomp_model_list = gcomp_model_list)}
    if(!is.null(iptw_model_list)){
      adapt_boot_result_loading <- analysis(adapt_boot, aces, iptw_model_list = iptw_model_list, iptw_model_list_st = iptw_model_list_st, ipcw_model_list = ipcw_model_list, ipcw_st = ipcw_st)
    }

    # bind the results for this ACE to the overall results
    adapt_boot_results <- rbind(adapt_boot_results, adapt_boot_result_loading)
    
  }
  
  # take the standard error of all replicates
  adapt_boot_results_se <- adapt_boot_results |>
    group_by(visit_number, ace) |>
    reframe(se_prev1 = sd(prev1),
            se_prev0 = sd(prev0),
            se_prev_diff = sd(prev_diff),
            se_log_prev_ratio = sd(log(prev1) - log(prev0)))
  
  # merge dataframe of standard errors from bootstrap with dataframe of point estimates
  point_estimates_se <- merge(point_estimates, adapt_boot_results_se, by = c("visit_number", "ace"))
  
  # calculate 95% confidence intervals
  adapt_results_ci <- point_estimates_se |>
    mutate(prev1_lower = prev1 - 1.96*se_prev1,
           prev1_upper = prev1 + 1.96*se_prev1,
           prev0_lower = prev0 - 1.96*se_prev0, 
           prev0_upper = prev0 + 1.96*se_prev0,
           
           prev1_lower = ifelse(prev1_lower < 0, 0, prev1_lower),
           prev0_lower = ifelse(prev0_lower < 0, 0, prev0_lower),
           
           prev_diff_lower = prev_diff - 1.96*se_prev_diff,
           prev_diff_upper = prev_diff + 1.96*se_prev_diff,
           
           prev_ratio_lower = exp(log(prev_ratio) - 1.96*se_log_prev_ratio),
           prev_ratio_upper = exp(log(prev_ratio) + 1.96*se_log_prev_ratio)
           )
  
  

  # return results dataframe
  return(adapt_results_ci)
  
}


## results_report: function to generate chart of estimates and

## adapt_results_ci: dataframe of point estimates and 95% confidence intervals
## file_path: character indicating file path
## label: character indicating label attached to file

## returns: NULL
results_report <- function(adapt_results_ci, file_path, label){
  
  # use model to predict missing data
  for(i in 1:length(aces)){
    
    # create plot for each individual ACE
    plot_loading <- adapt_results_ci |>
      filter(ace == aces[i]) |>
      ggplot() + 
      geom_line(aes(x = visit_number, y = prev0, color = "No ACE")) +
      geom_point(aes(x = visit_number, y = prev0, color = "No ACE")) +
      geom_ribbon(aes(x = visit_number, ymin = prev0_lower, ymax = prev0_upper), fill = "#332288", alpha = 0.25) +
      geom_line(aes(x = visit_number, y = prev1, color = "ACE")) +
      geom_point(aes(x = visit_number, y = prev1, color = "ACE")) +
      geom_ribbon(aes(x = visit_number, ymin = prev1_lower, ymax = prev1_upper), fill = "#CC6677", alpha = 0.25) +
      xlab(NULL) +
      ggtitle(paste0("Depression prevalence by experience of ", ace_names[i])) +
      ylab("Depression prevalence") +
      scale_x_continuous(breaks = c(1:3,5,7), labels = c("Year 1", "Year 2", "Year 3", "Year 5", "Year 7")) +
      ylim(0, 0.7) +
      scale_color_manual(name = "Exposure", labels = c("ACE", "No ACE"), values = c("#CC6677", "#332288")) +
      theme_classic()
    
    
    variable_name <- paste0("plot", i)
    assign(variable_name, plot_loading)
    
  }
    
    # create plot for joint effect
    plot7 <- plot_loading <- adapt_results_ci |>
      filter(ace == "ipv_sexual_violence") |>
      ggplot() + 
      geom_line(aes(x = visit_number, y = prev0, color = "Neither ACE")) +
      geom_point(aes(x = visit_number, y = prev0, color = "Neither ACE")) +
      geom_ribbon(aes(x = visit_number, ymin = prev0_lower, ymax = prev0_upper), fill = "#332288", alpha = 0.25) +
      geom_line(aes(x = visit_number, y = prev1, color = "Both ACEs")) +
      geom_point(aes(x = visit_number, y = prev1, color = "Both ACEs")) +
      geom_ribbon(aes(x = visit_number, ymin = prev1_lower, ymax = prev1_upper), fill = "#CC6677", alpha = 0.25) +
      xlab(NULL) +
      ggtitle(paste0("Depression prevalence by experience of both physical IPV and sexual violence")) +
      ylab("Depression prevalence") +
      scale_x_continuous(breaks = c(1:3,5,7), labels = c("Year 1", "Year 2", "Year 3", "Year 5", "Year 7")) +
      ylim(0, 0.7) +
      scale_color_manual(name = "Exposure", labels = c("Both ACEs", "Neither ACE"), values = c("#CC6677", "#332288")) +
      theme_classic()
    
  # Arrange plots with `arrangeGrob`
  adapt_plots <- arrangeGrob(plot1, plot2, plot3, plot4, plot5, plot7, plot6, nrow = 4)
  
  # Save the Grob to a PDF file
  ggsave(adapt_plots, file = paste0(file_path, "Charts/adapt_plots_", label, ".pdf"), width = 14, height = 16, units = "in")
  
  # round prevalence difference
  adapt_results_ci$pd_percent <- round(100*adapt_results_ci$prev_diff, digits = 1)
  
  # round prevalence difference lower confidence limit
  adapt_results_ci$pd_lcl_percent <- round(100*adapt_results_ci$prev_diff_lower, digits = 1)
  
  # ound  prevalence difference upper confidence limit 
  adapt_results_ci$pd_ucl_percent <- round(100*adapt_results_ci$prev_diff_upper, digits = 1)
  
  # make nice format for table
  adapt_results_ci$pd_report <- paste0(adapt_results_ci$pd_percent, " (", adapt_results_ci$pd_lcl_percent, ", ", adapt_results_ci$pd_ucl_percent, ")")
  
  # pivot wider just the `pd_report` column and arrange in order of table in manuscript
  adapt_report_pd <- adapt_results_ci |>
    select(visit_number, ace, pd_report) |>
    pivot_wider(names_from = visit_number, values_from = pd_report) |>
    arrange(factor(ace, levels = c("double_orphan", "not_care", "food_insecure", "ipv_physical", "sex_violence", "ipv_sexual_violence", "school_violence")))
  
  
  
  # write to csv
  write_csv(adapt_report_pd, paste0(file_path, "Output Data/adapt_pd_results_", label, ".csv"))
  
  
  # round prevalence ratio
  adapt_results_ci$pr <- round(adapt_results_ci$prev_ratio, digits = 2)
  
  # round prevalence difference lower confidence limit
  adapt_results_ci$pr_lcl <- round(adapt_results_ci$prev_ratio_lower, digits = 2)
  
  # ound  prevalence difference upper confidence limit 
  adapt_results_ci$pr_ucl <- round(adapt_results_ci$prev_ratio_upper, digits = 2)
  
  # make nice format for table
  adapt_results_ci$pr_report <- paste0(adapt_results_ci$pr, " (", adapt_results_ci$pr_lcl, ", ", adapt_results_ci$pr_ucl, ")")
  
  
  # pivot wider just the `pr_report` column
  adapt_report_pr <- adapt_results_ci |>
    select(visit_number, ace, pr_report) |>
    pivot_wider(names_from = visit_number, values_from = pr_report) |>
    arrange(factor(ace, levels = c("double_orphan", "not_care", "food_insecure", "ipv_physical", "sex_violence", "ipv_sexual_violence", "school_violence")))
  
  # write to csv
  write_csv(adapt_report_pr, paste0(file_path, "Output Data/adapt_pr_results_", label, ".csv"))
  
}

