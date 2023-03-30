### PURPOSE OF THIS SCRIPT
## Analyze pulse-chase data where pulse is such that fn != 1 for all transcripts

# Load dependencies ------------------------------------------------------------
library(dplyr)
library(bakR)
library(purrr)


# Helper functions that I will use on multiple occasions
logit <- function(x) log(x/(1-x))
inv_logit <- function(x) exp(x)/(1+exp(x))


# Function to analyze data -----------------------------------------------------

#' Helper function for Monte Carlo analysis strategy
MC_samples <- function(pulse_mean = 0, pulse_sd = 0.2,
                   chase_mean = -1, chase_sd = 0.2,
                   tchase = 2, samps = 1000,
                   XF, Exp_ID){

  # Draw from logit(fs4U) distributions
  MC_lfn_pulse <- rnorm(samps, mean = pulse_mean,
                        sd = pulse_sd)

  MC_lfn_chase <- pmin(rnorm(samps, mean = chase_mean,
                             sd = chase_sd),
                       MC_lfn_pulse)

  # Make sure lfn_pulse > lfn_chase
  keeps <- which(MC_lfn_chase < MC_lfn_pulse)

  final_MC_lfn_pulse <- MC_lfn_pulse[keeps]
  final_MC_lfn_chase <- MC_lfn_chase[keeps]


  # Make sure there are enough valid samples
  redos <- 0

  while(length(keeps) < samps & redos < 10){

    MC_lfn_pulse <- rnorm(samps, mean = pulse_mean,
                          sd = pulse_sd)

    MC_lfn_chase <- pmin(rnorm(samps, mean = chase_mean,
                               sd = chase_sd),
                         MC_lfn_pulse)

    new_keeps <- which(MC_lfn_chase < MC_lfn_pulse)
    final_MC_lfn_pulse <- c(final_MC_lfn_pulse, MC_lfn_pulse[new_keeps])
    final_MC_lfn_chase <- c(final_MC_lfn_chase, MC_lfn_chase[new_keeps])

    keeps <- c(keeps, new_keeps)

    redos <- redos + 1

  }


  MC_kdeg <- -((log(inv_logit(final_MC_lfn_chase)) - log(inv_logit(final_MC_lfn_pulse)))/tchase)
  MC_log_kdeg <- log(MC_kdeg)

  kdeg <- mean(MC_kdeg)
  log_kdeg <- mean(MC_log_kdeg)

  kdeg_sd <- sd(MC_kdeg)
  log_kdeg_sd <- sd(MC_log_kdeg)


  out_df <- data.frame(Exp_ID = Exp_ID,
             XF = XF,
             kdeg = kdeg,
             log_kdeg = log_kdeg,
             kdeg_sd = kdeg_sd,
             log_kdeg_sd = log_kdeg_sd)


  return(out_df)

}

#' Documentation for function:
#'
#' INPUT:
#' fit_pulse: bakRFit object from analysis of 4sU pulse data. For this and all
#' other bakRFit objects given to the function as input, the NSS parameter in the
#' bakRFit() function should be set to TRUE.
#'
#' fit_chase_1: bakRFit object from analysis of chase time points from experimental
#' condition 1. metadf for this should assign a different Exp_ID for each
#' chase time (default option expects Exp_ID 1 to be chase of 1 hour, Exp_ID 2
#' to be a chase of 2 hours, Exp_ID 3 to be a chase of 4 hours, and Exp_ID 4 to
#' be a chase of 8 hours). !!!Do NOT!!! set chase = TRUE in bakRFit(). Only set
#' NSS = TRUE and all other parameters left as defaults.
#'
#' fit_chase_2: same as fit_chase_1 but for the second experimental condition.
#'
#' chase_dict: data frame that maps chase times to Exp_IDs in the fit_chase objects.
#' By default this is created such that Exp_ID 1 represents a chase of 1 hour, Exp_ID 2 represents
#' a chase of 2 hours, Exp_ID 3 represents a chase 4 hours, and Exp_ID 4 represents
#' a chase of 8 hours
#'
#' Hybrid: Boolean (TRUE or FALSE) indicating whether Fast_Fit or Hybrid_Fit
#' should be used for analysis. Default it to use Fast_Fit.
#'
#' ztest: Boolean indicating whether a z-test or moderated t-test will be used
#' for assessing differences in kdeg. z-test will be less conservative.
#'
#' conservative: Boolean indicating whether to reduce the standard deviation of
#' the L2FC(kdeg) standard errors. Setting this to FALSE can cause some very small
#' effect sizes to be called signficant.
#'
#' reg_factor: If conservative is TRUE, then this determines the extent to which the
#' L2FC(kdeg) standard error distribution has its standard deviation reduced.
#'
#' null_cutoff: |L2FC(kdeg)| cutoff for null hypothesis testing. The null hypothesis
#' tested will be |L2FC(kdeg)| < |null_cutoff|. Setting null_cutoff to 0 is the standard
#' t/z-test.
#'
#' Debug: Boolean that if TRUE will launch the environment browswer for debugging purposes
#'
#' MC: Boolean that if TRUE will use Monte Carlo to quantify uncertainty rather
#' than the less accurate delta approximation.
#'
#' samps: number of Monte Carlo samples to use
#'
#' OUTPUT:
#' A list with 4 objects in it:
#' 1) Fit: A fake bakRFit object that will allow you to make Volcano plots with
#' bakR's plotVolcano() function. The important estimates from this analysis (kdeg
#' estimates from each chase time and differences between experimental conditions),
#' are available as separate objects to be discussed next.
#' 2) Effects_df: Looks just like the Effects_df of a bakRFit object. Exp_IDs represent
#' each of the different chase times, and the L2FC(kdeg)'s and p values are for comparing
#' experimental condition 2 vs. 1 for that particular chase time.
#' 3) Chase_kdegs_E1: kdeg and log(kdeg) estimates and uncertainties from the
#' chase time points in experimental condition 1.
#' 4) Chase_kdegs_E2: kdeg and log(kdeg) estimates and uncertainties from the
#' chase time points in experimental condition 2.
#'
analyze_pulse_chase <- function(fit_pulse, fit_chase_1, fit_chase_2,
                                chase_dict = tibble(Exp_ID = 1:4,
                                                    tchase = c(1, 2, 4, 8)),
                                Hybrid = FALSE, ztest = FALSE, conservative = FALSE,
                                reg_factor = 2, null_cutoff = 0, Debug = FALSE,
                                MC = TRUE, samps = 1000){


  if(Debug){
    browser()
  }

  nreps <- max(fit_pulse$Fast_Fit$Fn_Estimates$Replicate)


  ### Extract Regularized estimates from the pulse and chase datasets
  avg_fits <- vector(mode = "list", length = 2)


  for(i in 1:2){

    if(Hybrid){
      # Extract logit(fraction new) estimates from Stan fit
      if(i == 1){
        reg_chase <- fit_chase_1$Hybrid_Fit$Fit_Summary

        lfns_chase <- reg_chase[grep("alpha", reg_chase$parameter),]

        # number of features, conditions, and replicates
        NF <- fit_chase_1$Data_lists$Stan_data$NF
        NC <- fit_chase_1$Data_lists$Stan_data$nMT


        essential_chase <- tibble(Exp_ID = rep(1:NC, times = NF),
                            lfn_sd_chase = lfns_chase$sd,
                            lfn_chase = lfns_chase$mean,
                            XF = fit_chase1$Fast_Fit$Regularized_ests$XF)

      }else{
        reg_chase <- fit_chase_2$Hybrid_Fit$Fit_Summary

        lfns_chase <- reg_chase[grep("alpha", reg_chase$parameter),]

        # number of features, conditions, and replicates
        NF <- fit_chase_2$Data_lists$Stan_data$NF
        NC <- fit_chase_2$Data_lists$Stan_data$nMT


        essential_chase <- tibble(Exp_ID = rep(1:NC, times = NF),
                                  lfn_sd_chase = lfns_chase$sd,
                                  lfn_chase = lfns_chase$mean,
                                  XF = fit_chase2$Fast_Fit$Regularized_ests$XF)
      }

      reg_pulse <- fit_pulse$Hybrid_Fit$Fit_Summary

      lfns_pulse <- reg_pulse[grep("alpha", reg_pulse$parameter),]

      # number of features, conditions, and replicates
      NF <- fit_pulse$Data_lists$Stan_data$NF
      NC <- fit_pulse$Data_lists$Stan_data$nMT


      essential_pulse <- tibble(Exp_ID = rep(1:NC, times = NF),
                                lfn_sd_pulse = lfns_pulse$sd,
                                lfn_pulse = lfns_pulse$mean,
                                XF = fit_pulse$Fast_Fit$Regularized_ests$XF)

      essential_pulse <- essential_pulse[essential_pulse$Exp_ID == i,c("lfn_sd_pulse", "lfn_pulse", "XF")]

    }else{
      if(i == 1){
        reg_chase <- fit_chase_1$Fast_Fit$Regularized_ests

      }else{
        reg_chase <- fit_chase_2$Fast_Fit$Regularized_ests
      }

      reg_pulse <- fit_pulse$Fast_Fit$Regularized_ests

      # Extract Exp_ID 1 from reg_pulse
      reg_pulse <- reg_pulse[reg_pulse$Exp_ID == i,]

      # Create simplified data frames that will be merged
      essential_chase <- reg_chase[,c("Exp_ID", "sd_post", "log_kdeg_post", "XF")]
      essential_pulse <- reg_pulse[,c("sd_post", "log_kdeg_post", "XF")]

      # Change column names to make more informative and facilitate merging
      colnames(essential_chase) <- c("Exp_ID", "lfn_sd_chase", "lfn_chase", "XF")
      colnames(essential_pulse) <- c("lfn_sd_pulse", "lfn_pulse", "XF")

    }


    ### Estimate kdeg and kdeg uncertainty

    # Merge pulse and chase data
    merged_fits <- inner_join(essential_chase, essential_pulse, by = c("XF"))

    # Estimate log(fn) for pulse and chase using delta approximation

    if(MC){

      ## initialize data.frame to add to merged_fits
      nr <- nrow(merged_fits)
      temp_data <- rep(0, times = nr)
      MC_df <- tibble(Exp_ID = temp_data,
                      XF = temp_data,
                      kdeg = temp_data,
                      log_kdeg = temp_data,
                      kdeg_sd = temp_data,
                      log_kdeg_sd = temp_data)

      # Add chase time information
      merged_fits <- dplyr::inner_join(merged_fits, chase_dict, by = "Exp_ID")

      # More efficient: pmap
      MC_list <- pmap(list(Exp_ID = merged_fits$Exp_ID,
                           XF = merged_fits$XF,
                           pulse_mean = merged_fits$lfn_pulse,
                           chase_mean = merged_fits$lfn_chase,
                           pulse_sd = merged_fits$lfn_sd_pulse,
                           chase_sd = merged_fits$lfn_sd_chase,
                           tchase = merged_fits$tchase,
                           samps = samps), MC_samples)


      MC_df <- bind_rows(MC_list)


      merged_fits <- inner_join(merged_fits, MC_df, by = c("Exp_ID", "XF"))

    }else{
      merged_fits <- merged_fits %>%
        dplyr::mutate(log_fn_chase = log(inv_logit(lfn_chase)),
                      log_fn_chase_sd = ((1 - inv_logit(lfn_chase))^2)*(lfn_sd_chase^2),
                      log_fn_pulse = log(inv_logit(lfn_pulse)),
                      log_fn_pulse_sd = ((1 - inv_logit(lfn_pulse))^2)*(lfn_sd_pulse^2))

      # Add chase time information
      merged_fits <- dplyr::inner_join(merged_fits, chase_dict, by = "Exp_ID")


      # Estimate kdeg and kdeg sd
      merged_fits <- merged_fits %>%
        dplyr::mutate(kdeg = (log_fn_pulse - log_fn_chase)/tchase,
                      kdeg_sd = sqrt( (log_fn_chase_sd)^2 + (log_fn_pulse_sd)^2 )/tchase)

    }



    # Filter out negative kdeg estimates (from when chase fn > pulse fn due to random variance)
    merged_fits <- merged_fits[merged_fits$kdeg > 0,]

    # Delta approximate log(kdeg) and log(kdeg) sds
    merged_fits <- merged_fits %>%
      dplyr::mutate(log_kdeg = log(kdeg),
             log_kdeg_sd = kdeg_sd/kdeg)

    # Keep necessary columns
    if(i == 1){
      avg_fits_1 <- merged_fits[,c("Exp_ID", "XF", "log_kdeg", "log_kdeg_sd", "kdeg", "kdeg_sd", "tchase")]

    }else{
      avg_fits_2 <- merged_fits[,c("Exp_ID", "XF", "log_kdeg", "log_kdeg_sd", "kdeg", "kdeg_sd", "tchase")]

    }


  }

  # Save chase kdeg estimates
  chase_estimates_E1 <- avg_fits_1
  chase_estimates_E2 <- avg_fits_2

  # Change column names so as to allow for inner joining and comparing
  colnames(avg_fits_1) <- c("Exp_ID", "XF", "log_kdeg_1", "log_kdeg_sd_1", "kdeg_1", "kdeg_sd_1")
  colnames(avg_fits_2) <- c("Exp_ID", "XF", "log_kdeg_2", "log_kdeg_sd_2", "kdeg_2", "kdeg_sd_2")

  compare_df <- inner_join(avg_fits_1, avg_fits_2, by = c("Exp_ID", "XF"))

  effects <- compare_df %>%
    dplyr::mutate(L2FC_kdeg = log_kdeg_2 - log_kdeg_1,
           effect = log_kdeg_2 - log_kdeg_1,
           se = sqrt(log_kdeg_sd_1^2 + log_kdeg_sd_2^2)) %>%
    dplyr::select(XF, Exp_ID, L2FC_kdeg, effect, se)

  if(conservative){

    chase_sd <- sd(effects$se)
    chase_mean <- mean(effects$se)

    # regularize chase ses
    effects$se <- (effects$se/reg_factor) + (chase_mean/reg_factor)
  }


  # add chase time to effect sizes
  if(ztest){
    effects <- inner_join(effects, chase_dict, by = "Exp_ID") %>%
      dplyr::mutate(pval = pmin(1, 2*pnorm((abs(effect) - null_cutoff)/se, lower.tail = FALSE)),
             padj = p.adjust(pval, method = "BH"))
  }else{
    effects <- inner_join(effects, chase_dict, by = "Exp_ID") %>%
      dplyr::mutate(pval = pmin(1, 2*pt((abs(effect) - null_cutoff)/se, df = 2*(nreps - 1) + 2*fit_pulse$Fast_Fit$Hyper_Parameters[1], lower.tail = FALSE)),
             padj = p.adjust(pval, method = "BH"))
  }



  # Create bakRFit object to facilitate visualizations
  Fit <- list(Data_lists = fit_pulse$Data_lists,
              Fast_Fit = fit_pulse$Fast_Fit)

  Fit$Fast_Fit$Effects_df <- effects
  class(Fit) <- c("bakRFit")

  # Merge output into list
  res <- list(Fit = Fit,
              Effects_df = effects,
              Chase_kdegs_E1 = chase_estimates_E1,
              Chase_kdegs_E2 = chase_estimates_E2)

}




################################################################################
######################### START OF TESTING #####################################
################################################################################



# Try out function on simulated data -------------------------------------------

####### SIMULATNG PULSE-CHASE DATA (not something bakR can do by itself) #######

# "pulse" (condition 2 much larger fraction news than condition 1)
sim_pulse <- Simulate_bakRData(500, nreps = 2,
                               eff_mean = 0, eff_sd = 1,
                               num_kd_DE = c(0, 500))

fit_pulse <- bakRFit(sim_pulse$bakRData, NSS = TRUE)


avg_pulse <- fit_pulse$Fast_Fit$Regularized_ests

tchase <- c(1, 2, 4, 8)
tpulse <- 4

avg_chase <- avg_pulse[,c("Feature_ID", "Exp_ID", "log_kdeg_post", "sd_post", "XF")]

avg_chase <- avg_chase[avg_chase$Exp_ID == 1,]

colnames(avg_chase) <- c("Feature_ID", "Exp_ID", "lfn_pulse", "lfn_pulse_sd", "XF")

# what is kdeg estimate?
avg_chase <- avg_chase %>%
  mutate(kdeg = -log(1-inv_logit(lfn_pulse))/tpulse)

count <- 1
lfn_chase <- c()
for(t in tchase){
  lfn_chase <- c(lfn_chase, rnorm(nrow(avg_chase),
                                  mean = logit(inv_logit(avg_chase$lfn_pulse)*exp(-avg_chase$kdeg*t)),
                                  sd = 0))
}

Feature_ID_sim <- rep(avg_chase$Feature_ID, times = length(tchase))
Exp_ID_sim <- rep(1:length(tchase), each = nrow(avg_chase))


# Create a simulated Regularized_ests for the chase
regest_chase <- tibble(Feature_ID = Feature_ID_sim,
                       Exp_ID = Exp_ID_sim,
                       avg_log_kdeg = 0,
                       sd_log_kdeg = 0,
                       nreads = 0,
                       sdp = 0,
                       theta_o = 0,
                       sd_post = rep(avg_chase$lfn_pulse_sd, times = length(tchase)),
                       log_kdeg_post = lfn_chase,
                       kdeg = 0,
                       kdeg_sd = 0,
                       XF = rep(avg_chase$XF, times = length(tchase)))


fit_chase_1 <- list(Fast_Fit = list(Regularized_ests = regest_chase))


avg_chase <- avg_pulse[,c("Feature_ID", "Exp_ID", "log_kdeg_post", "sd_post", "XF")]

avg_chase <- avg_chase[avg_chase$Exp_ID == 2,]

colnames(avg_chase) <- c("Feature_ID", "Exp_ID", "lfn_pulse", "lfn_pulse_sd", "XF")

# what is kdeg estimate?
avg_chase <- avg_chase %>%
  mutate(kdeg = -log(1-inv_logit(lfn_pulse))/tpulse)

count <- 1
lfn_chase <- c()
for(t in tchase){
  lfn_chase <- c(lfn_chase, rnorm(nrow(avg_chase),
                                  mean = logit(inv_logit(avg_chase$lfn_pulse)*exp(-avg_chase$kdeg*t)),
                                  sd = 0))
}

Feature_ID_sim <- rep(avg_chase$Feature_ID, times = length(tchase))
Exp_ID_sim <- rep(1:length(tchase), each = nrow(avg_chase))


# Create a simulated Regularized_ests for the chase
regest_chase <- tibble(Feature_ID = Feature_ID_sim,
                       Exp_ID = Exp_ID_sim,
                       avg_log_kdeg = 0,
                       sd_log_kdeg = 0,
                       nreads = 0,
                       sdp = 0,
                       theta_o = 0,
                       sd_post = rep(avg_chase$lfn_pulse_sd, times = length(tchase)),
                       log_kdeg_post = lfn_chase,
                       kdeg = 0,
                       kdeg_sd = 0,
                       XF = rep(avg_chase$XF, times = length(tchase)))


fit_chase_2 <- list(Fast_Fit = list(Regularized_ests = regest_chase))


######################## END OF SIMULATION CODE ################################


# Run analysis
results_MC <- analyze_pulse_chase(fit_pulse = fit_pulse,
                               fit_chase_1 = fit_chase_1,
                               fit_chase_2 = fit_chase_2,
                               ztest = FALSE,
                               conservative = FALSE,
                               MC = TRUE,
                               Debug = FALSE)

results <- analyze_pulse_chase(fit_pulse = fit_pulse,
                                  fit_chase_1 = fit_chase_1,
                                  fit_chase_2 = fit_chase_2,
                                  ztest = FALSE,
                                  conservative = TRUE,
                                  MC = FALSE)



# Are standard errors correlated between the two analysis strategies?
plot(results$Effects_df$se, results_MC$Effects_df$se)
abline(0,1)

plot(results$Effects_df$se[results$Effects_df$Exp_ID == 1], results_MC$Effects_df$se[results_MC$Effects_df$Exp_ID == 1])
abline(0,1)

plot(results$Effects_df$se[results$Effects_df$Exp_ID == 2], results_MC$Effects_df$se[results_MC$Effects_df$Exp_ID == 2])
abline(0,1)

plot(results$Effects_df$se[results$Effects_df$Exp_ID == 3], results_MC$Effects_df$se[results_MC$Effects_df$Exp_ID == 3])
abline(0,1)

plot(results$Effects_df$se[results$Effects_df$Exp_ID == 4], results_MC$Effects_df$se[results_MC$Effects_df$Exp_ID == 4])
abline(0,1)


# Make volcano plot
plotVolcano(results$Fit$Fast_Fit)
plotVolcano(results_MC$Fit$Fast_Fit)



### I want to make a volcano plot colored by reference kdeg.

Volcano_MC <- results_MC$Effects_df

pulse_reg <- fit_pulse$Fast_Fit$Regularized_ests

pulse_ests <- pulse_reg %>%
  group_by(XF) %>%
  summarise(reads = mean(nreads),
            lfn_se = mean(sd_post),
            lfn = mean(log_kdeg_post))


Volcano_MC <- inner_join(Volcano_MC, pulse_ests, by = "XF")


ggplot(Volcano_MC[order(-Volcano_MC$lfn),], aes(x = L2FC_kdeg, y = -log10(pval), color = lfn)) +
  geom_point() +
  scale_color_viridis_b() +
  xlab("L2FC(kdeg)") +
  ylab("-log10(p-value")




### Check that kdeg estimates agree
plot(results$Chase_kdegs_E1$log_kdeg, results_MC$Chase_kdegs_E1$log_kdeg)
abline(0,1)

plot(results$Chase_kdegs_E1$log_kdeg_sd, results_MC$Chase_kdegs_E1$log_kdeg_sd)
abline(0,1)

# How well does pulse estimate uncertainty correlate with L2FC(kdeg) uncertainty
pulse_reg <- fit_pulse$Fast_Fit$Regularized_ests

pulse_ests <- pulse_reg %>%
  group_by(XF) %>%
  summarise(reads = mean(nreads),
            lfn_se = mean(sd_post))

### Assess if correlation of uncertainty with read coverage is maintained

E1_pulse_ests <- pulse_reg[pulse_reg$Exp_ID == 1,] %>%
  dplyr::select(nreads, sd_post, XF)

E1_MC_ests <- results_MC$Chase_kdegs_E1
E1_MC_ests <- inner_join(E1_MC_ests, E1_pulse_ests, by = "XF")
E1_chase_ests <- results$Chase_kdegs_E1
E1_chase_ests <- inner_join(E1_chase_ests, E1_pulse_ests, by = "XF")


plot(log10(E1_chase_ests$nreads), E1_chase_ests$log_kdeg_sd,
     ylim = c(0, 0.5))
plot(log10(E1_MC_ests$nreads), E1_MC_ests$log_kdeg_sd)

plot(log10(E1_MC_ests$sd_post), E1_MC_ests$log_kdeg_sd)
plot(log10(E1_MC_ests$sd_post), E1_chase_ests$log_kdeg_sd,
     ylim = c(0, 0.5))


# Prep for MA plot
plot(log10(pulse_ests$reads), pulse_ests$lfn_se)


Effects_MC <- results_MC$Effects_df
Effects <- results$Effects_df

Effects_MC <- inner_join(Effects_MC, pulse_ests, by = "XF")
Effects <- inner_join(Effects, pulse_ests, by = "XF")

plot(Effects_MC$lfn_se, Effects_MC$se)
plot(Effects$lfn_se, Effects$se, ylim = c(0, 0.5))

# Make MA plot
Effects_MC <- Effects_MC %>%
  mutate(sig = as.factor(ifelse(pval < 0.01, ifelse(effect < 0, "Stabilized", "Destabilized"), "Not Sig.")))

Effects <- Effects %>%
  mutate(sig = as.factor(ifelse(pval < 0.01, ifelse(effect < 0, "Stabilized", "Destabilized"), "Not Sig.")))


library(ggplot2)

ggplot(Effects_MC, aes(x = log10(reads), y = L2FC_kdeg, color = sig)) +
  geom_point() +
  xlab("log10(reads)") +
  ylab("L2FC(kdeg)") +
  scale_color_manual(values = c("orange", "darkgray"))

ggplot(Effects, aes(x = log10(reads), y = L2FC_kdeg, color = sig)) +
  geom_point() +
  xlab("log10(reads)") +
  ylab("L2FC(kdeg)") +
  scale_color_manual(values = c("orange", "darkgray", "blue"))



##### TEST OUT MC STRATEGY -----------------------------------------------------

MC_sim <- function(pulse_mean = 0, pulse_sd = 0.2,
                   chase_mean = -1, chase_sd = 0.2,
                   tchase = 2, samps = 1000){

  MC_lfn_pulse <- rnorm(samps, mean = pulse_mean,
                        sd = pulse_sd)

  MC_lfn_chase <- pmin(rnorm(samps, mean = chase_mean,
                             sd = chase_sd),
                       MC_lfn_pulse)

  # Make sure lfn_pulse > lfn_chase
  keeps <- which(MC_lfn_chase < MC_lfn_pulse)

  final_MC_lfn_pulse <- MC_lfn_pulse[keeps]
  final_MC_lfn_chase <- MC_lfn_chase[keeps]


  # Make sure there are enough valid samples
  redos <- 0

  while(length(keeps) < 1000 & redos < 10){

    MC_lfn_pulse <- rnorm(samps, mean = pulse_mean,
                          sd = pulse_sd)

    MC_lfn_chase <- pmin(rnorm(samps, mean = chase_mean,
                               sd = chase_sd),
                         MC_lfn_pulse)

    new_keeps <- which(MC_lfn_chase < MC_lfn_pulse)
    final_MC_lfn_pulse <- c(final_MC_lfn_pulse, MC_lfn_pulse[new_keeps])
    final_MC_lfn_chase <- c(final_MC_lfn_chase, MC_lfn_chase[new_keeps])

    keeps <- c(keeps, new_keeps)

    redos <- redos + 1

  }






  MC_kdeg <- -((log(inv_logit(final_MC_lfn_chase)) - log(inv_logit(final_MC_lfn_pulse)))/tchase)
  MC_log_kdeg <- log(MC_kdeg)

  kdeg <- mean(MC_kdeg)
  log_kdeg <- mean(MC_log_kdeg)

  kdeg_sd <- sd(MC_kdeg)
  log_kdeg_sd <- sd(MC_log_kdeg)

  return(log_kdeg_sd)
}

### Let's try varying the standard deviations
sds <- seq(from = 0.05, to = 2, length.out = 30)

log_kdeg_sd <- rep(0, times = length(sds))
for(i in 1:length(sds)){
  log_kdeg_sd[i] <- MC_sim(pulse_sd = sds[i], chase_sd = sds[i], samps =  10000)
}

plot(sds, log_kdeg_sd)
abline(0,1)

### Let's try varying the chase mean
means <- seq(from = -10, to = -1, length.out = 30)

log_kdeg_sd <- rep(0, times = length(means))
for(i in 1:length(sds)){
  log_kdeg_sd[i] <- MC_sim(pulse_mean = -0.5, chase_mean = means[i], samps =  10000)
}

plot(means, log_kdeg_sd)

### What if I vary the L2FC(kdeg)? So fix pulse-chase fn and vary chase accordingly?

# Reference kdeg
kdeg_E1 <- 0.25

# chase time
tchase <- 2

# pulse time
tpulse <- 4

# L2FC(kdeg)
L2FC_kdeg <- seq(from = -2, to = 2, length.out = 30)

# Experimental kdeg
kdeg_E2 <- 2^(L2FC_kdeg + log2(kdeg_E1))

# Reference fn after pulse
fn_p_E1 <- 1 - exp(-kdeg_E1*tpulse)

# Experimental fn after pulse
fn_p_E2 <- 1 - exp(-kdeg_E2*tpulse)


# Reference fn after chase
fn_c_E1 <- fn_p_E1*exp(-kdeg_E1*tchase)

# Experimental fn after chase
fn_c_E2 <- fn_p_E2*exp(-kdeg_E2*tchase)


log_kdeg_sd_E1 <- MC_sim(pulse_mean = logit(fn_p_E1),
                         chase_mean = logit(fn_c_E1))

log_kdeg_sd_E2 <- rep(0, times = length(L2FC_kdeg))


for(i in 1:length(L2FC_kdeg)){

  log_kdeg_sd_E2[i] <- MC_sim(pulse_mean = logit(fn_p_E2[i]),
                           chase_mean = logit(fn_c_E2[i]),
                           samps =  1000)
}


plot(L2FC_kdeg, log_kdeg_sd_E2)


##### WHAT IF WE ONLY CONSIDER KDEGS ON NATURAL SCALES? ------------------------

MC_sim <- function(pulse_mean = 0, pulse_sd = 0.2,
                   chase_mean = -1, chase_sd = 0.2,
                   tchase = 2, samps = 1000){

  MC_lfn_pulse <- rnorm(samps, mean = pulse_mean,
                        sd = pulse_sd)

  MC_lfn_chase <- pmin(rnorm(samps, mean = chase_mean,
                             sd = chase_sd),
                       MC_lfn_pulse)

  # Make sure lfn_pulse > lfn_chase
  keeps <- which(MC_lfn_chase < MC_lfn_pulse)

  final_MC_lfn_pulse <- MC_lfn_pulse[keeps]
  final_MC_lfn_chase <- MC_lfn_chase[keeps]


  # Make sure there are enough valid samples
  redos <- 0

  while(length(keeps) < 1000 & redos < 10){

    MC_lfn_pulse <- rnorm(samps, mean = pulse_mean,
                          sd = pulse_sd)

    MC_lfn_chase <- pmin(rnorm(samps, mean = chase_mean,
                               sd = chase_sd),
                         MC_lfn_pulse)

    new_keeps <- which(MC_lfn_chase < MC_lfn_pulse)
    final_MC_lfn_pulse <- c(final_MC_lfn_pulse, MC_lfn_pulse[new_keeps])
    final_MC_lfn_chase <- c(final_MC_lfn_chase, MC_lfn_chase[new_keeps])

    keeps <- c(keeps, new_keeps)

    redos <- redos + 1

  }






  MC_kdeg <- -((log(inv_logit(final_MC_lfn_chase)) - log(inv_logit(final_MC_lfn_pulse)))/tchase)
  MC_log_kdeg <- log(MC_kdeg)

  kdeg <- mean(MC_kdeg)
  log_kdeg <- mean(MC_log_kdeg)

  kdeg_sd <- sd(MC_kdeg)
  log_kdeg_sd <- sd(MC_log_kdeg)

  return(kdeg_sd)
}

### Let's try varying the standard deviations
sds <- seq(from = 0.05, to = 2, length.out = 30)

log_kdeg_sd <- rep(0, times = length(sds))
for(i in 1:length(sds)){
  log_kdeg_sd[i] <- MC_sim(pulse_sd = sds[i], chase_sd = sds[i], samps =  10000)
}

plot(sds, log_kdeg_sd)
abline(0,1)

### Let's try varying the chase mean
means <- seq(from = -3, to = -1, length.out = 30)

log_kdeg_sd <- rep(0, times = length(means))
for(i in 1:length(sds)){
  log_kdeg_sd[i] <- MC_sim(pulse_mean = -0.5, chase_mean = means[i], samps =  10000)
}

plot(means, log_kdeg_sd)

### What if I vary the L2FC(kdeg)? So fix pulse-chase fn and vary chase accordingly?

# Reference kdeg
kdeg_E1 <- 0.25

# chase time
tchase <- 2

# pulse time
tpulse <- 4

# L2FC(kdeg)
L2FC_kdeg <- seq(from = -2, to = 2, length.out = 30)

# Experimental kdeg
kdeg_E2 <- 2^(L2FC_kdeg + log2(kdeg_E1))

# Reference fn after pulse
fn_p_E1 <- 1 - exp(-kdeg_E1*tpulse)

# Experimental fn after pulse
fn_p_E2 <- 1 - exp(-kdeg_E2*tpulse)


# Reference fn after chase
fn_c_E1 <- fn_p_E1*exp(-kdeg_E1*tchase)

# Experimental fn after chase
fn_c_E2 <- fn_p_E2*exp(-kdeg_E2*tchase)


log_kdeg_sd_E1 <- MC_sim(pulse_mean = logit(fn_p_E1),
                         chase_mean = logit(fn_c_E1))

log_kdeg_sd_E2 <- rep(0, times = length(L2FC_kdeg))


for(i in 1:length(L2FC_kdeg)){

  log_kdeg_sd_E2[i] <- MC_sim(pulse_mean = logit(fn_p_E2[i]),
                              chase_mean = logit(fn_c_E2[i]),
                              samps =  1000)
}


plot(L2FC_kdeg, log_kdeg_sd_E2)



##### EXPLORE CORRELATIONS -----------------------------------------------------

MC_corr <- function(pulse_mean = 0, pulse_sd = 0.2,
                   chase_mean = -0.3, chase_sd = 0.2,
                   tchase = 2, samps = 1000){

  MC_lfn_pulse <- rnorm(samps, mean = pulse_mean,
                        sd = pulse_sd)

  MC_lfn_chase <- pmin(rnorm(samps, mean = chase_mean,
                             sd = chase_sd),
                       MC_lfn_pulse)

  # Make sure lfn_pulse > lfn_chase
  keeps <- which(MC_lfn_chase < MC_lfn_pulse)

  final_MC_lfn_pulse <- MC_lfn_pulse[keeps]
  final_MC_lfn_chase <- MC_lfn_chase[keeps]


  # Make sure there are enough valid samples
  redos <- 0

  while(length(keeps) < 1000 & redos < 10){

    MC_lfn_pulse <- rnorm(samps, mean = pulse_mean,
                          sd = pulse_sd)

    MC_lfn_chase <- pmin(rnorm(samps, mean = chase_mean,
                               sd = chase_sd),
                         MC_lfn_pulse)

    new_keeps <- which(MC_lfn_chase < MC_lfn_pulse)
    final_MC_lfn_pulse <- c(final_MC_lfn_pulse, MC_lfn_pulse[new_keeps])
    final_MC_lfn_chase <- c(final_MC_lfn_chase, MC_lfn_chase[new_keeps])

    keeps <- c(keeps, new_keeps)

    redos <- redos + 1

  }

  MC_df <- tibble(lfn_pulse = final_MC_lfn_pulse,
                  lfn_chase = final_MC_lfn_chase)

}

MC_data <- MC_corr()


plot(MC_data$lfn_pulse, MC_data$lfn_chase)
abline(0,1)




# # Run Hybrid analysis
# results <- analyze_pulse_chase(fit_pulse = fit_pulse,
#                                fit_chase_1 = fit_chase_1,
#                                fit_chase_2 = fit_chase_2,
#                                Hybrid = TRUE)
#
#
# # Make Hybrid volcano plot
# plotVolcano(results$Fit$Fast_Fit)
#
#
#
# ### Simulate independent data to make sure it works on actual bakRFits
# sim_pulse <- Simulate_bakRData(500, nreps = 2,
#                                num_kd_DE = c(0, 0),
#                                fn_mean = 1, fn_sd = 0.05)
#
# sim_chase_1 <- Simulate_bakRData(500, nreps = 2,
#                                num_kd_DE = c(0, 500),
#                                eff_mean = -0.75,
#                                eff_sd = 0.05,
#                                fn_mean = 0, fn_sd = 0.05)
#
# sim_chase_2 <- Simulate_bakRData(500, nreps = 2,
#                                  num_kd_DE = c(0, 500),
#                                  eff_mean = -0.75,
#                                  eff_sd = 0.05,
#                                  fn_mean = 0, fn_sd = 0.05)
#
#
# fit_pulse <- bakRFit(sim_pulse$bakRData, NSS = TRUE)
# fit_chase_1 <- bakRFit(sim_chase_1$bakRData, NSS = TRUE)
# fit_chase_2 <- bakRFit(sim_chase_2$bakRData, NSS = TRUE)
#
#
# # analyze
# results <- analyze_pulse_chase(fit_pulse = fit_pulse,
#                                fit_chase_1 = fit_chase_1,
#                                fit_chase_2 = fit_chase_2,
#                                chase_dict = data.frame(Exp_ID = 1:2, tchase = c(2, 4)),
#                                ztest = TRUE)
#
# plotVolcano(results$Fit$Fast_Fit)
#
# ### Hybrid
#
# # Get sampling warnings usually because data is weird, but that's to be expected.
# fit_pulse <- bakRFit(fit_pulse, HybridFit = TRUE, NSS = TRUE)
# fit_chase_1 <- bakRFit(fit_chase_1, HybridFit = TRUE, NSS = TRUE)
# fit_chase_2 <- bakRFit(fit_chase_2, HybridFit = TRUE, NSS = TRUE)
#
# # analyze
# results <- analyze_pulse_chase(fit_pulse = fit_pulse,
#                                fit_chase_1 = fit_chase_1,
#                                fit_chase_2 = fit_chase_2,
#                                chase_dict = data.frame(Exp_ID = 1:2, tchase = c(2, 4)),
#                                Hybrid = TRUE)
#
