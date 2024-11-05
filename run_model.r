# Start #####################################################################
# Title: Discrete Event Simulation (DES) model for Economic Evaluation in Schizophrenia
# Description: This is the script including all the key codes of the model to generate teh base-case evaluation results shown in the manuscript

# Initialization ----

# Clean
rm(list = ls())

# Set path (Need to set the path where the script and data are stored)
# setwd("J:\\DES Schizo Main\\Illustration")
path_wd <- getwd()

# Package
library(tidyverse)
library(foreach)
library(snowfall)
library(doSNOW)
library(parallel)
library(doParallel)

# Function ----

f_prepare_input <- function(para){
  
  # para <- readRDS(file.path(path_output, "input_base.rds"))
  # para <- input_base
  
  # Available treatment ----
  
  name_trt0 <- c("hal", 
                 "ami", "ari", "car", "lur","ola", "pal", "que", "ris",
                 "clo",
                 "ari2", "pal2", "ris2")
  
  name_trt <- name_trt0[str_c("discon_", name_trt0, "_prob") %in% names(para)]
  
  # Treatment comparative effect ----
  
  name_cf_mbse <- c("bmi", "tch", "hdl", "tri", "sbp", "glu")
  
  cf_trt <- list()
  for(trt_i in name_trt){
    # trt_i <- name_trt[1]
    
    tmp_cf_mbse_md <- para[str_c(name_cf_mbse, "_md_", trt_i, "_pla")] %>% unlist()
    
    tmp_cf_ltae_prob <- para[str_c(c("td", "discon"), "_", trt_i, "_prob")]
    tmp_cf_ltae_rate <- map(tmp_cf_ltae_prob, ~(-log(1-.x)/1)) %>% unlist()
    
    tmp_cf_rlps_rr_pla <- para[str_c("rlps_rr_", trt_i, "_pla")] %>% unlist()
    tmp_cf_rlps_prob <- para$rlps_prob
    tmp_cf_rlps_ar_pla <- 1/(log(1- min(0.99, tmp_cf_rlps_prob* tmp_cf_rlps_rr_pla)) / log(1- tmp_cf_rlps_prob))
    tmp_cf_rlps_adh_ar_pla <- tmp_cf_rlps_ar_pla * c(1, 1/para$rlps_hr_cmp_part_full, 1/para$rlps_hr_cmp_none_full)
    
    tmp_cf <- c(tmp_cf_mbse_md, tmp_cf_ltae_rate, tmp_cf_rlps_adh_ar_pla)
    names(tmp_cf) <- c(name_cf_mbse, "td", "discon", "rlps_full", "rlps_part", "rlps_none")
    if(trt_i == "clo"){
      tmp_cf <- c(tmp_cf, "agc" = para$agc_clo_prop)
    }
    
    # Fully adherent when the patients are on depot treatment
    if(str_detect(trt_i, "2")){
      tmp_cf[["rlps_part"]] <- tmp_cf[["rlps_full"]] 
      tmp_cf[["rlps_none"]] <- tmp_cf[["rlps_full"]] 
    }
    
    cf_trt[[trt_i]] <- tmp_cf
  }
  
  # Survival equation (baseline) ----
  
  cf_b <- with(para, 
               list(rcvr = rcvr_time,
                    rlps = -log(1-rlps_prob),
                    switch = switch_rlps_prop,
                    db = c("int" = db_logit_intercept, 
                           log(c("age_50t64" = db_or_age_50t64,
                                 "age_65" = db_or_age_65, 
                                 "male" = db_or_male, 
                                 "pardb" = db_or_pardb))),
                    chd = list(m1 = c("int" = chd1_m_aft_weibull_intercept, 
                                      log(c("age" = chd1_m_ar_age,
                                            "smk" = chd1_m_ar_smk)),
                                      "scale" = chd1_m_aft_weibull_scale),
                               m2 = c("int" = chd2_m_aft_weibull_intercept, 
                                      "age" = log(chd2_m_ar_age),
                                      "scale" = chd2_m_aft_weibull_scale),
                               f1 = c("int" = chd1_f_aft_weibull_intercept, 
                                      log(c("age" = chd1_f_ar_age,
                                            "meno" = chd1_f_ar_meno,
                                            "age_int_meno" = chd1_f_ar_age_int_meno,
                                            "smk" =  chd1_f_ar_smk,
                                            "alco" = chd1_f_ar_alcohol)),
                                      "scale" = chd1_f_aft_weibull_scale),
                               f2 = c("int" = chd2_f_aft_weibull_intercept, 
                                      log(c("age" = chd2_f_ar_age,
                                            "smk" = chd2_f_ar_smk)),
                                      "scale" = chd2_f_aft_weibull_scale)
                    ),
                    strk = list(m = log(c("int" = strk_m_aft_exp_rate,
                                          "age_sub_50" = strk_m_hr_age, 
                                          "lvh" = strk_m_hr_lvh, 
                                          "smk" = strk_m_hr_smk, 
                                          "af" = strk_m_hr_af)),
                                f = log(c("int" = strk_f_aft_exp_rate,
                                          "age_sub_50" = strk_f_hr_age, 
                                          "lvh" = strk_f_hr_lvh, 
                                          "smk" = strk_f_hr_smk, 
                                          "af" = strk_f_hr_af))
                    ),
                    death = list(m = lifetable %>% 
                                   select(age, male) %>%
                                   mutate(prob = pmin(1, male * death_m_smr)) %>%  
                                   select(-male),
                                 f = lifetable %>% 
                                   select(age, female) %>%
                                   mutate(prob = pmin(1, female * death_f_smr)) %>% 
                                   select(-female)),
                    death_agc = death_agc_prop, 
                    death_chd = list(f = death_chd_f_prop, 
                                     m = death_chd_m_prop),
                    death_strk = list(f = death_strk_f_prop, 
                                      m = death_strk_m_prop)
               )
  )
  
  # Acceleration ratios (Time varying) ----
  
  cf_t <- list()
  
  cf_t$db <- with(para, 
                  tibble(or = c(db_or_bmi_25t30, 
                                db_or_bmi_30,
                                db_or_sbp_int_atht,
                                db_or_hdl_int_sex,
                                db_or_tri_150,
                                db_or_glu_100),
                         odd_base = exp(db_logit_intercept)) %>% 
                    mutate(odd = or * odd_base,
                           prob = odd / (1 + odd),
                           rate = -log(1 - prob) / 1 ,
                           hr = rate / (-log( 1 - (odd_base / (1 + odd_base)) ) / 1),
                           logar = log(1 / hr)) %>% 
                    pull(logar))
  names(cf_t$db) <- c("bmi_25t30", "bmi_30", "sbp_int_atht", "hdl_int_sex","tri_150", "glu_100")
  
  cf_t$chd <- with(para, 
                   list(m1 = log(c(chd1_m_ar_tch_int_hdl, 
                                   chd1_m_ar_sbp,
                                   chd1_m_ar_sbp_int_atht,
                                   chd1_m_ar_db)),
                        f1 = log(c(chd1_f_ar_tch_int_hdl, 
                                   chd1_f_ar_sbp,
                                   chd1_f_ar_sbp_int_atht,
                                   chd1_f_ar_db,
                                   chd1_f_ar_tri)),
                        m2 = log(c(chd2_m_ar_tch_int_hdl, 
                                   chd2_m_ar_db)),
                        f2 = log(c(chd2_f_ar_tch_int_hdl, 
                                   chd2_f_ar_sbp,
                                   chd2_f_ar_db)))
  )
  names(cf_t$chd$m1) <- c("tch_int_hdl", "lnsbp", "sbp_int_atht2", "db")
  names(cf_t$chd$f1) <- c("tch_int_hdl", "lnsbp", "sbp_int_atht2", "db", "lntri")
  names(cf_t$chd$m2) <- c("tch_int_hdl", "db")
  names(cf_t$chd$f2) <- c("tch_int_hdl", "lnsbp", "db")
  
  cf_t$strk <- with(para,
                    list(m = log(1/c(strk_m_hr_sbp, 
                                     strk_m_hr_sbp_int_atht,
                                     strk_m_hr_cvd,
                                     strk_m_hr_db)),
                         f = log(1/c(strk_f_hr_sbp, 
                                     strk_f_hr_sbp_int_atht,
                                     strk_f_hr_cvd,
                                     strk_f_hr_db))
                    )
  )
  names(cf_t$strk$m) <- c("sbp_ctr", "sbp_int_atht2", "cvd", "db")
  names(cf_t$strk$f) <- c("sbp_ctr", "sbp_int_atht2", "cvd", "db")
  
  # Others ----
  
  # > QOL 
  
  qol_all <- with(para, 
                  list(cf = c("rcvr" = qol_rcvr, 
                              "rlps" = qol_rlps, 
                              "td" = qol_md_td, 
                              "db" = qol_md_db, 
                              "chd" = qol_md_chd, 
                              "strk" = qol_md_strk),
                       agc = (-1) * (1 - qol_rd_agc))
  )
  
  name_stse <- c("eps", "wg", "sdtn", "sxdf")
  tmp_cf <- with(para, c(qol_md_eps, qol_md_sdtn, qol_md_sxdf, qol_md_wg) * 3 / 12)
  
  qol_all$stse <- map(name_trt %>% set_names(), 
                      ~tmp_cf %*% (para[str_c(name_stse, "_", .x, "_prop")] %>% unlist())) %>% 
    unlist()
  
  # > Cost 
  
  cost_all <- with(para,
                   list(ru = c("rcvr" = cost_y_rcvr,
                               "rlps" = cost_y_rlps,
                               "rlps_acute" = cost_rlps),
                        agc = cost_agc,
                        td = cost_td,
                        chd = list(cvd0 = c("y0" = cost_y_chd_same_cvd0,
                                            "y1" = cost_y_chd_1y_cvd0,
                                            "y2" = cost_y_chd_2y_cvd0,
                                            "y3" = cost_y_chd_ge3y_cvd0),
                                   cvd1 = c("y0" = cost_y_chd_same_cvd1,
                                            "y1" = cost_y_chd_1y_cvd1,
                                            "y2" = cost_y_chd_2y_cvd1,
                                            "y3" = cost_y_chd_ge3y_cvd1)
                        ),
                        strk = list(cvd0 = c("y0" = cost_y_strk_same_cvd0,
                                             "y1" = cost_y_strk_1y_cvd0,
                                             "y2" = cost_y_strk_2y_cvd0,
                                             "y3" = cost_y_strk_ge3y_cvd0),
                                    cvd1 = c("y0" = cost_y_strk_same_cvd1,
                                             "y1" = cost_y_strk_1y_cvd1,
                                             "y2" = cost_y_strk_2y_cvd1,
                                             "y3" = cost_y_strk_ge3y_cvd1)
                        ),
                        db = list(cvd0 = c("y0_10" = cost_y_db_lt10y_cvd0,
                                           "y10p" = cost_y_db_ge10y_cvd0),
                                  cvd1 = c("y0_10" = cost_y_db_lt10y_cvd1,
                                           "y10p" = cost_y_db_ge10y_cvd1)
                        ),
                        death = list(cvd0 = cost_nvd_cvd0,
                                     cvd1 = cost_nvd_cvd1),
                        death_cv = list(cvd0 = cost_vd_cvd0,
                                        cvd1 = cost_vd_cvd1)
                   )
  )
  
  name_stse <- c("eps", "wg", "sdtn", "sxdf")
  tmp_cf <- with(para, c(cost_eps, cost_wg, cost_sdtn, cost_sxdf))
  cost_all$stse <- map(name_trt %>% set_names(), 
                       ~tmp_cf %*% (para[str_c(name_stse, "_", .x, "_prop")] %>% unlist())) %>% 
    unlist()
  
  cost_all$trt <- map_dbl(name_trt %>% set_names(), 
                          ~unlist(para[str_c("cost_y_", .x)],use.names = FALSE))
  
  output <- list(cf_trt = cf_trt, 
                 cf_b = cf_b,
                 cf_t = cf_t,
                 qol_all = qol_all,
                 cost_all = cost_all)
  return(output)
  
}
f_disc_con <- function(val, disc, t1, t2){
  
  # For compound continuous discounting - use INSTANTANEOUS rate.
  if(disc > 0){
    disc2 <- log(1 + disc)
    output <- ((val)/(0 - disc2)) * (exp(t2 * (0 - disc2)) - exp(t1 * (0 - disc2)))
  } else {
    output <- val * (t2 - t1)
  }
  return(output)
}
f_disc_ins <- function(val, disc, t){
  
  #discrete time discoutning for instantaneous costs and benefits
  output <- val * ((1+disc)^(-t)) 
  return(output)
}
f_simulate_progress <- function(pats, para, scn, mod, nsim = 10, ncore = 2, opt_ipd = FALSE){
# 
  # scn <- c("ari", "ari2", "clo")
  # mod <- list(tf = expr(100 - pat_i[['age']]), disc_cost = 0.025, disc_qaly = 0.025)
  # para <- readRDS(file.path(path_output, "para_base.rds"))
  # pats <- readRDS(file.path(path_output, "mod_pat.rds"))
  # nsim <- 10
  # ncore <- 2
  
  pats[,"alco"] <- pats[,"alco"] * 0.282 # convert uk unit to us ounces for using one Framingham equation for CHD prediction for female
  
  # Prepare global parameter ----
  
  cf_trts_all <- para$cf_trt
  cf_b_all <- para$cf_b
  cf_t_all <- para$cf_t
  qol_all <- para$qol_all
  cost_all <- para$cost_all
  
  # Prepare initial outcome 
  out_all <- list(
    qol = map_dbl(c("tot", "tot_disc", "dis", "stse", "ltse") %>% set_names(), ~0),
    cost = map_dbl(c("tot", "tot_disc", "drug", "ru","ltse", "death") %>% set_names(), ~0),
    other = map_dbl(c("dur_rcvr", "n_rlps", "chd10", "db10", "strk10") %>% set_names(), ~0)
  )
  
  # Create number of lines to stop switching
  nline_all <- length(scn)
  
  # update treatment impact 
  tmp_trt <- scn[[1]]
  tmp_cf_trt <- cf_trts_all[[tmp_trt]]
  name_mbse <- c("bmi", "tch", "hdl", "tri", "sbp", "glu")
  for(mbse in name_mbse){
    pats[,mbse] <- pats[,mbse] + tmp_cf_trt[mbse]
  }
  
  # Short-term side effect after receiving treatment 
  # - Never start from clozapine, so no agranulocytosis
  out_all$qol[c("tot","tot_disc","stse")] <- qol_all$stse[[tmp_trt]]
  out_all$cost[c("tot","tot_disc","stse")] <- cost_all$stse[[tmp_trt]]
  
  # Simulate each individual ----
  npat <- nrow(pats)
  
  set.seed(1234)
  random_seed <- matrix(sample(1:(npat * nsim), replace = FALSE), nrow = npat)
 
  #   # initiate parallel
  cl <- makeCluster(ncore)
  clusterExport(cl, c("f_disc_con", "f_disc_ins")) 
  clusterEvalQ(cl, library(tidyverse))
  registerDoParallel(cl)

  # output <- rep(list(NULL), nrow(pats))

  # for(i in 1:nrow(pats)){
  
  output <- foreach (i = 1:nrow(pats)) %dopar% {
    
    ## Prepare individual parameter ----
    
    # i = 1
    pat_i <- pats[i,]
    tf_pat <- eval(mod$tf) # for lifetime simulation, tf is related to the age of the patient
    
    ### Treatment ----
    
    # Only related to baseline compliance level
    tmp_cmp <- ifelse(pat_i[["cmp_Full"]] == 1, 1, ifelse(pat_i[["cmp_Partial"]] == 1, 2, 3))
    cf_trts_pat <- map(cf_trts_all,
                       ~c(.x[setdiff(names(.x), c("rlps_full", "rlps_part", "rlps_none"))], 
                          "rlps" = .x[c("rlps_full", "rlps_part", "rlps_none")][[tmp_cmp]]))
    
    ### Event ----
    
    cf_b_pat <- cf_b_all
    cf_t_pat <- cf_t_all
    
    # Create temporal patient profiles for profiles extraction
    tmp_age <- pat_i[["age"]]
    tmp_male <- pat_i[["male"]]
    tmp_alco <- pat_i[["alco"]]
    tmp_smk <- pat_i[["smk"]]
    tmp_af <- pat_i[["af"]]
    tmp_chd <- pat_i[["chd"]]
    tmp_cvd <- ifelse(pat_i[["strk"]] == 1|tmp_chd == 1,1,0)
    tmp_lvh<- pat_i[["lvh"]]
    tmp_pardb <- pat_i[["par_db"]]
    
    # Diabetes 
    tmp_cf <- (cf_b_all$db %*% c("int" = 1, 
                                 "age_50t64" = ifelse(tmp_age >= 50 & tmp_age < 65, 1,0),
                                 "age_65" = ifelse(tmp_age >= 65, 1,0),
                                 "male" = tmp_male,
                                 "pardb" = tmp_pardb)) %>% 
      as.vector() %>% exp()
    cf_b_pat$db <- -log(1 - (tmp_cf / (tmp_cf + 1))) / 1 
    # Odd to rate
    # - Odd to prob: prob = odd / (1 + odd)
    # - Prob to Rate (assuming constant rate): rate = -ln(1-prob) / t
    
    cf_t_pat$db <- expr(
      (cf_t_all$db %*% c("bmi_25t30" = ifelse(pat[["bmi"]] >= 25 & pat[["bmi"]] < 30, 1, 0),
                         "bmi_30" = ifelse(pat[["bmi"]] >= 30, 1,0),
                         "sbp_int_atht" = ifelse(pat[["atht"]] == 1 | pat[["sbp"]] >= 130, 1, 0),
                         "hdl_int_sex" = ifelse(pat[["male"]] == 1, ifelse(pat[["hdl"]] < 40, 1, 0),
                                                ifelse(pat[["hdl"]] < 50, 1, 0)),
                         "tri_150" = ifelse(pat[["tri"]] >= 150, 1, 0) ,
                         "glu_100" = ifelse(pat[["glu"]] >= 100, 1, 0))) %>%
        as.vector() %>% exp() 
    )
    
    # CHD 
    if(tmp_male == 1 & tmp_chd == 0){
      tmp_cf <- (cf_b_all$chd$m1[1:3] %*% c("int" = 1, 
                                            "age" = tmp_age,
                                            "smk" = tmp_smk)) %>% 
        as.vector() %>% exp()
      tmp_scale <- cf_b_all$chd$m1[[4]]
      cf_b_pat$chd <- c(tmp_cf, tmp_scale)
      cf_t_pat$chd <- expr(
        (cf_t_all$chd$m1 %*% c("tch_int_hdl" = log(max(0.01,pat[["tch"]]/pat[["hdl"]])),
                               "lnsbp" = log(max(0.01,pat[["sbp"]])),
                               "sbp_int_atht2" = ifelse(between(pat[["sbp"]], 110, 200) & pat[["atht"]] == 1, 
                                                        (200-pat[["sbp"]]) * (pat[["sbp"]] - 110) / 100, 0),
                               "db" = pat[["db"]])) %>% 
          as.vector() %>% exp()
      )
    } else if (tmp_male == 1 & tmp_chd == 1){
      tmp_cf <- (cf_b_all$chd$m2[1:2] %*% c("int" = 1, 
                                            "age" = tmp_age)) %>% 
        as.vector() %>% exp()
      tmp_scale <- cf_b_all$chd$m2[[3]]
      cf_b_pat$chd <- c(tmp_cf, tmp_scale)
      cf_t_pat$chd <- expr(
        (cf_t_all$chd$m2 %*% c("tch_int_hdl" = log(max(0.01,pat[["tch"]]/pat[["hdl"]])),
                               "db" = pat[["db"]])) %>% 
          as.vector() %>% exp()
      )
    } else if (tmp_male == 0 & tmp_chd == 0) {
      tmp_cf <- (cf_b_all$chd$f1[1:6] %*% c("int" = 1,
                                            "age" = tmp_age,
                                            "meno" = ifelse(tmp_age >= 50, 1, 0), # Assume meno at age 50
                                            "age_int_meno" = ifelse(tmp_age >= 50, tmp_age, 0),
                                            "smk" = tmp_smk,
                                            "alco" = tmp_alco)) %>% 
        as.vector() %>% exp()
      tmp_scale <- cf_b_all$chd$f1[[7]]
      cf_b_pat$chd <- c(tmp_cf, tmp_scale)
      cf_t_pat$chd <- expr(
        (cf_t_all$chd$f1 %*% c("tch_int_hdl" = log(max(0.01,pat[["tch"]]/pat[["hdl"]])),
                               "lnsbp" = log(max(0.01,pat[["sbp"]])),
                               "sbp_int_atht2" = ifelse(between(pat[["sbp"]], 110, 200) & pat[["atht"]] == 1, 
                                                        (200-pat[["sbp"]]) * (pat[["sbp"]] - 110) / 100, 0),
                               "db" = pat[["db"]],
                               "lntri" = log(max(0.01,pat[["tri"]])))) %>% 
          as.vector() %>% exp()
      )
    } else {
      tmp_cf <- (cf_b_all$chd$f2[1:3] %*% c("int" = 1,
                                            "age" = tmp_age,
                                            "smk" = tmp_smk)) %>% 
        as.vector() %>% exp()
      tmp_scale <- cf_b_all$chd$f2[[4]]
      cf_b_pat$chd <- c(tmp_cf, tmp_scale)
      cf_t_pat$chd <- expr(
        (cf_t_all$chd$f2 %*% c("tch_int_hdl" = log(max(0.01,pat[["tch"]]/pat[["hdl"]])),
                               "lnsbp" = log(max(0.01,pat[["sbp"]])),
                               "db" = pat[["db"]])) %>% 
          as.vector() %>% exp()
      )
    }
    
    # Stroke
    # - Cannot do it like death, because we have the expr needing global parameter
    if(tmp_male == 1){ 
      cf_b_pat$strk <- (cf_b_all$strk$m %*% c("int" = 1,
                                              "age_sub_50" = (tmp_age - 50) / 10,
                                              "lvh" = tmp_lvh,
                                              "smk" = tmp_smk,
                                              "af" = tmp_af)) %>% 
        as.vector() %>% exp()
      cf_t_pat$strk <- expr(
        (cf_t_all$strk$m %*% c("sbp_ctr" = (pat[["sbp"]] - 110) / 10, 
                               "sbp_int_atht2" = ifelse(between(pat[["sbp"]], 110, 200) & pat[["atht"]] == 1, 
                                                        (200-pat[["sbp"]]) * (pat[["sbp"]] - 110) / 100, 0),
                               "cvd" = pat[["chd"]] * pat[["strk"]],
                               "db" = pat[["db"]])) %>%
          as.vector() %>% exp()
      )
    } else {
      cf_b_pat$strk <- (cf_b_all$strk$f %*% c("int" = 1,
                                              "age_sub_50" = (tmp_age - 50) / 10,
                                              "lvh" = tmp_lvh,
                                              "smk" = tmp_smk,
                                              "af" = tmp_af)) %>% 
        as.vector() %>% exp()
      cf_t_pat$strk <- expr(
        (cf_t_all$strk$f %*% c("sbp_ctr" = (pat[["sbp"]] - 110) / 10, 
                               "sbp_int_atht2" = ifelse(between(pat[["sbp"]], 110, 200) & pat[["atht"]] == 1, 
                                                        (200-pat[["sbp"]]) * (pat[["sbp"]] - 110) / 100, 0),
                               "cvd" = pat[["chd"]] * pat[["strk"]],
                               "db" = pat[["db"]])) %>%
          as.vector() %>% exp()
      ) 
    }
    
    # Death
    if(tmp_male == 1){ tmp_cf <- "m" 
    } else { tmp_cf <- "f"} 
    cf_b_pat$death <- cf_b_all$death[[tmp_cf]] %>% 
      dplyr::filter(age >= tmp_age) %>% 
      dplyr::mutate(survprob = 1- prob, 
                    cumsurvprob = survprob * lag(survprob, default = 1)) %>% 
      dplyr::select(age, cumsurvprob) %>% 
      dplyr::bind_rows(tibble(age = 101, cumsurvprob=0)) %>% 
      dplyr::mutate(t = age - min(age)) %>% 
      dplyr::select(t, cumsurvprob)
    cf_b_pat$death_chd <- cf_b_all$death_chd[[tmp_cf]]
    cf_b_pat$death_strk <- cf_b_all$death_strk[[tmp_cf]]
    
    ### QOL ----
    
    # no impact on qol_all
    # qol_pat <- qol_all
    
    ### Cost ----
    
    cost_pat <- cost_all
    if(tmp_cvd == 1){ tmp_cf <- "cvd0" 
    } else { tmp_cf <- "cvd1"} 
    
    cost_pat$db <- cost_all$db[[tmp_cf]]
    cost_pat$chd <- cost_all$chd[[tmp_cf]]
    cost_pat$strk <- cost_all$strk[[tmp_cf]]
    cost_pat$death <- cost_all$death[[tmp_cf]]
    cost_pat$death_cv <- cost_all$death_cv[[tmp_cf]]
    
    ## Simulate first order ----
    
    outcome <- lapply(1:nsim, function(j) {
      
      # j <- 1
      set.seed(random_seed[i,j])
      
      # Initialize ----
      
      # > Attributes ----
      
      # Patient profiles (time varying one)
      pat <- pat_i[c("male", "atht", # Int term
                     "bmi","glu", "hdl", "tch","tri", "sbp",
                     "rcvr",
                     "td", "db", "chd", "strk",
                     "death", 
                     "line", "trt")]
      
      # Treatment profiles
      trt <- scn[[1]]
      cf_trt_pat <- cf_trts_pat[[trt]]
      
      # Outcomes
      out_qol <- out_all$qol
      out_cost <- out_all$cost
      out_other <- out_all$other
      
      # Supporting information 
      t_ltse <- c("db" = 0, "chd" = 0, "strk" = 0)
      anyltse <- c("db" = pat[["db"]],"chd" = 0, "strk" = 0)
      
      # > Time to event ----
      
      # Disease
      # - Start from stable
      evtar_rlps <- cf_trt_pat[["rlps"]]
      evttime_rlps <- rexp(1,cf_b_pat$rlps) * evtar_rlps
      evttime_rcvr <- Inf
      
      # Treatment sequence
      evttime_discon <- rexp(1, cf_trt_pat[["discon"]]) 
      
      # Side effect modelled directly
      evttime_td <- rexp(1, cf_trt_pat[["td"]]) 
      
      # Side effect modelled via metabolic impact
      
      # > DB
      if(pat["db"] == 0){
        evtar_db <- eval(cf_t_pat$db)
        evttime_db <- rexp(1, cf_b_pat$db) * evtar_db
      } else {
        evtar_db <- 1
        evttime_db <- Inf
      }
      
      # > CHD
      evtar_chd <- eval(cf_t_pat$chd)
      evttime_chd <- (-log(1-runif(1)))^ cf_b_pat$chd[[2]] * cf_b_pat$chd[[1]] * evtar_chd
      
      # AFT weibull 
      # - P = 1 - exp(-exp(u))
      # - u = (ln(t) - m)/ s
      # - m = cf %*% x
      # Reverse AFT weibull
      # - u = ln(-ln(1 - P))
      # - ln(t) = ln(-ln(1 - P)) * s + m
      # - t = exp(ln(-ln(1-P))* s) * exp(m)  
      #     = [exp(ln(-ln(1-P)))^s] * exp(m) 
      #     = [-ln(1-P)]^s * exp(m)
      #     = [-ln(1-P)]^Scale * exp(cf_b %*% x_b) * exp(cf_t %*% x_t)
      
      # > Stroke
      evtar_strk <- eval(cf_t_pat$strk)
      evttime_strk <- rexp(1, cf_b_pat$strk) * evtar_strk
      
      # > Death 
      tmp_ran <- runif(1)
      tmp_cf <- cf_b_pat$death$t[min(which(cf_b_pat$death$cumsurvprob <= tmp_ran))] 
      evttime_death <- tmp_cf + 0.5 
      
      # > Summary
      evtar <- c("rlps" = evtar_rlps,
                 "db" = evtar_db,
                 "chd" = evtar_chd,
                 "strk" = evtar_strk)
      evttime <- c("rlps" = evttime_rlps, 
                   "rcvr" = evttime_rcvr,
                   "discon" = evttime_discon,
                   "td" = evttime_td, 
                   "db" = evttime_db,
                   "chd" = evttime_chd,
                   "strk" = evttime_strk,
                   "death" = evttime_death) %>% sort()
      
      # Loop of reaction to event ----
      
      # Initialize simulation condition
      curtime <- 0
      reachtf <- 0
      
      while(curtime < tf_pat){
        
        ##  Update accrual ----
        
        ### Time ----
        
        # Get next time
        nexttime <- evttime[[1]]
        
        if(curtime + nexttime > tf_pat){ 
          nexttime <- tf_pat - curtime
          reachtf <- 1
        } 
        new_curtime <- curtime + nexttime
        
        if(pat[["rcvr"]] == 1) { out_other[["dur_rcvr"]] <- out_other[["dur_rcvr"]] + nexttime }
        
        ### QALY ----
        
        tmp_qol <- rep(0, 5)
        names(tmp_qol) <- c("dis", "stse", "ltse", "tot", "tot_disc")
        
        # Disease 
        tmp_qol[1] <- qol_all$cf[c("rcvr", "rlps")][[ifelse(pat[["rcvr"]] == 1, 1, 2)]]
        
        # TD
        if(pat[["td"]] == 1){ tmp_qol[2] <- qol_all$cf[["td"]] }
        
        # LTSE
        tmp_qol[3] <- as.vector(pat[c("db", "chd", "strk")] %*% qol_all$cf[c("db", "chd", "strk")])
        
        # TOTAL
        tmp_qol[4] <- sum(tmp_qol[1:3])
        
        # Discounted
        tmp_qol[5] <- f_disc_con(val = tmp_qol[4], disc = mod$disc_qaly, t1 = curtime, t2 = new_curtime)
        
        # Summary 
        # - tmp_qol * nexttime: QALY until next event
        for(k in c("dis", "stse", "ltse", "tot")){
          out_qol[[k]] <- out_qol[[k]] + tmp_qol[[k]] * nexttime
        }
        # - tmp_qol[5]: already QALY until next event
        out_qol[["tot_disc"]] <- out_qol[["tot_disc"]] + tmp_qol[5]
        
        ### Cost ----
        
        tmp_cost <- rep(0, 4)
        names(tmp_cost) <- c("drug", "ru", "tot", "tot_disc")
        
        # Drug cost 
        tmp_cost[1] <- pat[["trt"]] * cost_pat$trt[[trt]] # pat[["trt]] indicating whether currently on treatment or not
        
        # Resource use (Excluding rlps_acute: to be added at the event occurrence)
        tmp_cost[2] <- cost_pat$ru[[ifelse(pat[["rcvr"]] == 1, 1, 2)]]
        
        # TOTAL(Excluding LTSE: to be added at the end)
        tmp_cost[3] <- sum(tmp_cost[1:2])
        
        # Discounted
        tmp_cost[4] <- f_disc_con(val = tmp_cost[3], disc = mod$disc_cost, t1 = curtime, t2 = new_curtime) 
        
        # Summary
        for(k in c("drug", "ru", "tot")){
          out_cost[[k]] <- out_cost[[k]] + tmp_cost[[k]] * nexttime
        }
        out_cost[["tot_disc"]] <- out_cost[["tot_disc"]] + tmp_cost[4]
        
        ### calendar ----
        
        curtime <- new_curtime
        
        if(reachtf == 1){ break }
        
        ## React to next event ----
        
        ### General ----
        
        # Update time
        evttime <- evttime - nexttime
        
        # Take next event
        nextevt <- names(evttime)[1]
        
        ### Event specific ----
        
        # Create new attribute and event time easier to see the change
        next_pat <- pat
        next_evttime <- evttime
        next_evtar <- evtar
        
        if(nextevt == "rcvr"){ # Recovery 
          
          # Patient attribute
          next_pat[["rcvr"]] <- 1
          
          # Time to event
          
          # > Remove event
          next_evttime["rcvr"] <- Inf
          
          # > Add new event
          next_evttime[["rlps"]] <- rexp(1,cf_b_pat$rlps) * evtar[["rlps"]]
          next_evttime[["discon"]] <- rexp(1,cf_trt_pat[["discon"]])
          
        } else if (nextevt == "rlps"){ # Relapse 
          
          # Patient attribute
          # - start recovering
          # - no more relapse until recover 
          next_pat[["rcvr"]] <- 0
          
          # Time to event
          
          # > Remove event
          # - No further relapse during relapse
          # - No discontinuation during relapse
          next_evttime[c("discon", "rlps")] <- Inf
          
          # > New event
          # - Start time to recover
          next_evttime["rcvr"] <- cf_b_pat$rcvr
          
          # Outcome (treatment for acute phase relapse)
          
          # > Cost of acute relapse treatment
          tmp_cost <- cost_pat$ru[["rlps_acute"]]
            
          # > Accumulate costs
          out_cost[c("tot", "ru")] <- out_cost[c("tot", "ru")] + tmp_cost
          out_cost[["tot_disc"]] <- out_cost[["tot_disc"]] + f_disc_ins(val = tmp_cost, disc = mod$disc_cost, t = curtime)
          
          out_other[["n_rlps"]] <- out_other[["n_rlps"]] + 1
          
          # Subsequent potential consequence: switch
          if(pat[["line"]] < nline_all) {
            if(runif(1) <= cf_b_pat$switch){
              toswitch <- 1 
            } else { 
              toswitch <- 0
            }
          } else { 
            toswitch <- 0
          }
          
          if(toswitch == 1){
            
            # Patient attribute status
            new_line <- pat[["line"]] + 1
            new_trt <- scn[[new_line]]
            new_cf_trt_pat <- cf_trts_pat[[new_trt]]
            
            next_pat[["line"]] <- new_line
            next_pat[["trt"]] <- 1
            
            # > Reversible side effect
            tmp_qol <- 0
            tmp_cost <- 0
            
              # > = Agranulocytosis
            if(new_trt == "clo"){ 
              if(runif(1) <= new_cf_trt_pat[["agc"]]){ 
                if(runif(1) <= cf_b_pat$death_agc) {
                  next_pat[["death"]] <- 1
                  # Directly drop out with death (with profiles updated)
                  pat <- next_pat
                  break
                }
                tmp_base <- (qol_all$cf[c("rcvr", "td", "db", "chd", "strk")] %*% 
                  pat[c("rcvr", "td", "db", "chd", "strk")]) %>% as.vector() + 
                  ifelse(pat[["rcvr"]] == 0, qol_all$cf[["rlps"]], 0)
                tmp_qol <- tmp_qol + tmp_base * qol_all$agc
                tmp_cost <- tmp_cost + cost_pat$agc
              }
            }
              # > = Others: EPS, weight gain, sedation, sexual dysfunction
            tmp_qol <- tmp_qol + qol_all$stse[[new_trt]]
            tmp_cost <- tmp_cost + cost_pat$stse[[new_trt]]
            
              # > = Accumulate QALY and costs
            out_qol[c("tot", "stse")] <- out_qol[c("tot", "stse")] + tmp_qol
            out_cost[c("tot", "stse")] <- out_cost[c("tot", "stse")] + tmp_cost
            
              # > = Accumulate discounted QALY and costs
            out_qol[["tot_disc"]] <- out_qol[["tot_disc"]] + f_disc_ins(val = tmp_qol, disc = mod$disc_qaly, t = curtime)
            out_cost[["tot_disc"]] <- out_cost[["tot_disc"]] + f_disc_ins(val = tmp_cost, disc = mod$disc_cost, t = curtime)

            # > metabolic impact (which leading to impact on diabetes, CHD, stroke)
            
            # Old metabolic impact
            tmp_mbse <- c("bmi", "tch", "hdl", "tri", "sbp", "glu")
            if(pat[["trt"]] == 0){
              for(mbse in tmp_mbse){ next_pat[[mbse]] <- pat[[mbse]]  + new_cf_trt_pat[[mbse]] }
            } else {
              for(mbse in tmp_mbse){ next_pat[[mbse]] <- pat[[mbse]] - cf_trt_pat[[mbse]] + new_cf_trt_pat[[mbse]] }
            }
            
            # Time to event
            
            # > update event
            next_evtar[["rlps"]] <- new_cf_trt_pat[["rlps"]] # Currently there is no relapse to happen during relapse
            
            for(ltse in c("db", "chd", "strk")){
              if(anyltse[[ltse]] == 0){
                new_evtar <-  eval(cf_t_pat[[ltse]])
                next_evttime[[ltse]] <- evttime[[ltse]] / evtar[[ltse]] * new_evtar
                next_evtar[[ltse]] <- new_evtar
              }
            }
            
            # > New event
            if(pat[["td"]] == 0){ next_evttime[["td"]] <- rexp(1, new_cf_trt_pat[["td"]]) } 
            
            # Update treatment 
            trt <- new_trt
            cf_trt_pat <- new_cf_trt_pat
            
          } else if (pat[["trt"]] == 0){ # previously not on treatment: return to treatment at acute phase
            
            # Patient attribute
            
            # > Treatment
            next_pat[["trt"]] <- 1
            
            # > Metabolic 
            for(mbse in c("bmi", "tch", "hdl", "tri", "sbp", "glu")){  
              next_pat[[mbse]] <- pat[[mbse]] + cf_trt_pat[[mbse]] 
            }
            
            # > Reversible side effect: already tested previously
            
            # Time to event
            
            # > Update event
            next_evtar[["rlps"]] <- cf_trt_pat[["rlps"]] # Currently there is no relapse to happen during relapse
            
            for(ltse in c("db", "chd", "strk")){
              if(anyltse[[ltse]] == 0){
                new_evtar <-  eval(cf_t_pat[[ltse]])
                next_evttime[[ltse]] <- evttime[[ltse]] / evtar[[ltse]] * new_evtar
                next_evtar[[ltse]] <- new_evtar
              }
            }
            
            # > New event
            if(pat[["td"]] == 0){ next_evttime[["td"]] <- rexp(1, cf_trt_pat[["td"]]) } 
            
          } # else nothing will change
          
        } else if (nextevt == "discon"){  # Discontinuation 
          
          # Patient attribute
          
          # > Treatment status
          next_pat[["trt"]] <- 0
          
          # > Metabolic impact
          for(mbse in c("bmi", "tch", "hdl", "tri", "sbp", "glu")){
            next_pat[[mbse]] <- pat[[mbse]] - cf_trt_pat[[mbse]] 
          }
          
          # Time to event
          
          # > Remove event
          next_evttime[["discon"]] <- Inf
          if(pat[["td"]] == 0){ next_evttime[["td"]] <- Inf }
          
          # > Update event
          new_evtar <- 1
          next_evttime[["rlps"]] <- evttime[["rlps"]] / evtar[["rlps"]] * new_evtar
          next_evtar[["rlps"]] <- new_evtar
          
          for(ltse in c("db", "chd", "strk")){
            if(anyltse[[ltse]] == 0){
              new_evtar <-  eval(cf_t_pat[[ltse]])
              next_evttime[[ltse]] <- evttime[[ltse]] /evtar[[ltse]] * new_evtar
              next_evtar[[ltse]] <- new_evtar
            }
          }
          
        } else if (nextevt == "td"){ # Tardive dyskinesia
          
          # Patient attribute
          next_pat[["td"]] <- 1
          
          # Time to event
          next_evttime[["td"]] <- Inf
          
          # Outcome (treatment for acute phase relapse)
          
          # > Cost of management of TD (one time cost)
          tmp_cost <- cost_pat$td
          
          # > Accumulate costs
          out_cost[c("tot", "ru")] <- out_cost[c("tot", "ru")] + tmp_cost
          out_cost[["tot_disc"]] <- out_cost[["tot_disc"]] + f_disc_ins(val = tmp_cost, disc = mod$disc_cost, t = curtime)
          
        } else if (nextevt == "db"){ # Diabetes
          
          # Patient attribute
          next_pat[["db"]] <- 1
          anyltse[["db"]] <- 1
          t_ltse[["db"]] <- curtime
          # Time to event
          
          # > Remove event
          next_evttime[["db"]] <- Inf
          
          # > Update other long-term SE
          for(ltse in c("chd", "strk")){ 
            if(anyltse[[ltse]] == 0){
              new_evtar <-  eval(cf_t_pat[[ltse]])
              next_evttime[[ltse]] <- evttime[[ltse]] / evtar[[ltse]] * new_evtar
              next_evtar[[ltse]] <- new_evtar
            }
          }
          
        } else if (nextevt == "chd"){ # CHD
          
          # Patient attribute
          next_pat[["chd"]] <- 1
          anyltse[["chd"]] <- 1
          t_ltse[["chd"]] <- curtime
          
          # Update death
          if(runif(1) <= cf_b_pat$death_chd){ 
            next_pat[["death"]] <- 1
            pat <- next_pat
            break
          }
          
          # Time to event
          
          # > Remove event
          next_evttime[["chd"]] <- Inf
          
          # > Update event
          if(anyltse[["strk"]] == 0){
            new_evtar <-  eval(cf_t_pat[["strk"]])
            next_evttime[["strk"]] <- evttime[["strk"]] / evtar[["strk"]] * new_evtar
            next_evtar[["strk"]] <- new_evtar
          }      
          
        } else if (nextevt == "strk"){ # Stroke
          
          # Remove event
          next_pat[["strk"]] <- 1
          anyltse[["strk"]] <- 1
          t_ltse[["strk"]] <- curtime
          
          # Update death
          if(runif(1) <= cf_b_pat$death_strk){ 
            next_pat[["death"]] <- 1
            pat <- next_pat
            break
          }
          
          next_evttime[["strk"]] <- Inf
          
        } else if (nextevt == "death"){ # Death 
          
          # Patient attribute
          next_pat[["death"]] <- 1
          pat <- next_pat
          break
        }
        
        # Update everything
        
        # > Update order
        evttime <- sort(next_evttime)
        pat <- next_pat
        evtar <- next_evtar
        
      }
      
      # Update special ----
      
      # > Costs ----
      
      # CHD and Stroke
      for(ltse in c("chd", "strk")){
        
        t_first <- t_ltse[[ltse]]
        
        if(t_first > 0){
          
          tmp_dur <- curtime - t_first
          tmp_n <- ifelse(tmp_dur > 3, 3, floor(tmp_dur))
          
          tmp_cf <- cost_pat[[ltse]][1:(tmp_n + 1)]
          tmp_t <- map(0:tmp_n, ~c(.x, ifelse(.x == tmp_n, tmp_dur, .x + 1)) + t_first)
          tmp_t2 <- map(tmp_t, ~.x[2] - .x[1])
          
          tmp_cost <- map2_dbl(tmp_t2, tmp_cf, ~.x * .y) %>% sum()
          tmp_cost_disc <- map2_dbl(tmp_t, tmp_cf, 
                                    ~f_disc_con(val = .y, 
                                                disc = mod$disc_cost, 
                                                t1 = .x[1],
                                                t2 = .x[2])) %>% sum()
          out_cost[c("tot", "ltse")] <- out_cost[c("tot", "ltse")] + tmp_cost
          out_cost[["tot_disc"]] <- out_cost[["tot_disc"]] + tmp_cost_disc
          
        }
      }
      
      # Diabetes
      t_first <- t_ltse[["db"]]
      
      if(t_first >0){
        tmp_dur <- curtime - t_first
        if(tmp_dur < 10){
          tmp_cost <- tmp_dur * cost_pat$db[[1]]
          tmp_cost_disc <- f_disc_con(val = cost_pat$db[[1]], disc = mod$disc_cost, t1 = t_first, t2 = curtime)
        } else {
          tmp_cost <- 10 *  cost_pat$db[[1]] + (tmp_dur - 10) * cost_pat$db[[2]]
          tmp_cost_disc <- f_disc_con(val = cost_pat$db[[1]], disc = mod$disc_cost, t1 = t_first, t2 = t_first + 10) + 
            f_disc_con(val = cost_pat$db[[2]], disc = mod$disc_cost, t1 = t_first + 10, t2 = curtime)
        }
        out_cost[c("tot", "ltse")] <- out_cost[c("tot", "ltse")] + tmp_cost
        out_cost["tot_disc"] <- out_cost["tot_disc"] + tmp_cost_disc
      } 
      
      # Death
      if(pat[["death"]] == 1){
        if(nextevt %in% c("chd", "strk")){
          tmp_cost <- cost_pat$death_cv
        } else {
          tmp_cost <- cost_pat$death
        }
        out_cost[c("tot", "death")] <- out_cost[c("tot", "death")] + tmp_cost
        out_cost["tot_disc"] <- out_cost["tot_disc"] + 
          f_disc_ins(val = tmp_cost, disc = mod$disc_cost, t = curtime)
      }
      
      # > Incidence ----
      
      for(ltse in c("db", "chd", "strk")){
        if(t_ltse[[ltse]] != 0 & t_ltse[[ltse]] <10){
          out_other[[str_c(ltse, 10)]] <- 1
        }
      }

      # Output first order ----
      
      names(out_cost) <- str_c("cost_", names(out_cost))
      names(out_qol) <- str_c("qol_", names(out_qol))
      names(out_other) <- str_c("other_", names(out_other))
      results <-c(out_cost, out_qol, out_other, "dur_ly" = curtime)
      
      return(results)
      
    })
    
    # Output individual ----
    if(opt_ipd == TRUE){
      output_i <- bind_rows(outcome, .id = "id_sim")
    } else {
      output_i <- bind_rows(outcome) %>% summarize_all(~mean(.))
    }
    
    return(output_i)
  }
  
  # Output all ----
  
  stopCluster(cl)
  
  # return(bind_rows(output))
  if(opt_ipd == TRUE){
    output_fin <- bind_rows(output, .id = "id_pat")
  } else {
    output_fin <- bind_rows(output) %>% summarize_all(~mean(.)) # Will need to modify other_db10 as some people may have db at baseline
  }
  return(output_fin)
}

# Perform simulation ----
 
## Prepare input ----

dat_model <- read.csv(file.path(path_wd, "dat_model.csv"))
dat_lifetable <- read.csv(file.path(path_wd, "dat_lifetable.csv"))
input <- map(dat_model$value,~.x)
names(input) <- dat_model$parameter
input$lifetable <- dat_lifetable
para <- f_prepare_input(input)
saveRDS(para, file.path(path_wd, "mod_para.rds"))

## Example simulation ----

pat <- tibble(id = 1, rcvr = 1, line = 1, trt = 1, td = 0, death = 0,
              age = 25, bmi = 25, tch = 150, hdl = 40,
              tri = 140, glu = 80, sbp = 120,
              alco = 0, male = 1, smk = 0, par_db = 0, 
              atht = 0, af = 0, lvh = 0, db = 0, chd = 0, strk = 0,
              cmp_Full = 0, cmp_None = 0, cmp_Partial = 1) %>%
  as.matrix()
para <- readRDS(file.path(path_wd, "mod_para.rds"))
mod <- list(tf = expr(10), disc_cost = 0.035, disc_qaly = 0.035)
scn <- list(ola = c("ami", "ola", "clo"),
            ris2 = c("ami", "ris2", "clo"))
nsim <- 1000
ncore <- 1
rst <- map(scn, ~f_simulate_progress(pat, para, .x, mod, nsim, ncore))

tmp <- bind_rows(rst, .id = "scn") %>% 
  select(scn, cost = cost_tot_disc, qaly = qol_tot_disc) 
output <- bind_cols(tmp, tmp %>% filter(scn == "ola") %>% select(cost0 = cost, qaly0 = qaly)) %>% 
  mutate(inc_cost = cost - cost0, 
         inc_qaly = qaly - qaly0,
         icer = inc_cost / inc_qaly,
         nmb20k = inc_qaly * 20000 - inc_cost) %>%
  select(scn, cost, qaly, inc_cost, inc_qaly, icer, nmb20k)

write.csv(output, file = file.path(path_wd, "tbl_example.csv"))

## Base-case simulation ----

name_ap <- c("ami", "ari", "car", "lur", "ola", "que", "ris", "ari2", "ris2", "pal2")
name_pair_ap <- expand_grid(l1 = name_ap, l2 = name_ap) %>% filter(l1 != l2)
scns <- map2(name_pair_ap$l1,name_pair_ap$l2, ~c(.x, .y, "clo"))
names(scns) <- name_pair_ap %>% transmute(l = str_c(l1, "_",l2)) %>% pull(l)
mod <- list(tf = expr(10), disc_cost = 0.035, disc_qaly = 0.035)
para <- readRDS(file.path(path_wd, "mod_para.rds"))
# In the paper, we simulated for 1000 participants and each participant was simulated 100 times, but it took hours to run
# For code illustration, we only set the default number of participant to be 100 and number of simulation to be 1 (will take 6 minutes to run)
pats <- read.csv(file.path(path_wd, "dat_pats.csv"))[1:100,] %>% as.matrix()
nsim <- 1 
ncore <- 1
rst <- rep(list(NULL), n = length(scns))
t <- proc.time()
for(i in 1:length(scns)){ # 
  rst[[i]] <- f_simulate_progress(pats, para, scns[[i]], mod, nsim, ncore)
  print(i)
}
names(rst) <- names(scns)
tp <- bind_rows(rst, .id = "scn")
t2 <- proc.time()
t2 - t

output <- tp %>% 
  rename(cost = cost_tot_disc, 
         qaly = qol_tot_disc) %>% 
  select(scn, cost, qaly) %>% 
  mutate(l1 = str_sub(scn, 1, str_locate(scn, "_")[,1] - 1)) %>% 
  select(-scn) %>% 
  group_by(l1) %>% 
  summarize_all(~mean(.)) %>%
  rename(scn = l1) 

write.csv(output, file = file.path(path_wd, "tbl_basecase.csv"))
