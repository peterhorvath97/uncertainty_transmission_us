create_indicator2 <- function(data, variable, prob){
  
  require(dplyr)
  require(rlang)
  
  data <- data %>%
    mutate(q = quantile(!!sym(variable), probs = prob, na.rm = TRUE), 
           indicator = ifelse(!!sym(variable) > q, 1, 0)) %>% 
    select(date, indicator)
  
  data
}

create_indicator3 <- function(data, variable, upper, lower){
  
  require(dplyr)
  require(rlang)
  
  data <- data %>%
    mutate(q_up = quantile(!!sym(variable), probs = upper, na.rm = TRUE), 
           q_dn = quantile(!!sym(variable), probs = lower, na.rm = TRUE),
           indicator_up = ifelse(!!sym(variable) > q_up, 1, 0),
           indicator_dn = ifelse(!!sym(variable) < q_dn, 1, 0),
           indicator_mid = ifelse(indicator_up == 0 & indicator_dn == 0, 1, 0)) %>% 
    select(date, indicator_up, indicator_mid, indicator_dn)
  
  data
}

sign_tvar_2regime <- function(data, indicator, constr_mat, seed,
                              nlags, draws, subdraws, nkeep, KMIN, KMAX, constant, steps){
  
  require(dplyr)
  paste('US_only/codes/varsignr/', list.files('US_only/codes/varsignr') , sep = '') %>% lapply(source)
  
  
  indicator_inv <- indicator %>% mutate(indicator = ifelse(indicator == 1, 0, 1)) %>% as.matrix()
  indicator <- as.matrix(indicator)
  data_dn <- data[,1:(ncol(data))] %>% as.matrix()
  data_up <- data_dn
  
  
  for(i in 1:ncol(data_dn)){
    colnames(data_dn)[i] <- paste(colnames(data_dn)[i], "dn", sep = "_")
    colnames(data_up)[i] <- paste(colnames(data_up)[i], "up", sep = "_")
    data_dn[,i] <- data_dn[,i]*indicator_inv
    data_up[,i] <- data_up[,i]*indicator
  }
  
  data_ts <- bind_cols(data_up, data_dn) %>% ts()
  
  set.seed(seed)
  constr_up <- constr_mat[1,]
  model_up <- rwz.reject(Y=data_ts, nlags=nlags, draws=draws, subdraws=subdraws, nkeep=nkeep, KMIN=KMIN,
                         KMAX=6, constrained=constr_up, constant=constant, steps=steps)
  
  set.seed(seed)
  constr_dn <- constr_mat[2,]+ncol(constr_mat)*constr_mat[2,]/abs(constr_mat[2,])
  model_dn <- rwz.reject(Y=data_ts, nlags=nlags, draws=draws, subdraws=subdraws, nkeep=nkeep, KMIN=KMIN,
                         KMAX=KMAX, constrained=constr_dn, constant=constant, steps=steps)
  
  list(model_up, model_dn)
  
}

sign_tvar_2regime_controls <- function(data, indicator, controls = NA, constr_mat, seed,
                              nlags, draws, subdraws, nkeep, KMIN, KMAX, constant, steps){
  
  require(dplyr)
  paste('US_only/codes/varsignr/', list.files('US_only/codes/varsignr') , sep = '') %>% lapply(source)
  
  
  indicator_inv <- indicator %>% mutate(indicator = ifelse(indicator == 1, 0, 1)) %>% as.matrix()
  indicator <- as.matrix(indicator)
  data_dn <- data[,1:(ncol(data))] %>% as.matrix()
  data_up <- data_dn
  
  
  for(i in 1:ncol(data_dn)){
    colnames(data_dn)[i] <- paste(colnames(data_dn)[i], "dn", sep = "_")
    colnames(data_up)[i] <- paste(colnames(data_up)[i], "up", sep = "_")
    data_dn[,i] <- data_dn[,i]*indicator_inv
    data_up[,i] <- data_up[,i]*indicator
  }
  
  if (!is.logical(controls)){
    data_ts <- bind_cols(data_up, data_dn, controls) %>% ts()
  } else {
    data_ts <- bind_cols(data_up, data_dn) %>% ts()
    }
  
  
  
  set.seed(seed)
  constr_up <- constr_mat[1,]
  model_up <- rwz.reject(Y=data_ts, nlags=nlags, draws=draws, subdraws=subdraws, nkeep=nkeep, KMIN=KMIN,
                         KMAX=6, constrained=constr_up, constant=constant, steps=steps)
  
  set.seed(seed)
  constr_dn <- constr_mat[2,]+ncol(constr_mat)*constr_mat[2,]/abs(constr_mat[2,])
  model_dn <- rwz.reject(Y=data_ts, nlags=nlags, draws=draws, subdraws=subdraws, nkeep=nkeep, KMIN=KMIN,
                         KMAX=KMAX, constrained=constr_dn, constant=constant, steps=steps)
  
  list(model_up, model_dn)
  
}

sign_tvar_3regime <- function(data, indicator, constr_mat, seed,
                              nlags, draws, subdraws, nkeep, KMIN, KMAX, constant, steps){
  
  require(dplyr)
  paste('IMF paper/codes/varsignr/', list.files('IMF paper/codes/varsignr') , sep = '') %>% lapply(source)
  
  
  indicator_up <- indicator[,1] %>% as.matrix()
  indicator_mid <- indicator[,2] %>% as.matrix()
  indicator_dn <- indicator[,3] %>% as.matrix()
  
  indicator <- as.matrix(indicator)
  data_dn <- data[,1:(ncol(data))] %>% as.matrix()
  data_mid <- data_dn
  data_up <- data_dn
  
  
  for(i in 1:ncol(data_dn)){
    colnames(data_dn)[i] <- paste(colnames(data_dn)[i], "dn", sep = "_")
    colnames(data_mid)[i] <- paste(colnames(data_mid)[i], "mid", sep = "_")
    colnames(data_up)[i] <- paste(colnames(data_up)[i], "up", sep = "_")
    data_dn[,i] <- data_dn[,i]*indicator_dn
    data_mid[,i] <- data_mid[,i]*indicator_mid
    data_up[,i] <- data_up[,i]*indicator_up
  }
  
  data_ts <- bind_cols(data_up, data_mid, data_dn) %>% ts()
  set.seed(seed)
  constr_dn <- constr_mat[3,]+2*ncol(constr_mat)*constr_mat[3,]/abs(constr_mat[3,])
  model_dn <- rwz.reject(Y=data_ts, nlags=nlags, draws=draws, subdraws=subdraws, nkeep=nkeep, KMIN=KMIN,
                         KMAX=KMAX, constrained=constr_dn, constant=constant, steps=steps)
  
  
  set.seed(seed)
  constr_mid <- constr_mat[2,]+ncol(constr_mat)*constr_mat[2,]/abs(constr_mat[2,])
  model_mid <- rwz.reject(Y=data_ts, nlags=nlags, draws=draws, subdraws=subdraws, nkeep=nkeep, KMIN=KMIN,
                          KMAX=KMAX, constrained=constr_mid, constant=constant, steps=steps)
  
  set.seed(seed)
  constr_up <- constr_mat[1,]
  model_up <- rwz.reject(Y=data_ts, nlags=nlags, draws=draws, subdraws=subdraws, nkeep=nkeep, KMIN=KMIN,
                         KMAX=6, constrained=constr_up, constant=constant, steps=steps)
  
  list(model_up, model_mid, model_dn)
  
}

get_irf <- function(dim, prob, irflength, irfdata){
  
  require(dplyr)
  
  prob <- prob
  x <- irfdata[ , ,dim] %>% t() 
  t <- matrix(nrow = nrow(x), ncol = 1) 
  med <- matrix(nrow = nrow(x), ncol = 1)
  up <- matrix(nrow = nrow(x), ncol = 1)
  dn <- matrix(nrow = nrow(x), ncol = 1)
  for(i in(1:irflength)){
    t[i] <- i-1
    med[i] <- median(x[i,])
    up[i] <- quantile(x[i,], probs = (1-prob))
    dn[i] <- quantile(x[i,], probs = prob)
  }
  irf <- bind_cols(t, med, up, dn) %>%
    rename(t = ...1,
           median = ...2,
           upper = ...3,
           lower = ...4)
  irf
}

get_all_irfs <- function(mod_list, varname, regimename,
                         prob, irflength){
  irfdata_list <- NULL
  pivot_list <- NULL
  for(i in 1:length(mod_list)){
    for(j in 1:length(varname)){
      pivot_list[[j]] <- get_irf(dim = j+(i-1)*length(varname), prob = prob, irflength = irflength, irfdata = mod_list[[i]]$IRFS) %>% 
        mutate(label = varname[j],
               regime = regimename[i])
    }
    irfdata_list[[i]] <- pivot_list 
  }
  bind_rows(irfdata_list)  
}

get_all_irfs_clean <- function(mod_list, varname, regimename,
                               prob, irflength){
  get_all_irfs(mod_list, varname, regimename,
               prob, irflength) %>% 
    left_join(get_all_irfs(mod_list, varname, regimename,
                           prob, irflength) %>% 
                filter(t == 0, label == varname[1]) %>% 
                mutate(mult = 1/median) %>% 
                select(regime, mult), 
              by = "regime") %>% 
    mutate(label = factor(label, levels = varname),
           regime = factor(regime, levels = regimename)) %>% 
    mutate(median = median * mult,
           upper = upper * mult, 
           lower = lower * mult)
  
  
}

get_shock2 <- function(model, prob, indicatordata, modlags){
  t(apply(model[[1]]$SHOCKS, 2, quantile, probs = c(0.5, prob, 
                                                    1-prob))) %>% as_tibble() %>% 
    rename(median_1 = `50%`,
           lower_1 = `16%`,
           upper_1 = `84%`) %>% 
    bind_cols(t(apply(model[[2]]$SHOCKS, 2, quantile, probs = c(0.5, 0.16, 
                                                                0.84))) %>% as_tibble() %>% 
                rename(median_2 = `50%`,
                       lower_2 = `16%`,
                       upper_2 = `84%`)) %>% 
    bind_cols(indicatordata %>% 
                slice(-(1:modlags))) %>% 
    mutate(indicatorinv = ifelse(indicator == 1, 0, 1)) %>% 
    mutate(median = median_1 * indicator + median_2 * indicatorinv,
           lower = lower_1 * indicator + lower_2 * indicatorinv,
           upper = upper_1 * indicator + upper_2 * indicatorinv) %>% 
    select(date, median, lower, upper, indicator) %>% 
    mutate(regime = ifelse(indicator == 1, 'High uncertainty', 'Low uncertainty')) %>% 
    select(-indicator)
}