#loading packages
#packages
library(fredr)
library(tidyverse)
library(gridExtra)
library(VARsignR)
library(rlang)
library(lubridate)

#Custom functions
get_data_macro <- function(key){
  
  require(fredr)
  require(dplyr)
  require(purrr)
  
#set api key
fredr_set_key(key)

#download data
params <- list(
  series_id = c("DFF", "USACPIALLMINMEI", "INDPRO"),
  frequency = "m",
  observation_start = as.Date("1950-01-01")
)


data_macro  <- pmap_dfr(
  .l = params,
  .f = ~ fredr(series_id = .x, frequency = .y)
) %>%
  dplyr::select(date, series_id, value) %>%
  spread(key = series_id, value = value) %>%
  drop_na() %>% rename(ffr = DFF,
                       cpi = USACPIALLMINMEI,
                       indpro = INDPRO
  ) 

data_macro
}

get_data_unc <- function(key){
  
  require(fredr)
  require(dplyr)
  require(purrr)
  
  #set api key
  fredr_set_key(key)
  
  #download data
  params <- list(
    series_id = c("USEPUINDXD", "TEDRATE", "VIXCLS"),
    frequency = "m",
    observation_start = as.Date("1950-01-01")
  )
  
  
  data_unc  <- pmap_dfr(
    .l = params,
    .f = ~ fredr(series_id = .x, frequency = .y)
  ) %>%
    select(date, series_id, value) %>%
    spread(key = series_id, value = value) %>%
    rename(epu = USEPUINDXD,
           ted = TEDRATE,
           vix = VIXCLS
    ) 
  
  data_unc
}

cortab_create <- function(data, digits){
  cortab <- data %>% select(-date) %>% drop_na() %>% cor()
  cortab <- round(cortab, digits = digits)
  cortab <- cortab*lower.tri(cortab)
  diag(cortab) <- 1
  for(i in 1:nrow(cortab)) {
    for(j in 1:ncol(cortab)){
      if(cortab[i,j] == 0) {
        cortab[i,j] <- " "
      }
    }
  }
  cortab
  rownames(cortab) <- c("MPU", "TED", "EPU", "VIX")
  colnames(cortab) <- c("MPU", "TED", "EPU", "VIX")
  cortab
}

create_indicator2 <- function(data, variable, prob){
  
  require(dplyr)
  require(rlang)
  
  data <- data %>%
    mutate(q = quantile(!!sym(variable), probs = prob), 
           indicator = ifelse(!!sym(variable) > q, 1, 0)) %>% 
    select(date, indicator)
  
  data
}

create_indicator3 <- function(data, variable, upper, lower){
  
  require(dplyr)
  require(rlang)
  
  data <- data %>%
    mutate(q_up = quantile(!!sym(variable), probs = upper), 
           q_dn = quantile(!!sym(variable), probs = lower),
           indicator_up = ifelse(!!sym(variable) > q_up, 1, 0),
           indicator_dn = ifelse(!!sym(variable) < q_dn, 1, 0),
           indicator_mid = ifelse(indicator_up == 0 & indicator_dn == 0, 1, 0)) %>% 
    select(date, indicator_up, indicator_mid, indicator_dn)
  
  data
}

sign_tvar_2regime <- function(data, indicator, constr_mat, seed,
                              nlags, draws, subdraws, nkeep, KMIN, KMAX, constant, steps){
  
  require(dplyr)
  require(VARsignR)
  
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

sign_tvar_3regime <- function(data, indicator, constr_mat, seed,
                              nlags, draws, subdraws, nkeep, KMIN, KMAX, constant, steps){
  
  require(dplyr)
  require(VARsignR)
  
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

#Hívjunk meg adatokat
data_unc <- get_data_unc("cda47ae66b38ed7988c0a9c2ec80c94f")
data_macro <- get_data_macro("cda47ae66b38ed7988c0a9c2ec80c94f")

data_macro <- data_macro %>% 
  inner_join(seq(from = max(min(data_macro$date), min(data_unc$date)), 
                 to = min(max(data_macro$date), max(data_unc$date), as.Date("2019-12-01")),
                 by = 1) %>% 
               as_tibble() %>% 
               mutate(day = day(value)) %>% filter(day == 1) %>% select(-day) %>% 
               rename(date = value),
             by = "date")

data_unc <- data_unc %>% 
  inner_join(seq(from = max(min(data_macro$date), min(data_unc$date)), 
                 to = min(max(data_macro$date), max(data_unc$date), as.Date("2019-12-01")),
                 by = 1) %>% 
               as_tibble() %>% 
               mutate(day = day(value)) %>% filter(day == 1) %>% select(-day) %>% 
               rename(date = value),
             by = "date")

data_macro <- data_macro %>% filter(year(date) <= 2007)
data_unc <- data_unc %>% filter(year(date) <= 2007)

#Rajzoljuk ki, hogymivel dolgozunk
#itt előre le kéne filterelni
data_macro %>%
  gather(key = "variable", value = "value", ffr, cpi, indpro) %>%
  mutate(variable = case_when(variable == "ffr" ~ "Federal Funds Rate",
                              variable == "cpi" ~ "Consumer Prices",
                              variable == "indpro" ~ "Industrial Production Index"),
         variable = factor(variable, levels = c("Federal Funds Rate", "Consumer Prices", "Industrial Production Index"))) %>%
  ggplot(aes(x = date, y = value)) +
  geom_line() +
  facet_wrap(~variable, scales = "free", nrow = 3) + theme_minimal() +
  labs(x = "",
       y = "") 

data_unc %>% 
  gather(key = "variable", value = "value", mpu, epu, ted, vix) %>% 
  mutate(variable = case_when(variable == "epu" ~ "EPU index",
                               variable == "ted" ~ "TED spread",
                               variable == "vix" ~ "VIX index"),
         variable = factor(variable, levels = c("EPU index", "TED spread", "VIX index"))) %>% 
  ggplot(aes(x = date, y = value, group = variable, color = factor(variable))) +
  geom_line() + 
  geom_smooth(method = "gam", se = FALSE, linetype = "dashed") +
  scale_y_log10() +
  scale_color_viridis_d(begin = 0.05, end = 0.4, name = "Index", option = "turbo") +
  theme_minimal() +
  labs(x = "",
       y = "") 

data_unc %>% 
  select(date, mpu) %>% 
  mutate(mpu_med = median(mpu),
         mpu_25 = quantile(mpu, 0.25),
         mpu_75 = quantile(mpu, 0.75)) %>% 
  gather(key = "variable", value = "value", mpu, mpu_med, mpu_25, mpu_75) %>% 
  mutate(variable = case_when(variable == "epu" ~ "EPU",
                              variable == "epu_med" ~ "Median",
                              variable == "epu_25" ~ "Lower quartile",
                              variable == "epu_75" ~ "Upper quartile")) %>% 
  ggplot(aes(x = date, y = value, group = variable, size = variable, color = variable)) +
  geom_line() +
  scale_size_manual(breaks = c("EPU", "Median", "Lower quartile", "Upper quartile"),
                    values = c(.5, 1, 1, 1)) +
  scale_color_manual(breaks = c("EPU", "Median", "Lower quartile", "Upper quartile"),
                     values = c("black", "#df2c14", "#0457ac", "#0457ac")) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = "",
       y = "")

#Korrelációkat ki kéne számolni az indexek között
cortab <- cortab_create(data_unc, 3)

#Vegyünk első differenciákat
data_macro <- data_macro %>%
  mutate(cpi = cpi - lag(cpi),
         indpro = indpro - lag(indpro)) %>% drop_na()

data_unc <- data_unc[2:nrow(data_unc), ]



#Rajzoljuk ki regime indicator szerint színezve
#nemakaromlátni
#####
p2 <- data_macro %>% 
  inner_join(create_indicator2(data_unc, "epu", 0.5), by = "date") %>% 
  mutate(indicator = ifelse(indicator == 1, "High uncertainty regime", "Low uncertainty regime")) %>% 
  gather(key = "variable", value = "value", ffr, cpi, indpro) %>%
  mutate(variable = case_when(variable == "ffr" ~ "Federal Funds Rate",
                              variable == "cpi" ~ "Consumer Prices",
                              variable == "indpro" ~ "Industrial Production Index"),
         variable = factor(variable, levels = c("Federal Funds Rate", "Consumer Prices", "Industrial Production Index"))) %>%
  ggplot(aes(x = date, y = value, color = indicator, group = 1)) +
  geom_line(size = .8) +
  facet_wrap(~variable, scales = "free", nrow = 3) + theme_minimal() +
  scale_color_manual(breaks = c("High uncertainty regime", "Low uncertainty regime"),
                     values = c("#E9002D", "#00B000")) +
  labs(x = "",
       y = "",
       caption = "Shaded areas indicate NBER recession preiods") + 
  theme(legend.title = element_blank()) +
  geom_rect(data = fredr(series_id = "USREC",
                         frequency = "m") %>% 
              select(date, recession = value) %>% 
              mutate(diff = recession - lag(recession)) %>% 
              filter(!is.na(diff)) %>% 
              mutate(recession_start = ifelse(diff == 1, as.character(date), NA),
                     recession_end = ifelse(diff == -1, as.character(date), NA)) %>% 
              filter(!is.na(recession_start) | !is.na(recession_end)) %>% 
              mutate(recession_end = ifelse(!is.na(recession_start), lead(recession_end), NA)) %>%
              filter(!is.na(recession_start)) %>% 
              mutate(across(.cols = c(recession_start, recession_end), .fns = ~as.Date(.x))) %>% 
              select(recession_start, recession_end) %>% 
              filter(recession_start >= min(data_unc$date),
                     recession_start <= max(data_unc$date)),
            inherit.aes = F,
            aes(xmin = recession_start, xmax = recession_end, ymin = -Inf, ymax = Inf), 
                fill = "grey50", alpha = 0.5)
  
  
p3 <- data_macro %>% 
  inner_join(create_indicator3(data_unc, "epu", 0.75, 0.25), by = "date") %>% 
  mutate(indicator = case_when(indicator_up == 1 ~ "High uncertainty regime",
                               indicator_mid == 1 ~ "Medium uncertainty regime",
                               indicator_dn == 1 ~ "Low uncertainty regime"),
         indicator = factor(indicator, levels = c("High uncertainty regime",
                                                  "Medium uncertainty regime",
                                                  "Low uncertainty regime"))) %>% 
  gather(key = "variable", value = "value", ffr, cpi, indpro) %>%
  mutate(variable = case_when(variable == "ffr" ~ "Federal Funds Rate",
                              variable == "cpi" ~ "Consumer Prices",
                              variable == "indpro" ~ "Industrial Production Index"),
         variable = factor(variable, levels = c("Federal Funds Rate", "Consumer Prices", "Industrial Production Index"))) %>%
  ggplot(aes(x = date, y = value, color = indicator, group = 1)) +
  geom_line(size = .8) +
  facet_wrap(~variable, scales = "free", nrow = 3) + theme_minimal() +
  scale_color_manual(breaks = c("High uncertainty regime", "Medium uncertainty regime", "Low uncertainty regime"),
                     values = c("#E9002D", "#9EC9E2" , "#00B000")) +
  labs(x = "",
       y = "",
       caption = "Shaded areas indicate NBER recession preiods") + 
  theme(legend.title = element_blank()) +
  geom_rect(data = fredr(series_id = "USREC",
                         frequency = "m") %>% 
              select(date, recession = value) %>% 
              mutate(diff = recession - lag(recession)) %>% 
              filter(!is.na(diff)) %>% 
              mutate(recession_start = ifelse(diff == 1, as.character(date), NA),
                     recession_end = ifelse(diff == -1, as.character(date), NA)) %>% 
              filter(!is.na(recession_start) | !is.na(recession_end)) %>% 
              mutate(recession_end = ifelse(!is.na(recession_start), lead(recession_end), NA)) %>%
              filter(!is.na(recession_start)) %>% 
              mutate(across(.cols = c(recession_start, recession_end), .fns = ~as.Date(.x))) %>% 
              select(recession_start, recession_end) %>% 
              filter(recession_start >= min(data_unc$date),
                     recession_start <= max(data_unc$date)),
            inherit.aes = F,
            aes(xmin = recession_start, xmax = recession_end, ymin = -Inf, ymax = Inf), 
            fill = "grey50", alpha = 0.5)
#####
grid.arrange(p2, p3)


#2Regime model
mod_2regime <- sign_tvar_2regime(data = data_macro %>% 
                               inner_join(create_indicator2(data_unc, "epu", 0.5), by = "date") %>% 
                               select(-date, -indicator), 
                  indicator = create_indicator2(data_unc, "epu", 0.5) %>% 
                    inner_join(data_macro, by = "date") %>% select(indicator), 
                  constr_mat = rbind(c(+1, -2, -3) ,
                                     c(+1, -2, -3)), seed = 2022,
                  nlags = 3, draws = 20000, subdraws = 200, nkeep = 1000, 
                  KMIN = 1, KMAX = 6, constant = TRUE, steps = 21)




#3Regime model
mod_3regime <- sign_tvar_3regime(data = data_macro %>% 
                               inner_join(create_indicator3(data_unc, "epu", 0.75, 0.25), by = "date") %>% 
                               select(-date, -starts_with("indicator")), 
                             indicator = create_indicator3(data_unc, "epu", 0.75, 0.25) %>% 
                               inner_join(data_macro, by = "date") %>% select(starts_with("indicator")), 
                             constr_mat = rbind(c(+1, -2, -3) ,
                                                c(+1, -2, -3),
                                                c(+1, -2, -3)), seed = 2022,
                             nlags = 3, draws = 20000, subdraws = 200, nkeep = 1000, 
                             KMIN = 1, KMAX = 6, constant = TRUE, steps = 21)




get_all_irfs(mod_list = mod_2regime, 
             varname = c("FFR", "CPI", "INDPRO"),
             regimename = c("High uncertainty", "Low uncertainty"),
             prob = 0.16,
             irflength = 21) %>% 
  left_join(get_all_irfs(mod_list = mod_2regime, 
                         varname = c("FFR", "CPI", "INDPRO"),
                         regimename = c("High uncertainty", "Low uncertainty"),
                         prob = 0.16,
                         irflength = 21) %>%
              filter(t == 0, label == "FFR") %>% 
              mutate(mult = 1/median) %>% 
              select(regime, mult), 
            by = "regime") %>% 
  mutate(label = factor(label, levels = c("FFR", "CPI", "INDPRO")),
         regime = factor(regime, levels = c("High uncertainty", "Low uncertainty"))) %>% 
  mutate(median = median * mult,
         upper = upper * mult, 
         lower = lower * mult) %>% 
  ggplot(aes(x = t, y = median, ymin = lower, ymax = upper)) +
  geom_line() +
  geom_hline(yintercept = 0, color = "red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") + 
  facet_wrap(regime~label, scales = "free")+
  theme_minimal() +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11)) +
  labs(x = "",
       y ="")


get_all_irfs(mod_list = mod_3regime, 
                          varname = c("FFR", "CPI", "INDPRO"),
                          regimename = c("High uncertainty", "Medium uncertainty", "Low uncertainty"),
                          prob = 0.16,
                          irflength = 21) %>% 
  left_join(get_all_irfs(mod_list = mod_3regime, 
                         varname = c("FFR", "CPI", "INDPRO"),
                         regimename = c("High uncertainty", "Medium uncertainty", "Low uncertainty"),
                         prob = 0.16,
                         irflength = 21) %>%
              filter(t == 0, label == "FFR") %>% 
     mutate(mult = 1/median) %>% 
              select(regime, mult), 
            by = "regime") %>% 
  mutate(label = factor(label, levels = c("FFR", "CPI", "INDPRO")),
         regime = factor(regime, levels = c("High uncertainty", "Medium uncertainty", "Low uncertainty"))) %>% 
  mutate(median = median * mult,
         upper = upper * mult, 
         lower = lower * mult) %>% 
  ggplot(aes(x = t, y = median, ymin = lower, ymax = upper)) +
  geom_line() +
  geom_hline(yintercept = 0, color = "red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") + 
  facet_wrap(regime~label, scales = "free")+
  theme_minimal() +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11)) +
  labs(x = "",
       y ="")

#retrieve coefficients
coefs <- function(model, regime){
coef <- model[[regime]]$BDraws
coef_mat <- coef[1,,] %>% as_tibble() %>% 
  bind_cols(coef[1,,] %>% rownames())
for(i in 2:807){

  coef_mat <- bind_rows(coef_mat, coef[i,,] %>% as_tibble() %>% 
                          bind_cols(coef[i,,] %>% rownames))

}

coef_mat
}

coefmat <- NULL

for(i in 1:3){
  coefmat[[i]] <- coefs(mod_3regime, i)
}

bind_rows(coefmat) %>% 
  rename(coef = ...10) %>% 
  group_by(coef) %>% 
  summarize(across(everything(), ~median(.x, na.rm = TRUE))) %>% 
  select(coef, contains("ffr")) %>% 
  print(n = 100)

indicator <- create_indicator3(data_unc, "epu", 0.75, 0.25)
indicator_up <- indicator[,2] %>% as.matrix()
indicator_mid <- indicator[,3] %>% as.matrix()
indicator_dn <- indicator[,4] %>% as.matrix()

indicator <- as.matrix(indicator %>% select(-date))
data <- data_macro %>% 
  inner_join(create_indicator3(data_unc, "epu", 0.75, 0.25), by = "date") %>% 
  select(-date, -starts_with("indicator"))
data_dn <- data[,1:(ncol(data))] %>% as.matrix()
data_mid <- data_dn
data_up <- data_dn


for(i in 1:(ncol(data_dn))){
  colnames(data_dn)[i] <- paste(colnames(data_dn)[i], "dn", sep = "_")
  colnames(data_mid)[i] <- paste(colnames(data_mid)[i], "mid", sep = "_")
  colnames(data_up)[i] <- paste(colnames(data_up)[i], "up", sep = "_")
  data_dn[,i] <- data_dn[,i]*indicator_dn
  data_mid[,i] <- data_mid[,i]*indicator_mid
  data_up[,i] <- data_up[,i]*indicator_up
}

data_ts <- bind_cols(data_up, data_mid, data_dn) 

bind_cols(data_macro$date, data_ts) %>% 
  rename(date = ...1) %>% 
  mutate(ffr_up_est = 0.0246 
         + -0.419*lag(cpi_dn) + 0.865*lag(cpi_dn, n = 2) + -0.598*lag(cpi_dn, n = 3)
         +-0.592*lag(cpi_mid) + 0.0805*lag(cpi_mid, n = 2) + -0.0369*lag(cpi_mid, n = 3) 
         +0.252*lag(cpi_up) + 0.108*lag(cpi_up, n = 2) + 0.112*lag(cpi_up, n = 3)
         + -0.438*lag(ffr_dn) + 0.681*lag(ffr_dn, n = 2) + -0.211*lag(ffr_dn, n = 3)
         +-0.364*lag(ffr_mid) + 0.625*lag(ffr_mid, n = 2) + -0.147*lag(ffr_mid, n = 3) 
         +0.184*lag(ffr_up) + 0.584*lag(ffr_up, n = 2) + -0.0763*lag(ffr_up, n = 3)
         + -0.122*lag(indpro_dn) + -0.147*lag(indpro_dn, n = 2) + -0.158*lag(indpro_dn, n = 3)
         +0.155*lag(indpro_mid) + -0.261*lag(indpro_mid, n = 2) + -0.0638*lag(indpro_mid, n = 3) 
         +0.000391*lag(indpro_up) + 0.0138*lag(indpro_up, n = 2) + -0.0994*lag(indpro_up, n = 3),
         ffr_mid_est = -0.125 
         + -0.419*lag(cpi_dn) + -0.297*lag(cpi_dn, n = 2) + 0.293*lag(cpi_dn, n = 3)
         +0.303*lag(cpi_mid) + 0.274*lag(cpi_mid, n = 2) + 0.590*lag(cpi_mid, n = 3) 
         +-0.260*lag(cpi_up) + 0.0159*lag(cpi_up, n = 2) + -0.385*lag(cpi_up, n = 3)
         + 1.62*lag(ffr_dn) + -1.92*lag(ffr_dn, n = 2) + 0.708*lag(ffr_dn, n = 3)
         +1.9*lag(ffr_mid) + -1.85*lag(ffr_mid, n = 2) + 0.669*lag(ffr_mid, n = 3) 
         +1.5*lag(ffr_up) + -1.8*lag(ffr_up, n = 2) + 0.664*lag(ffr_up, n = 3)
         + -0.420*lag(indpro_dn) + -0.305*lag(indpro_dn, n = 2) + -0.661*lag(indpro_dn, n = 3)
         +-0.119*lag(indpro_mid) + 0.260*lag(indpro_mid, n = 2) + 0.105*lag(indpro_mid, n = 3) 
         +0.120*lag(indpro_up) + -0.0275*lag(indpro_up, n = 2) + 0.132*lag(indpro_up, n = 3),
         ffr_dn_est = 0.0974 
         + 0.959*lag(cpi_dn) +  -0.602*lag(cpi_dn, n = 2) + 0.320*lag(cpi_dn, n = 3)
         +0.228*lag(cpi_mid) + -0.303*lag(cpi_mid, n = 2) + -0.605*lag(cpi_mid, n = 3) 
         +-0.0639*lag(cpi_up) + -0.110*lag(cpi_up, n = 2) + 0.205*lag(cpi_up, n = 3)
         + 0.197*lag(ffr_dn) + 0.936*lag(ffr_dn, n = 2) + -0.583*lag(ffr_dn, n = 3)
         +-0.181*lag(ffr_mid) + 0.950*lag(ffr_mid, n = 2) + 0.669*lag(ffr_mid, n = 3) 
         +-0.325*lag(ffr_up) + 0.913*lag(ffr_up, n = 2) + -0.633 *lag(ffr_up, n = 3)
         +  0.581 *lag(indpro_dn) + 0.447*lag(indpro_dn, n = 2) + 0.832*lag(indpro_dn, n = 3)
         +0.00263*lag(indpro_mid) + 0.0263*lag(indpro_mid, n = 2) +  0.00669*lag(indpro_mid, n = 3) 
         +-0.0325 *lag(indpro_up) + 0.00431*lag(indpro_up, n = 2) + -0.0375*lag(indpro_up, n = 3)) %>% 
  bind_cols(indicator) %>% 
  mutate(ffr_est = indicator_up * ffr_up_est + indicator_mid * ffr_mid_est + indicator_dn * ffr_dn_est,
         ffr = ffr_up + ffr_mid + ffr_dn) %>% 
  select(date, ffr, ffr_est) %>% 
  mutate(innov = ffr-ffr_est) %>% 
  ggplot(aes(x = date, y = innov)) +
  geom_line(size = .6) +
  theme_minimal() +
  labs(y = "",
       x = "") 

bind_cols(data_macro$date, data_ts) %>% 
  rename(date = ...1) %>% 
  mutate(ffr_up_est = 0.0246 
         + -0.419*lag(cpi_dn) + 0.865*lag(cpi_dn, n = 2) + -0.598*lag(cpi_dn, n = 3)
         +-0.592*lag(cpi_mid) + 0.0805*lag(cpi_mid, n = 2) + -0.0369*lag(cpi_mid, n = 3) 
         +0.252*lag(cpi_up) + 0.108*lag(cpi_up, n = 2) + 0.112*lag(cpi_up, n = 3)
         + -0.438*lag(ffr_dn) + 0.681*lag(ffr_dn, n = 2) + -0.211*lag(ffr_dn, n = 3)
         +-0.364*lag(ffr_mid) + 0.625*lag(ffr_mid, n = 2) + -0.147*lag(ffr_mid, n = 3) 
         +0.184*lag(ffr_up) + 0.584*lag(ffr_up, n = 2) + -0.0763*lag(ffr_up, n = 3)
         + -0.122*lag(indpro_dn) + -0.147*lag(indpro_dn, n = 2) + -0.158*lag(indpro_dn, n = 3)
         +0.155*lag(indpro_mid) + -0.261*lag(indpro_mid, n = 2) + -0.0638*lag(indpro_mid, n = 3) 
         +0.000391*lag(indpro_up) + 0.0138*lag(indpro_up, n = 2) + -0.0994*lag(indpro_up, n = 3),
         ffr_mid_est = -0.125 
         + -0.419*lag(cpi_dn) + -0.297*lag(cpi_dn, n = 2) + 0.293*lag(cpi_dn, n = 3)
         +0.303*lag(cpi_mid) + 0.274*lag(cpi_mid, n = 2) + 0.590*lag(cpi_mid, n = 3) 
         +-0.260*lag(cpi_up) + 0.0159*lag(cpi_up, n = 2) + -0.385*lag(cpi_up, n = 3)
         + 1.62*lag(ffr_dn) + -1.92*lag(ffr_dn, n = 2) + 0.708*lag(ffr_dn, n = 3)
         +1.9*lag(ffr_mid) + -1.85*lag(ffr_mid, n = 2) + 0.669*lag(ffr_mid, n = 3) 
         +1.5*lag(ffr_up) + -1.8*lag(ffr_up, n = 2) + 0.664*lag(ffr_up, n = 3)
         + -0.420*lag(indpro_dn) + -0.305*lag(indpro_dn, n = 2) + -0.661*lag(indpro_dn, n = 3)
         +-0.119*lag(indpro_mid) + 0.260*lag(indpro_mid, n = 2) + 0.105*lag(indpro_mid, n = 3) 
         +0.120*lag(indpro_up) + -0.0275*lag(indpro_up, n = 2) + 0.132*lag(indpro_up, n = 3),
         ffr_dn_est = 0.0974 
         + 0.959*lag(cpi_dn) +  -0.602*lag(cpi_dn, n = 2) + 0.320*lag(cpi_dn, n = 3)
         +0.228*lag(cpi_mid) + -0.303*lag(cpi_mid, n = 2) + -0.605*lag(cpi_mid, n = 3) 
         +-0.0639*lag(cpi_up) + -0.110*lag(cpi_up, n = 2) + 0.205*lag(cpi_up, n = 3)
         + 0.197*lag(ffr_dn) + 0.936*lag(ffr_dn, n = 2) + -0.583*lag(ffr_dn, n = 3)
         +-0.181*lag(ffr_mid) + 0.950*lag(ffr_mid, n = 2) + 0.669*lag(ffr_mid, n = 3) 
         +-0.325*lag(ffr_up) + 0.913*lag(ffr_up, n = 2) + -0.633 *lag(ffr_up, n = 3)
         +  0.581 *lag(indpro_dn) + 0.447*lag(indpro_dn, n = 2) + 0.832*lag(indpro_dn, n = 3)
         +0.00263*lag(indpro_mid) + 0.0263*lag(indpro_mid, n = 2) +  0.00669*lag(indpro_mid, n = 3) 
         +-0.0325 *lag(indpro_up) + 0.00431*lag(indpro_up, n = 2) + -0.0375*lag(indpro_up, n = 3)) %>% 
  bind_cols(indicator) %>% 
  mutate(ffr_est = indicator_up * ffr_up_est + indicator_mid * ffr_mid_est + indicator_dn * ffr_dn_est,
         ffr = ffr_up + ffr_mid + ffr_dn) %>% 
  select(date, ffr, ffr_est) %>% 
  mutate(innov = ffr-ffr_est) %>% 
  gather(key = "key", value = "value", ffr, ffr_est) %>% 
  ggplot(aes(x = date, y = value, color = key)) +
  geom_line(size = .6, )+
  theme_minimal() +
  labs(y = "",
       x = "") 
  

#Innentől robustness check
#Rob check a lag-okra
list_lags <- NULL
for(i in 1:3){
  list_lags[[i]] <- 
    sign_tvar_3regime(data = data_macro %>% 
                        inner_join(create_indicator3(data_unc, "epu", 0.75, 0.25), by = "date") %>% 
                        select(-date, -starts_with("indicator")), 
                      indicator = create_indicator3(data_unc, "epu", 0.75, 0.25) %>% 
                        inner_join(data_macro, by = "date") %>% select(starts_with("indicator")), 
                      constr_mat = rbind(c(+1, -2, -3) ,
                                         c(+1, -2, -3),
                                         c(+1, -2, -3)), seed = 2022,
                      nlags = i, draws = 20000, subdraws = 200, nkeep = 1000, 
                      KMIN = 1, KMAX = 6, constant = TRUE, steps = 21)

}
irf_list <- NULL
for(i in 1:3){
  irf_list[[i]] <-
    get_all_irfs(mod_list = list_lags[[i]], 
                 varname = c("FFR", "CPI", "INDPRO"),
                 regimename = c("High uncertainty", "Medium uncertainty", "Low uncertainty"),
                 prob = 0.16,
                 irflength = 21) %>% 
    left_join(get_all_irfs(mod_list = list_lags[[i]], 
                           varname = c("FFR", "CPI", "INDPRO"),
                           regimename = c("High uncertainty", "Medium uncertainty", "Low uncertainty"),
                           prob = 0.16,
                           irflength = 21) %>%
                filter(t == 0, label == "FFR") %>% 
                mutate(mult = 1/median) %>% 
                select(regime, mult), 
              by = "regime") %>% 
    mutate(label = factor(label, levels = c("FFR", "CPI", "INDPRO")),
           regime = factor(regime, levels = c("High uncertainty", "Medium uncertainty", "Low uncertainty"))) %>%
    mutate(median = median * mult,
           upper = upper * mult, 
           lower = lower * mult) %>% 
    mutate(lag = paste("Lags = ", i, sep = ""))
}

bind_rows(irf_list) %>% 
  ggplot(aes(x = t, y = median, color = lag, group = lag)) +
  geom_line( size = 1) +
  geom_hline(yintercept = 0, color = "red") +
  facet_wrap(regime~label, scales = "free")+
  scale_color_viridis_d(begin = 0.05, end = 0.8, name = "Index", option = "magma")+
  theme_minimal() +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11)) +
  labs(x = "",
       y ="") +
  theme(legend.title = element_blank())




#Rob_check az indexekre
data_unc <- drop_na(data_unc)

data_macro <- data_macro %>% 
  inner_join(seq(from = max(min(data_macro$date), min(data_unc$date)), 
                 to = min(max(data_macro$date), max(data_unc$date)),
                 by = 1) %>% 
               as_tibble() %>% 
               mutate(day = day(value)) %>% filter(day == 1) %>% select(-day) %>% 
               rename(date = value),
             by = "date")

data_unc <- data_unc %>% 
  inner_join(seq(from = max(min(data_macro$date), min(data_unc$date)), 
                 to = min(max(data_macro$date), max(data_unc$date)),
                 by = 1) %>% 
               as_tibble() %>% 
               mutate(day = day(value)) %>% filter(day == 1) %>% select(-day) %>% 
               rename(date = value),
             by = "date")




colnames <- colnames(data_unc %>% select(-date))
legendnames <- c("TED spread", "EPU", "VIX")
index_list <- NULL
for(i in 1:length(colnames)){
  index_list[[i]] <-
    sign_tvar_3regime(data = data_macro %>% 
                        inner_join(create_indicator3(data_unc, colnames[i], 0.75, 0.25), by = "date") %>% 
                        select(-date, -starts_with("indicator")), 
                      indicator = create_indicator3(data_unc, colnames[i], 0.75, 0.25) %>% 
                        inner_join(data_macro, by = "date") %>% select(starts_with("indicator")), 
                      constr_mat = rbind(c(+1, -2, -3) ,
                                         c(+1, -2, -3),
                                         c(+1, -2, -3)), seed = 2022,
                      nlags = 3, draws = 20000, subdraws = 200, nkeep = 1000, 
                      KMIN = 1, KMAX = 6, constant = TRUE, steps = 21)
}
irf_list <- NULL
for(i in 1:length(colnames)){
  irf_list[[i]] <-
    get_all_irfs(mod_list = index_list[[i]], 
                 varname = c("FFR", "CPI", "INDPRO"),
                 regimename = c("High uncertainty", "Medium uncertainty", "Low uncertainty"),
                 prob = 0.16,
                 irflength = 21) %>% 
    left_join(get_all_irfs(mod_list = index_list[[i]], 
                           varname = c("FFR", "CPI", "INDPRO"),
                           regimename = c("High uncertainty", "Medium uncertainty", "Low uncertainty"),
                           prob = 0.16,
                           irflength = 21) %>%
                filter(t == 0, label == "FFR") %>% 
                mutate(mult = 1/median) %>% 
                select(regime, mult), 
              by = "regime") %>% 
    mutate(label = factor(label, levels = c("FFR", "CPI", "INDPRO")),
           regime = factor(regime, levels = c("High uncertainty", "Medium uncertainty", "Low uncertainty"))) %>%
    mutate(median = median * mult,
           upper = upper * mult, 
           lower = lower * mult) %>% 
    mutate(index = legendnames[i])
}
bind_rows(irf_list) %>% 
  ggplot(aes(x = t, y = median, color = index, group = index)) +
  geom_line(size = 1) +
  geom_hline(yintercept = 0, color = "red") +
  facet_wrap(regime~label, scales = "free")+
  scale_color_viridis_d(begin = 0.05, end = 0.8, name = "Index", option = "viridis")+
theme_minimal() +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11)) +
  labs(x = "",
       y ="") +
  theme(legend.title = element_blank())











