#packages
suppressPackageStartupMessages({
  library(readxl)
  library(tidyverse)
  library(tseries)
  library(forecast)
  library(TSA)
  library(fpp2)
  library(lmtest)
  library(vars)
  library(mFilter)
  library(ggplot2)
  library(tsDyn)
  library(fredr)
  library(VARsignR)
  library(lubridate)
  library(cowplot)
  library(gridExtra)
  library(patchwork)
  library(data.table)
  library(tsibble)
  library(mvnfast)
})

folder <- "C:/Users/horva/Documents/Egyetem/PhD/Research Material/Uncertainty project/?j mappa/"


#set api key
fredr_set_key("cda47ae66b38ed7988c0a9c2ec80c94f")

#download data
params <- list(
  series_id = c("USEPUINDXD", "DFF", "USACPIALLMINMEI", "INDPRO"),
  frequency = "m",
  observation_start = as.Date("1950-01-01")
)


import  <- pmap_dfr(
  .l = params,
  .f = ~ fredr(series_id = .x, frequency = .y)
) %>%
  dplyr::select(date, series_id, value) %>%
  spread(key = series_id, value = value) %>%
  drop_na() %>% rename(epu = USEPUINDXD,
                                        ffr = DFF,
                                        cpi = USACPIALLMINMEI,
                                        indpro = INDPRO
                                        ) 

#select data range
data <- import %>%
  dplyr::select(date, epu, ffr, cpi, indpro) %>%
  drop_na() %>% filter(date <= as.Date("2019-12-01"))

#plotting the data series
data %>%
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

#generate first differences / growth rates  
data <- data %>%
  mutate(cpi = cpi - lag(cpi),
         indpro = indpro - lag(indpro)) %>% drop_na()

#creating threshold values
data <- data %>%
  mutate(epu50 = median(epu),
         epu75 = quantile(epu, probs = 0.75),
         epu25 = quantile(epu, probs = 0.25))

#plot epu index with threshold values
data %>% dplyr::select(date, epu, epu50, epu25, epu75) %>%
  ggplot(aes(x = date, y = epu)) +
  geom_line() +
  geom_hline(aes(yintercept = epu50), color = "darkred") +
  geom_hline(aes(yintercept = epu25), color = "darkgreen") +
  geom_hline(aes(yintercept = epu75), color = "darkgreen") +
  theme_minimal() +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11)) +
  labs(x = "",
       y = "",
       title = "BBD economic policy uncertainty index")

#function to retrieve impulse responses
get_irf <- function(dim, prob, irflength, irfdata){
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

#some parameters for estimation and irfs
nlag <- 3
seed <- 2022
prob <- 0.16
irflength <- 21
draw <- 20000
subdraw <- 200

#plot variables coloured by regime indicator
p2 <- data %>%
  dplyr::select(date, ffr, cpi, indpro, epu, epu50) %>%
  mutate(indicator = ifelse(epu > epu50, 1, 2),
         indicator = case_when(indicator == 1 ~ "High uncertainty regime",
                               indicator == 2 ~ "Low uncertainty regime")) %>%
  gather(key = "variable", value = "value", ffr, cpi, indpro) %>%
  mutate(variable = case_when(variable == "ffr" ~ "Federal Funds Rate",
                              variable == "cpi" ~ "Consumer Prices",
                              variable == "indpro" ~ "Industrial Production Index"),
         variable = factor(variable, levels = c("Federal Funds Rate", "Consumer Prices", "Industrial Production Index"))) %>%
  ggplot(aes(x = date, y = value, color = indicator, group = 1)) +
  geom_line(size = 1) +
  facet_wrap(~variable, scales = "free", nrow = 3) + theme_minimal() +
  labs(x = "",
       y = "") + 
  theme(legend.title = element_blank()) 

# fit 2 regime model
data_2regime <- data %>%
  dplyr::select(date, ffr, cpi, indpro, epu, epu50) %>%
  mutate(indicator = ifelse(epu > epu50, 1, 0),
         indicator_inv = ifelse(indicator == 1, 0, 1),
         ffr_up = ffr*indicator,
         cpi_up = cpi*indicator,
         indpro_up = indpro*indicator,
         ffr_dn = ffr*indicator_inv,
         cpi_dn = cpi*indicator_inv,
         indpro_dn = indpro*indicator_inv) %>% as_tibble() %>%
  dplyr::select(ffr_up, cpi_up, indpro_up, ffr_dn, cpi_dn, indpro_dn) %>% ts()


set.seed(seed)
constr_up <- c(+1, -2, -3)
model_up <- rwz.reject(Y=data_2regime, nlags=nlag, draws=draw, subdraws=subdraw, nkeep=1000, KMIN=1,
                       KMAX=6, constrained=constr_up, constant=TRUE, steps=21)
irfs1 <- model_up$IRFS

set.seed(seed)
constr_dn <- c(+4, -5, -6)
model_dn <- rwz.reject(Y=data_2regime, nlags=nlag, draws=draw, subdraws=subdraw, nkeep=1000, KMIN=1,
                       KMAX=6, constrained=constr_dn, constant=TRUE, steps=21)
irfs2 <- model_dn$IRFS



ffr_up <- get_irf(1, prob, irflength, irfs1) %>% mutate(label = "FFR",
                                                        regime = "High uncertainty")
cpi_up <- get_irf(2, prob, irflength, irfs1) %>% mutate(label = "CPI",
                                                        regime = "High uncertainty")
indpro_up <- get_irf(3, prob, irflength, irfs1) %>% mutate(label = "INDPRO",
                                                           regime = "High uncertainty")

ffr_dn <- get_irf(4, prob, irflength, irfs2) %>% mutate(label = "FFR",
                                                        regime = "Low uncertainty")
cpi_dn <- get_irf(5, prob, irflength, irfs2) %>% mutate(label = "CPI",
                                                        regime = "Low uncertainty")
indpro_dn <- get_irf(6, prob, irflength, irfs2) %>% mutate(label = "INDPRO",
                                                           regime = "Low uncertainty")

bind_rows(ffr_up, cpi_up, indpro_up,
          ffr_dn, cpi_dn, indpro_dn) %>%
  mutate(label = factor(label, levels = c("FFR", "CPI", "INDPRO")),
         regime = factor(regime, levels = c("High uncertainty", "Low uncertainty"))) %>% 
  left_join(bind_rows(ffr_up, cpi_up, indpro_up,
                      ffr_dn, cpi_dn, indpro_dn) %>%
              mutate(label = factor(label, levels = c("FFR", "CPI", "INDPRO")),
                     regime = factor(regime, levels = c("High uncertainty", "Low uncertainty"))) %>% 
              filter(t == 0) %>% 
              mutate(mult = 1/median) %>% 
              mutate(mult = ifelse(mult < 0, NA, mult)) %>% 
              group_by(regime) %>% 
              mutate(mult = min(mult, na.rm = TRUE)) %>% 
              ungroup() %>% 
              mutate(median = median * mult,
                     upper = upper * mult,
                     lower = lower * mult) %>% 
              dplyr::select(label, regime, mult), 
            by = c("regime", 'label')) %>% 
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

#fit 3 regime model
data_3regime <- data %>%
  dplyr::select(date, ffr, cpi, indpro, epu, epu25, epu75) %>%
  mutate(indicator_up = ifelse(epu > epu75, 1,0),
         indicator_dn = ifelse(epu < epu25, 1, 0),
         indicator_mid = ifelse(indicator_up == 0 & indicator_dn == 0, 1, 0),
         ffr_up = ffr*indicator_up,
         cpi_up = cpi*indicator_up,
         indpro_up = indpro*indicator_up,
         ffr_mid = ffr*indicator_mid,
         cpi_mid = cpi*indicator_mid,
         indpro_mid = indpro*indicator_mid,
         ffr_dn = ffr*indicator_dn,
         cpi_dn = cpi*indicator_dn,
         indpro_dn = indpro*indicator_dn) %>% 
  as_tibble() %>%
  dplyr::select(ffr_up, cpi_up, indpro_up, ffr_mid, cpi_mid, indpro_mid, ffr_dn, cpi_dn, indpro_dn) %>% ts()


set.seed(seed)
constr_up <- c(+1, -2, -3)
model_up <- rwz.reject(Y=data_3regime, nlags=nlag, draws=draw, subdraws=subdraw, nkeep=1000, KMIN=1,
                       KMAX=6, constrained=constr_up, constant=TRUE, steps=21)
irfs1 <- model_up$IRFS

set.seed(seed)
constr_mid <- c(+4, -5, -6)
model_mid <- rwz.reject(Y=data_3regime, nlags=nlag, draws=draw, subdraws=subdraw, nkeep=1000, KMIN=1,
                        KMAX=6, constrained=constr_mid, constant=TRUE, steps=21)
irfs2 <- model_mid$IRFS

set.seed(seed)
constr_dn <- c(+7, -8, -9)
model_dn <- rwz.reject(Y=data_3regime, nlags=nlag, draws=draw, subdraws=subdraw, nkeep=1000, KMIN=1,
                       KMAX=6, constrained=constr_dn, constant=TRUE, steps=21)
irfs3 <- model_dn$IRFS


prob <- 0.16
irflength <- 21

ffr_up <- get_irf(1, prob, irflength, irfs1) %>% mutate(label = "FFR",
                                                        regime = "High uncertainty")
cpi_up <- get_irf(2, prob, irflength, irfs1) %>% mutate(label = "CPI",
                                                        regime = "High uncertainty")
indpro_up <- get_irf(3, prob, irflength, irfs1) %>% mutate(label = "INDPRO",
                                                           regime = "High uncertainty")

ffr_mid <- get_irf(4, prob, irflength, irfs2) %>% mutate(label = "FFR",
                                                         regime = "Medium uncertainty")
cpi_mid <- get_irf(5, prob, irflength, irfs2) %>% mutate(label = "CPI",
                                                         regime = "Medium uncertainty")
indpro_mid <- get_irf(6, prob, irflength, irfs2) %>% mutate(label = "INDPRO",
                                                            regime = "Medium uncertainty")

ffr_dn <- get_irf(7, prob, irflength, irfs3) %>% mutate(label = "FFR",
                                                        regime = "Low uncertainty")
cpi_dn <- get_irf(8, prob, irflength, irfs3) %>% mutate(label = "CPI",
                                                        regime = "Low uncertainty")
indpro_dn <- get_irf(9, prob, irflength, irfs3) %>% mutate(label = "INDPRO",
                                                           regime = "Low uncertainty")


bind_rows(ffr_up, cpi_up, indpro_up,
          ffr_mid, cpi_mid, indpro_mid,
          ffr_dn, cpi_dn, indpro_dn) %>%
  mutate(label = factor(label, levels = c("FFR", "CPI", "INDPRO")),
         regime = factor(regime, levels = c("High uncertainty", "Medium uncertainty", "Low uncertainty"))) %>% 
  left_join(bind_rows(ffr_up, cpi_up, indpro_up,
                      ffr_mid, cpi_mid, indpro_mid,
                      ffr_dn, cpi_dn, indpro_dn) %>%
              mutate(label = factor(label, levels = c("FFR", "CPI", "INDPRO")),
                     regime = factor(regime, levels = c("High uncertainty", "Medium uncertainty", "Low uncertainty"))) %>% 
              filter(t == 0) %>% 
              mutate(mult = 1/median) %>% 
              mutate(mult = ifelse(mult < 0, NA, mult)) %>% 
              group_by(regime) %>% 
              mutate(mult = min(mult, na.rm = TRUE)) %>% 
              ungroup() %>% 
              mutate(median = median * mult,
                     upper = upper * mult,
                     lower = lower * mult) %>% 
              dplyr::select(label, regime, mult), 
            by = c("regime", 'label')) %>% 
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













