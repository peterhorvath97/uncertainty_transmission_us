#loading packages
#packages
library(fredr)
library(tidyverse)
library(gridExtra)
library(fredr)
library(VARsignR)
library(rlang)
library(lubridate)

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
  select(date, series_id, value) %>%
  spread(key = series_id, value = value) %>%
  drop_na() %>% rename(epu = USEPUINDXD,
                       ffr = DFF,
                       cpi = USACPIALLMINMEI,
                       indpro = INDPRO
  ) 

#select data range
data <- import %>%
  select(date, epu, ffr, cpi, indpro) %>%
  drop_na() %>% filter(date <= as.Date("2019-12-01"))



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

#plot the threshold values
data %>% 
  select(date, starts_with("epu")) %>% 
  gather(key = "variable", value = "value", starts_with("epu")) %>% 
  mutate(variable = case_when(variable == "epu" ~ "EPU",
                              variable == "epu50" ~ "Median",
                              variable == "epu25" ~ "Lower quartile",
                              variable == "epu75" ~ "Upper quartile")) %>% 
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



unc_measures <- bind_rows(
  fredr(series_id = "USEPUINDXD",
        frequency = "m") %>%
    dplyr::select(date, value) %>%
    mutate(id = "EPU index",
           id = as.factor(id)),
  fredr(series_id = "TEDRATE",
        frequency = "m") %>%
    dplyr::select(date, value)%>%
    mutate(id = "TED spread",
           id = as.factor(id)),
  fredr(series_id = "VIXCLS",
        frequency = "m") %>%
    dplyr::select(date, value)%>%
    mutate(id = "VIX index",
           id = as.factor(id))
) %>% filter(date <= as.Date("2019-12-01"))

#Plot uncertainty measures
unc_measures %>% 
  ggplot(aes(x = date, y = value, group = id, color = factor(id))) +
  geom_line() + 
  geom_smooth(method = "gam", se = FALSE, linetype = "dashed") +
  scale_y_log10() +
  scale_color_viridis_d(begin = 0.05, end = 0.4, name = "Index", option = "turbo") +
  theme_minimal() +
  labs(x = "",
       y = "")

unc_measures <- unc_measures %>% 
  mutate(id = case_when(id == "EPU index" ~ "epu",
                        id == "TED spread" ~ "ted",
                        id == "VIX index" ~ "vix")) %>% 
  pivot_wider(names_from = id, values_from = value) %>% 
  drop_na()

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
  rownames(cortab) <- c("EPU", "TED", "VIX")
  colnames(cortab) <- c("EPU", "TED", "VIX")
  cortab
}


cortab_create(unc_measures, 3)


#variables by regime indicator
#####
#plot variables coloured by regime indicator - 2
p2 <- data %>%
  select(date, ffr, cpi, indpro, epu, epu50) %>%
  mutate(indicator = ifelse(epu > epu50, 1, 2),
         indicator = case_when(indicator == 1 ~ "High uncertainty regime",
                               indicator == 2 ~ "Low uncertainty regime")) %>%
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
              filter(recession_start >= min(data$date),
                     recession_start <= max(data$date)),
            inherit.aes = F,
            aes(xmin = recession_start, xmax = recession_end, ymin = -Inf, ymax = Inf), 
            fill = "grey50", alpha = 0.5)

#plot variables coloured by regime indicator - 3
p3 <- data %>%
  select(date, ffr, cpi, indpro, epu, epu25, epu75) %>%
  mutate(indicator_up = ifelse(epu > epu75, 1,0),
         indicator_dn = ifelse(epu < epu25, 1, 0),
         indicator_mid = ifelse(indicator_up == 0 & indicator_dn == 0, 1, 0),
         indicator = case_when(indicator_up == 1 ~ "High uncertainty regime",
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
              filter(recession_start >= min(data$date),
                     recession_start <= max(data$date)),
            inherit.aes = F,
            aes(xmin = recession_start, xmax = recession_end, ymin = -Inf, ymax = Inf), 
            fill = "grey50", alpha = 0.5)
#####
grid.arrange(p2, p3)


#some parameters for estimation and irfs
nlag <- 3
seed <- 2022
prob <- 0.16
irflength <- 21
draw <- 20000
subdraw <- 200

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





# fit 2 regime model
data_2regime <- data %>%
  select(date, ffr, cpi, indpro, epu, epu50) %>%
  mutate(indicator = ifelse(epu > epu50, 1, 0),
         indicator_inv = ifelse(indicator == 1, 0, 1),
         ffr_up = ffr*indicator,
         cpi_up = cpi*indicator,
         indpro_up = indpro*indicator,
         ffr_dn = ffr*indicator_inv,
         cpi_dn = cpi*indicator_inv,
         indpro_dn = indpro*indicator_inv) %>% as_tibble() %>%
  select(ffr_up, cpi_up, indpro_up, ffr_dn, cpi_dn, indpro_dn) %>% ts()

est_irfs_2regime <- function(data_2regime, nlag, draw, subdraw, irflength, prob, seed){
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
              filter(t == 0, label == "FFR") %>% 
              mutate(mult = 1/median) %>% 
              select(regime, mult), 
            by = "regime") %>% 
  mutate(median = median * mult,
         upper = upper * mult, 
         lower = lower * mult)
}

irfs_2regime <- est_irfs_2regime(data_2regime, nlag, draw, subdraw, irflength, prob, seed)

irfs_2regime %>% 
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
  select(date, ffr, cpi, indpro, epu, epu25, epu75) %>%
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
  select(ffr_up, cpi_up, indpro_up, ffr_mid, cpi_mid, indpro_mid, ffr_dn, cpi_dn, indpro_dn) %>% ts()

est_irfs_3regime <- function(data_3regime, nlag, draw, subdraw, irflength, prob, seed){
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
              filter(t == 0, label == "FFR") %>% 
              mutate(mult = 1/median) %>% 
              select(regime, mult), 
            by = "regime") %>% 
  mutate(median = median * mult,
         upper = upper * mult, 
         lower = lower * mult)

}

irfs_3regime <- est_irfs_3regime(data_3regime, nlag, draw, subdraw, irflength, prob, seed)

irfs_3regime %>% 
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


#Robcheck for nlags
nlags_list <- NULL
for(i in 1:3){
  nlags_list[[i]] <- est_irfs_3regime(data_3regime, nlag = i, draw, subdraw, irflength, prob, seed) %>% 
    mutate(lag = paste("Lags = ", i, sep = ""))
}
bind_rows(nlags_list) %>% 
  ggplot(aes(x = t, y = median, color = lag, group = lag)) +
  geom_line(size = 1) +
  geom_hline(yintercept = 0, color = "red") +
  facet_wrap(regime~label, scales = "free")+
  scale_color_viridis_d(begin = 0.05, end = 0.8, name = "Index", option = "magma")+
  theme_minimal() +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11)) +
  labs(x = "",
       y ="") +
  theme(legend.title = element_blank())  

#Robcheck for index
data <- data %>% select(date, ffr, cpi, indpro) %>% 
  inner_join(unc_measures, by = "date")

data_list <- NULL
indices <- c("epu", "ted", "vix")
for(i in 1:length(indices)){
  data_list[[i]] <- 
    data %>% 
    select(ffr, cpi, indpro, !!sym(indices[i])) %>% 
    mutate(indicator_up = ifelse(!!sym(indices[i]) > quantile(!!sym(indices[i]), 0.75), 1, 0),
           indicator_dn = ifelse(!!sym(indices[i]) < quantile(!!sym(indices[i]), 0.25), 1, 0),
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
    select(ffr_up, cpi_up, indpro_up, ffr_mid, cpi_mid, indpro_mid, ffr_dn, cpi_dn, indpro_dn) %>% ts()
}
index_list <- NULL
for(i in 1:length(indices)){
  index_list[[i]] <- est_irfs_3regime(data_list[[i]], nlag, draw, subdraw, irflength, prob, seed) %>% 
    mutate(index = str_to_upper(indices[i]))
}
bind_rows(index_list) %>% 
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




