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
  library(lpirfs)
})

folder <- "C:/Users/horva/Documents/Egyetem/PhD/Research Material/Uncertainty project/Új mappa/"


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
  drop_na() %>% as_tsibble() %>% rename(epu = USEPUINDXD,
                                        ffr = DFF,
                                        cpi = USACPIALLMINMEI,
                                        indpro = INDPRO
  ) 

#select data range
data <- import %>%
  dplyr::select(date, epu, ffr, cpi, indpro) %>%
  drop_na() %>% filter(date <= as.Date("2019-12-01")) %>%
  mutate(cpi = cpi - lag(cpi),
         indpro = indpro - lag(indpro)) %>% drop_na()



switching_data <- ifelse(data$epu > median(data$epu), 1, 0)

data <- data %>% as_tibble() %>% dplyr::select(-date, -epu)
  

results_lin <- lp_lin(endog_data = data,
                      lags_endog_lin = 1,
                      trend = 0,
                      shock_type = 1,
                      confint = 1.96,
                      hor = 20)
                      
plot(results_lin)                      
                      

results_nl <- lp_nl(data,
                    lags_endog_lin = 1, lags_endog_nl = 1,
                    trend = 0, shock_type = 1,
                    confint = 1.67, hor = 20,
                    switching = switching_data, lag_switching = FALSE,
                    use_logistic = FALSE)
plot(results_nl)                      
                      
                        
                        