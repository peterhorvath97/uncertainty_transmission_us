data <- function(){
#loading packages
#packages
require(fredr)
require(tidyverse)
require(rlang)
require(lubridate)
require(readxl)
  
  x13adj <- function(value, freq, date){
    require(forecast)
    require(seasonal)
    require(tidyverse)
    require(lubridate)
    tryCatch(
      value %>% 
        ts(frequency = freq,
           start = c(date %>% min %>% year, 
                     date %>% min %>% month)) %>% 
        seas() %>% 
        seasadj() %>% 
        as.numeric(),
      error = function(error){value}
    )
  }
  
fredr_set_key("cda47ae66b38ed7988c0a9c2ec80c94f")


macro_codes <- c(
  'DFF', #FFR
  'USACPIALLMINMEI', #CPI
  'GDPC1' #GDP
  )

unc_codes <- c(
  'USEPUINDXD', #EPU
  'EPUMONETARY', #EPU monetary
  'EPUFISCAL', #EPU fiscal
  'EPUTRADE', #EPU trade
  'EPUFINREG', #EPU financial regulation
  'EPUSOVDEBT', #EPU debt
  'EPUREG' #EPU regulation
)

macro_codes_alt <- c(
  'DGS1', #1 year rate
  'CORESTICKM159SFRBATL', #Sticky CPI (less food and energy)
  'PCEPI' #PCE
)

rob_codes <- c(
  'TEDRATE', #TED spread
  'WLEMUINDXD', #Equity market uncertainty
  'VIXCLS', #VIX index
  'WUIUSA' #WUI
)

fred_query <- function(ids, freq, start){
  
  require(fredr)
  require(dplyr)
  require(purrr)
  
  #download data
  params <- list(
    series_id = ids,
    frequency = freq,
    observation_start = as.Date(start)
  )
  
  
  data  <- pmap_dfr(
    .l = params,
    .f = ~ fredr(series_id = .x, frequency = .y) ) %>%
    dplyr::select(date, series_id, value) %>%
    spread(key = series_id, value = value) %>%
    drop_na() 
  
  data
}

macro <- fred_query(
  ids = macro_codes,
  freq = 'q',
  start = '1950-01-01') %>% 
  rename(r = DFF,
         pi = USACPIALLMINMEI,
         y = GDPC1) %>% 
  select(date, r, pi, y) %>% 
  mutate(date = as_date(date)) %>% 
  mutate(pi = x13adj(pi, 4, date))

macro <- macro %>% 
  mutate(
    year = year(date)) %>% 
  group_by(year) %>% 
  mutate(mean_y = mean(y),
         is.2015 = ifelse(year == 2015, 1, NA),
         mean_y = mean_y*is.2015) %>%
  ungroup() %>% 
  mutate(mean_y = mean(mean_y, na.rm = T),
         y = 100*y/mean_y) %>% 
  select(date, r, pi, y)
  
  

uncertainty <- fred_query(
  ids = unc_codes,
  freq = 'q',
  start = '1950-01-01') %>% 
  rename(epu = USEPUINDXD,
         epu_mon = EPUMONETARY,
         epu_fis = EPUFISCAL,
         epu_trade = EPUTRADE,
         epu_finreg = EPUFINREG,
         epu_debt = EPUSOVDEBT,
         epu_reg = EPUREG) %>% 
  mutate(date = as_date(date))

macro_alt <- fred_query(
  ids = macro_codes_alt,
  freq = 'q',
  start = '1950-01-01') %>% 
  rename(r1y = DGS1,
         pi_stick = CORESTICKM159SFRBATL,
         pi_pce = PCEPI) %>% 
  mutate(date = as_date(date)) %>% 
  left_join(read_excel('US/data/shadowrate.xlsx') %>% 
              drop_na() %>% 
              mutate(date = as_date(date),
                     year = year(date),
                     quarter = quarter(date)) %>% 
              group_by(year, quarter) %>% 
              mutate(ffr_shadow = mean(ffr_shadow)) %>% 
              ungroup() %>%
              distinct(year, quarter, .keep_all = T) %>% 
              select(date, ffr_shadow))


fin_stress <- fred_query(
  ids = 'STLFSI4',  #Financial stress
  freq = 'q',
  start = '1950-01-01') %>% 
  rename(fsi = STLFSI4) %>% 
  mutate(date = as_date(date))


forrob <- fred_query(
  ids = rob_codes,
  freq = 'q',
  start = '1950-01-01') %>% 
  rename(ted = TEDRATE,
         vix = VIXCLS,
         equnc = WLEMUINDXD,
         wui = WUIUSA) %>% 
  mutate(date = as_date(date)) %>% 
  left_join(read_excel("US/data/data_gpr_export.xls", 
                       col_types = c("date", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "text", "text")) %>% 
              select(date = month,
                     gpr = GPRC_USA) %>% 
              drop_na() %>% 
              mutate(date = as_date(date),
                     year = year(date),
                     quarter = quarter(date)) %>% 
              group_by(year, quarter) %>% 
              mutate(gpr = mean(gpr)) %>% 
              ungroup() %>%
              distinct(year, quarter, .keep_all = T) %>% 
              select(date, gpr)) %>% 
  left_join(    read_excel("US/data/tpu_web_latest.xlsx", 
                           sheet = "TPU_MONTHLY") %>% 
                  mutate(date = as_date(DATE),
                         year = year(date),
                         quarter = quarter(date)) %>% 
                  group_by(year, quarter) %>% 
                  mutate(tpu = mean(TPU)) %>% 
                  ungroup() %>%
                  distinct(year, quarter, .keep_all = T) %>% 
                  select(date, tpu))





mindate <- max(
  min(macro$date),
  min(uncertainty$date)
)

maxdate <- min(
  max(macro$date),
  max(uncertainty$date)
)



data <- mget(ls()) %>% keep(is.data.frame)

for(i in 1:length(data)){
  data[[i]] <- data[[i]] %>% 
    filter(date >= mindate,
           date <= maxdate) 
  }


data
}

data <- data() 