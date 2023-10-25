library(fredr)
library(tidyverse)
library(rlang)
library(lubridate)
library(lmtest)

source('US/codes/data.R')
invisible(paste('US/codes/varsignr/', list.files('US/codes/varsignr') , sep = '') %>% lapply(source))


residuals <- read_rds('US/data/estout_clean.Rds')[['innovs']]

uncertainty <- data[['uncertainty']]
base_lin_data <- data[['macro']] 

base_lin <- rwz.reject(base_lin_data %>% 
                         select(-date) %>% 
                         ts(),
                       constrained = c(+1, -2, -3),
                       nlags = 4, 
                       draws = 20000, 
                       subdraws = 200, 
                       nkeep = 1000, 
                       KMIN = 1, 
                       KMAX = 6, 
                       constant = TRUE, 
                       steps = 21)

res_lin <- bind_cols(
  base_lin_data %>% 
  select(date) %>% 
  slice(-(1:4)), 
  t(apply(base_lin[['SHOCKS']], 
          2, 
          quantile, 
          probs = c(0.5, .16, 1-.16))) %>% 
    as_tibble() %>% 
    rename(median = `50%`,
           lower = `16%`,
           upper = `84%`)
)




granger_clean <- function(residuals, uncertainty){
df <- inner_join(residuals, uncertainty)

shock <-df %>% 
  select(median) %>% 
  ts() 

unc <- df %>% 
  select(epu) %>% 
  ts()

#Does the shock granger cause uncertainty
test_shocks <- lmtest::grangertest(shock,
                    unc,
                    order = 4)
test_shocks <- test_shocks['Pr(>F)'] %>% 
  as_tibble() %>% 
  drop_na() %>% 
  rename(p = `Pr(>F)`) %>% 
  mutate(test = 'Test 1')

#Does uncertainty granger cause the shock
test_unc <- lmtest::grangertest(unc,
                    shock,
                    order = 4)
test_unc<- test_unc['Pr(>F)'] %>% 
  as_tibble() %>% 
  drop_na() %>% 
  rename(p = `Pr(>F)`) %>% 
  mutate(test = 'Test 2')

out <- bind_rows(test_unc, test_shocks)

out
}

granger_out <- bind_rows(
granger_clean(residuals = res_lin,
              uncertainty = uncertainty) %>% 
  spread(test, p) %>% 
  mutate(model = 'Linear'),

granger_clean(residuals = residuals %>% 
                filter(spec == 'main',
                       control == 'none',
                       thresh == '0.5',
                       index == 'epu',
                       restr == 'full'),
              uncertainty = uncertainty) %>% 
  spread(test, p) %>% 
  mutate(model = 'Base TVAR, Median'),

granger_clean(residuals = residuals %>% 
                filter(spec == 'main',
                       control == 'none',
                       thresh == '0.7',
                       index == 'epu',
                       restr == 'full'),
              uncertainty = uncertainty) %>% 
  spread(test, p) %>% 
  mutate(model = 'Base TVAR, 70th perc.'),

granger_clean(residuals = residuals %>% 
                filter(spec == 'main',
                       control == 'both',
                       thresh == '0.5',
                       index == 'epu',
                       restr == 'full'),
              uncertainty = uncertainty) %>% 
  spread(test, p) %>% 
  mutate(model = 'Full control TVAR, Median'),

granger_clean(residuals = residuals %>% 
                filter(spec == 'main',
                       control == 'both',
                       thresh == '0.7',
                       index == 'epu',
                       restr == 'full'),
              uncertainty = uncertainty) %>% 
  spread(test, p) %>% 
  mutate(model = 'Full control TVAR, 70th perc.'),

granger_clean(residuals = residuals %>% 
                filter(spec == 'shadow',
                       control == 'both',
                       thresh == '0.5',
                       index == 'epu',
                       restr == 'full'),
              uncertainty = uncertainty) %>% 
  spread(test, p) %>% 
  mutate(model = 'Full control TVAR with shadow rate, Median'),

granger_clean(residuals = residuals %>% 
                filter(spec == 'shadow',
                       control == 'both',
                       thresh == '0.7',
                       index == 'epu',
                       restr == 'full'),
              uncertainty = uncertainty) %>% 
  spread(test, p) %>% 
  mutate(model = 'Full control TVAR with shadow rate, 70th perc.')
) %>% 
  select(Model = model,
         `Test 1`,
         `Test 2`) %>% 
  mutate(across(-Model, ~round(.x, digits = 3)))


saveRDS(granger_out, 'US/data/granger.rds')

