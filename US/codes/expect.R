library(tidyverse)
library(fredr)

fredr_set_key("cda47ae66b38ed7988c0a9c2ec80c94f")

epu <- fredr('USEPUINDXD') %>% 
  select(date, epu = value)

expinf_cleveland <- fredr('EXPINF1YR') %>% 
  select(date, expinf = value)

expinf_mich <- fredr('MICH') %>% 
  select(date, expinf = value)


test1 <- inner_join(epu, expinf_cleveland) %>% 
  drop_na() %>% 
  mutate(med = median(epu),
         regime = ifelse(epu > med, 'High Uncertainty', 'Low Uncertainty')) %>% 
  select(regime, expinf) %>% 
  var.test(expinf ~ regime, data = .,
           alternative = 'less')


test2 <- inner_join(epu, expinf_mich) %>% 
  drop_na() %>% 
  mutate(med = median(epu),
         regime = ifelse(epu > med, 'High Uncertainty', 'Low Uncertainty')) %>% 
  select(regime, expinf) %>% 
  var.test(expinf ~ regime, data = .,
           alternative = 'less')

test1
test2


inner_join(epu, expinf_mich) %>% 
  drop_na() %>% 
  mutate(med = median(epu),
         regime = ifelse(epu > med, 'High Uncertainty', 'Low Uncertainty')) %>% 
  select(regime, expinf) %>% 
  t.test(expinf ~ regime, data = .,
           alternative = 'less')


inner_join(epu, expinf_cleveland) %>% 
  drop_na() %>% 
  mutate(med = median(epu),
         regime = ifelse(epu > med, 'High Uncertainty', 'Low Uncertainty')) %>% 
  select(regime, expinf) %>% 
  t.test(expinf ~ regime, data = .,
           alternative = 'less')


