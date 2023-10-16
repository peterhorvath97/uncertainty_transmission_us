prep <- function(){
library(fredr)
library(tidyverse)
library(rlang)
library(lubridate)
source('US/codes/data.R')
source('US/codes/functions.R')

data[['macro']] <- data[['macro']] %>% 
  mutate(pi = pi - lag(pi),
         y = y - lag(y)) %>% 
  drop_na()

data[['uncertainty']] <- data[['uncertainty']] %>% 
  full_join(data[['forrob']])

data[['forrob']] <- NULL

prep1 <- function(data, index, threshold, controls = NA){
#Define specs
data_prepped <- NULL

#Main 
data_prepped[['main']] <- 
  list(
    data = 
      data[['macro']] %>% 
      inner_join(create_indicator2(data[['uncertainty']], index, threshold), by = "date") %>% 
      select(-indicator),
    indicator = 
      create_indicator2(data[['uncertainty']], index, threshold) %>% 
      inner_join(data[['macro']], by = "date") %>% 
      select(date, indicator),
    controls = controls
  )

#Lag
data_prepped[['lag']] <- 
  list(
    data = 
      data[['macro']] %>% 
      inner_join(create_indicator2(data[['uncertainty']], index, threshold), by = "date") %>% 
      select(-indicator),
    indicator = 
      create_indicator2(data[['uncertainty']], index, threshold) %>% 
      inner_join(data[['macro']], by = "date") %>% 
      mutate(indicator = lag(indicator)) %>% 
      select(date, indicator),
    controls = controls
  )

#1y interest rate swap
data_prepped[['r1y']] <- 
  list(
    data = 
      data[['macro']] %>% 
      select(-r) %>% 
      left_join(data[['macro_alt']] %>% 
                  select(date, r = r1y)) %>% 
      select(date, r, pi, y) %>% 
      inner_join(create_indicator2(data[['uncertainty']], index, threshold), by = "date") %>% 
      select(-indicator),
    indicator = 
      create_indicator2(data[['uncertainty']], index, threshold) %>% 
      inner_join(data[['macro']], by = "date") %>% 
      select(date, indicator),
    controls = controls
  )

#shadow rate swap
data_prepped[['shadow']] <- 
  list(
    data = 
      data[['macro']] %>% 
      select(-r) %>% 
      left_join(data[['macro_alt']] %>% 
                  select(date, r = ffr_shadow)) %>% 
      select(date, r, pi, y) %>% 
      inner_join(create_indicator2(data[['uncertainty']], index, threshold), by = "date") %>% 
      select(-indicator),
    indicator = 
      create_indicator2(data[['uncertainty']], index, threshold) %>% 
      inner_join(data[['macro']], by = "date") %>% 
      select(date, indicator),
    controls = controls
  )



data_prepped
}
prep_controls <- function(data, index, threshold){

c1 <- prep1(data, index, threshold) 
c2 <- c1
c3 <- c1
c4 <- c1

for(i in 1:length(c1)){
c2[[i]][['controls']] <- c2[[i]][['data']] %>% 
  left_join(data[['uncertainty']] %>% 
              select(date, sym(index))) %>% 
  inner_join(create_indicator2(data[['uncertainty']], index, threshold), by = "date") %>% 
  select(-indicator, -r, -pi, -y)

c3[[i]][['controls']] <- c3[[i]][['data']] %>% 
  left_join(data[['fin_stress']]) %>% 
  inner_join(create_indicator2(data[['uncertainty']], index, threshold), by = "date") %>% 
  select(-indicator, -r, -pi, -y)


c4[[i]][['controls']] <- c4[[i]][['data']] %>% 
  left_join(data[['fin_stress']]) %>% 
  left_join(data[['uncertainty']] %>% 
              select(date, sym(index))) %>% 
  inner_join(create_indicator2(data[['uncertainty']], index, threshold), by = "date") %>% 
  select(-indicator, -r, -pi, -y)
}

controls <- list(none = c1, 
     uncertainty = c2,
     financial_stress = c3, 
     both = c4)

controls
}

prep_thresholds <- function(data, index, thresholds){
  thresholds <- thresholds
  x <- NULL
  for(i in 1:length(thresholds)){
    x[[i]] <-  suppressMessages(prep_controls(data, index, thresholds[i]))
  }
  names(x) <- thresholds
  
  x
}
prep_indices <- function(data, thresholds){
  indices <- data[['uncertainty']] %>% 
    select(-date) %>% 
    colnames()
  
  x <- NULL
  for(i in 1:length(indices)){
    x[[i]] <- prep_thresholds(data, indices[[i]], c(0.5, 0.7))
  }
  names(x) <- indices
  
  x
}

prep_indices(data, c(0.5, 0,7))


}


data <- suppressMessages(prep())

