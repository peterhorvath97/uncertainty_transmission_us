library(fredr)
library(tidyverse)
library(rlang)
library(lubridate)
library(foreach)
library(doParallel)
library(doSNOW)

source('US_only/codes/dataprep.R')
source('US_only/codes/functions.R')
clean <- function(estdata){
  estdata$data <- estdata$data %>% 
    drop_na()
  estdata$indicator <- estdata$indicator %>% 
    drop_na()
  estdata$controls <- if(!is.logical(estdata$controls)){
    estdata$controls %>% 
      drop_na()} else{
        NA}
  

  
  mindate <- max(
    min(estdata$data$date),
    min(estdata$indicator$date),
    if(!is.logical(estdata$controls)){
      min(estdata$controls$date) 
    } else{
      NA
    },
    na.rm = T
  )
  
  maxdate <- min(
    max(estdata$data$date),
    max(estdata$indicator$date),
    if(!is.logical(estdata$controls)){
      max(estdata$controls$date) 
    } else{
      NA
    },
    na.rm = T
  )
  
  estdata$data <- estdata$data %>% 
    filter(date >= mindate,
           date <= maxdate) 
  estdata$indicator <- estdata$indicator %>% 
    filter(date >= mindate,
           date <= maxdate) 
  estdata$controls <- if(!is.logical(estdata$controls)){
    estdata$controls %>% 
      filter(date >= mindate,
             date <= maxdate)
  } else {
    NA
  }
  
  estdata
}

data2

#Sign restrictions - just follow Uhlig paper with 1-6
restr <- list(full = list(KMIN = 1,
                          KMAX = 6)#,
              #cont = list(KMIN = 1,
              #            KMAX = 1),
              #shape = list(KMIN = 2,
              #             KMAX = 6)
              )

restrs <- names(restr)
indices <- names(data)
#Indices - Focus on EPU only
indices <- indices[str_detect(indices, 'epu')]
thresholds <- names(data[[1]])
controlstruct <- names(data[[1]][[1]])
specs <- names(data[[1]][[1]][[1]])

ncores <- parallel::detectCores(logical = F)-1
cl <- makeCluster(ncores)
#registerDoParallel(ncores)
registerDoSNOW(cl)


estout <- NULL
for(m in 1:length(restrs)){
  outs3 <- NULL
  rest <- restrs[m]
for(n in 1:length(indices)){
outs2 <- NULL
index <- indices[n]
for(k in 1:length(thresholds)){
outs1 <- NULL
thresh <- thresholds[k]
for(j in 1:length(controlstruct)){
cont <- controlstruct[j]
pb <- txtProgressBar(max = length(specs), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
out <- 
  foreach(i = 1:length(specs),
          .packages = c('tidyverse'),
          .export = c("data", "sign_tvar_2regime_controls", "get_shock2", "get_all_irfs_clean", 'clean'),
          .errorhandling = 'pass',
          .options.snow = opts) %dopar% {
            
            estdata <- data[[index]][[thresh]][[cont]][[i]]
            estdata <- clean(estdata)
            est <- sign_tvar_2regime_controls(data = estdata$data %>% select(-date),
                                              indicator = estdata$indicator %>% select(-date),
                                              controls = if(!is.logical(estdata$controls)){
                                                estdata$controls %>% select(-date)
                                              } else {
                                                NA
                                              },
                                              constr_mat = rbind(c(+1, -2, -3) ,
                                                                 c(+1, -2, -3)), 
                                              seed = 2022,
                                              nlags = 4, 
                                              draws = 20000, 
                                              subdraws = 200, 
                                              nkeep = 1000, 
                                              KMIN = restr[[rest]][['KMIN']], 
                                              KMAX = restr[[rest]][['KMAX']], 
                                              constant = TRUE, 
                                              steps = 21)
            
            irfs <- get_all_irfs_clean(mod_list = est, 
                                       varname = c("r", "pi", "y"), 
                                       regimename = c("High uncertainty", "Low uncertainty"),
                                       prob = .16, 
                                       irflength = 21) %>% 
              mutate(label = toupper(label)) %>%
              mutate(label = factor(label, levels = c("R", "PI", "Y")),
                     regime = factor(regime, levels = c("High uncertainty", "Low uncertainty")))
            
            
            
            shocks <- get_shock2(model = est,
                                 prob = .16,
                                 indicatordata = estdata$indicator,
                                 modlags = 4)
            
            
            result <- list(shocks, irfs) 
          }
close(pb)
names(out) <- specs
outs1[[j]] <- out
}
names(outs1) <- controlstruct
outs2[[k]] <- outs1
}
names(outs2) <- thresholds
outs3[[n]] <- outs2
}
names(outs3) <- indices
estout[[m]] <- outs3
}
names(estout) <- restrs
stopCluster(cl)

saveRDS(estout, 'US_only/data/estout.Rds')

source('US_only/codes/estout_cleanup.R')
