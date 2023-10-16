library(tidyverse)
library(rlang)
library(lubridate)

estout <- read_rds('US/data/estout.Rds')

catch <- 'start'
for(m in 1:length(restrs)){
  restr <- restrs[m]
  for(n in 1:length(indices)){
    index <- indices[n]
    for(k in 1:length(thresholds)){
      thresh <- thresholds[k]
      for(j in 1: length(controlstruct)){
        cont <- controlstruct[j]
        for(i in 1:length(specs)){
          spec <- specs[i]
          valid <- estout[[restr]][[index]][[thresh]][[cont]][[spec]][[2]] %>% is.data.frame()
          if(valid){
            
          } else {
            estout[[restr]][[index]][[thresh]][[cont]][[spec]] <- NULL
            catch <- c(catch, paste(restr, index, thresh, cont, spec, sep = '-'))
          }
        }
        len <- (length(estout[[restr]][[index]][[thresh]][[cont]]) > 0)
        if(len){ } else {
          estout[[restr]][[index]][[thresh]][[cont]] <- NULL
        }
      }
      len <- (length(estout[[restr]][[index]][[thresh]]) > 0)
      if(len){ } else {
        estout[[restr]][[index]][[thresh]] <- NULL
      }
    }
    len <- (length(estout[[restr]][[index]]) > 0)
    if(len){ } else {
      estout[[restr]][[index]] <- NULL
    }
  }
  len <- (length(estout[[restr]]) > 0)
  if(len){ } else {
    estout[[restr]] <- NULL
  }
}


catch <- catch[2:length(catch)] %>% 
  as_tibble() %>% 
  separate(value,
           into = c('restr', 'index', 'threshold', 'control', 'spec'),
           sep = '-') %>% 
  select(-restr)


namesform <- names(estout)
for(m in 1:length(namesform)){
  restr <- namesform[m]
  
  irf_out3 <- NULL
  res_out3 <- NULL
  namesforn <- names(estout[[restr]])
  for(n in 1:length(namesforn)){
    index <- namesforn[n]
    
    irf_out2 <- NULL
    res_out2 <- NULL
    namesfork <- names(estout[[restr]][[index]])
    for(k in 1:length(namesfork)){
      thresh <- namesfork[k]
      
      irf_out1 <- NULL
      res_out1 <- NULL
      namesforj <- names(estout[[restr]][[index]][[thresh]])
      for(j in 1:length(namesforj)){
        cont <- namesforj[j]
        
        irfs <- NULL
        res <- NULL
        namesfori <- names(estout[[restr]][[index]][[thresh]][[cont]])
        for(i in 1: length(namesfori)){
          spec <- namesfori[i]
          irfs[[i]] <- estout[[restr]][[index]][[thresh]][[cont]][[spec]][[2]] %>% 
            mutate(spec = namesfori[i],
                   control = namesforj[j],
                   thresh = namesfork[k],
                   index = namesforn[n],
                   restr = namesform[m])
          res[[i]] <- estout[[restr]][[index]][[thresh]][[cont]][[spec]][[1]] %>% 
            mutate(spec = namesfori[i],
                   control = namesforj[j],
                   thresh = namesfork[k],
                   index = namesforn[n],
                   restr = namesform[m])
        }
        irf_out1[[j]] <- irfs %>% bind_rows()
        res_out1[[j]] <- res %>% bind_rows()
      }
      irf_out2[[k]] <- irf_out1 %>% bind_rows()
      res_out2[[k]] <- res_out1 %>% bind_rows()
    }
    irf_out3[[n]] <- irf_out2 %>% bind_rows()
    res_out3[[n]] <- res_out2 %>% bind_rows()
  }
  if(length(namesform) > 1){
    irf_out4[[m]] <- irf_out3 %>% bind_rows()
    res_out4[[m]] <- res_out3 %>% bind_rows()
  } else {
    irf_out4 <- irf_out3 %>% bind_rows()
    res_out4 <- res_out3 %>% bind_rows()
  }
  
}
if(length(namesform) > 1){
  irf_out <- irf_out4 %>% bind_rows()
  res_out <- res_out4 %>% bind_rows()
} else {
  irf_out <- irf_out4 
  res_out <- res_out4 
}

estout_clean <- list(irfs = irf_out,
                     innovs = res_out,
                     invalids = catch)

saveRDS(estout_clean, 'US/data/estout_clean.Rds')
