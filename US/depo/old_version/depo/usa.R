#packages
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

#set api key
fredr_set_key("cda47ae66b38ed7988c0a9c2ec80c94f")

#retrieve data from fred database
ffr <- fredr(series_id = "DFF",
             frequency = "m") %>%
  dplyr::select(date, value) %>%
  rename(ffr = value)

cpi <- fredr(series_id = "USACPIALLMINMEI",
             frequency = "m",
             units = "chg") %>%
  dplyr::select(date, value) %>%
  rename(cpi = value)

indpro <- fredr(series_id = "INDPRO",
                frequency = "m",
                units = "chg") %>%
  dplyr::select(date, value) %>%
  rename(indpro = value)

r3m <- fredr(series_id = "DTB3",
             frequency = "m") %>%
  dplyr::select(date, value) %>%
  rename(r3m = value)

r10y <- fredr(series_id = "DGS10",
              frequency = "m") %>%
  dplyr::select(date, value) %>%
  rename(r10y = value)

m1 <- fredr(series_id = "M1REAL",
            frequency = "m",
            units = "chg") %>%
  dplyr::select(date, value) %>%
  rename(m1 = value)

m2 <- fredr(series_id = "M2REAL",
            frequency = "m",
            units = "chg") %>%
  dplyr::select(date, value) %>%
  rename(m2 = value)

depint <- fredr(series_id = "IR3TCD01USM156N",
                frequency = "m") %>%
  dplyr::select(date, value) %>%
  fill(value, .direction = "downup")%>%
  rename(depint = value)

oil <- fredr(series_id = "WTISPLC",
             frequency = "m",
             units = "chg") %>%
  dplyr::select(date, value) %>%
  rename(oil = value)

totres <- fredr(series_id = "TRESEGUSM052N",
                frequency = "m",
                units = "chg") %>%
  dplyr::select(date, value) %>%
  mutate(value = value/1000000000) %>%
  rename(totres = value)

borrow <- fredr(series_id = "BORROW",
                frequency = "m",
                units = "chg") %>%
  dplyr::select(date, value) %>%
  rename(borrow = value)

#setup ts dataframe
datalist <- list(ffr, cpi, indpro, r3m, r10y, m1, m2, 
                 depint, oil, totres, borrow)
mindate <- vector()
maxdate <- vector()

for (i in 1:length(datalist)){
  mindate[i] <- min(datalist[[i]]$date)
  maxdate[i] <- max(datalist[[i]]$date)
}
mindate <- as.Date(max(mindate))
maxdate <- as.Date(min(maxdate))
maxdate <- as.Date("2019-12-01")

for (i in 1:length(datalist)){
  datalist[[i]] <- datalist[[i]] %>%
    filter(date>mindate,
           date<maxdate) 
}

vardata <- datalist[[1]]
for (i in 2:length((datalist))){
  vardata <- vardata %>% left_join(datalist[[i]], by = "date")
}
vardata <- vardata %>%
  dplyr::select(-date) %>%
  mutate(m2 = m2-m1,
         nonborrow = totres - borrow) %>%
  dplyr::select(-totres)%>%
  relocate(borrow, .after = ffr)%>%
  relocate(nonborrow, .after = borrow)%>%
  relocate(m1, .after = nonborrow) %>%
  relocate(m2, .after = m1) %>%
  relocate(depint, .after = m2)%>%
  relocate(r3m, .after = depint)%>% 
  relocate(r10y, .after = r3m)%>%
  relocate(indpro, .after = cpi)%>%
  relocate(oil, .before = cpi)%>%
  ts(start = format(mindate, "%Y"), end = format(maxdate, "%Y"), 
     frequency = 12)

vardata <- vardata %>%
  relocate(borrow, .after = ffr)%>%
  relocate(nonborrow, .after = borrow)%>%
  relocate(m1, .after = nonborrow) %>%
  relocate(m2, .after = m1) %>%
  relocate(depint, .after = m2)%>%
  relocate(r3m, .after = depint)%>% 
  relocate(r10y, .after = r3m)%>%
  relocate(indpro, .after = cpi)%>%
  relocate(oil, .before = cpi)%>%
  ts(start = format(mindate, "%Y"), end = format(maxdate, "%Y"), 
     frequency = 12)


rm(ffr, cpi, indpro, r3m, r10y, m1, m2, 
   depint, oil, totres, borrow, i, datalist, mindate, maxdate)

lagopt <- VARselect(vardata, lag.max = 24) %>%
  data.frame()
lagopt <- data.frame(t(lagopt$selection))
colnames(lagopt) <- c("AIC", "HQ", "SC", "FPE")
lagopt <- lagopt %>% dplyr::select(HQ) %>% as.numeric()


varest <- VAR(vardata, p = lagopt)


irf <- irf(varest, impulse = "ffr",
           response = c("cpi", "indpro"),
           n.ahead= 20, boot = TRUE)
plot(irf)


tvarest <- TVAR(data = vardata, lag = lagopt, nthresh = 1, 
                thVar = vardata[,ncol(vardata)], thDelay=1, trim=0.15, mTh=1)
plot(tvarest)

rm(lagopt)

#source("https://raw.githubusercontent.com/MatthieuStigler/tsDyn_GIRF/master/GIRF2")

#library(remotes)
#install_github("angusmoore/tvarGIRF", ref = "v0.1.1")
#library(tvarGIRF)


