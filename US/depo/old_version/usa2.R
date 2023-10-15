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
library(VARsignR)
library(lubridate)
library(cowplot)
library(gridExtra)
library(patchwork)
library(data.table)
#https://cran.r-project.org/web/packages/VARsignR/vignettes/VARsignR-vignette.html

folder <- "C:/Users/horva/Documents/Egyetem/PhD/Research Material/Uncertainty project/Új mappa/"

#set api key
fredr_set_key("cda47ae66b38ed7988c0a9c2ec80c94f")

#retrieve data from fred database
#####
epu <- fredr(series_id = "USEPUINDXD",
             frequency = "m") %>%
  dplyr::select(date, value) %>%
  rename(epu = value)

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

vix <- fredr(series_id = "VIXCLS",
             frequency = "m") %>%
  dplyr::select(date, value) %>%
  rename(vix = value)

ted <- fredr(series_id = "TEDRATE",
             frequency = "m") %>%
  dplyr::select(date, value) %>%
  rename(ted = value)
#####


#setup ts dataframe
#####
datalist <- list(epu, ffr, cpi, indpro)
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


vardata <- vardata %>% dplyr::select(-epu, -date)%>%
  ts(start = format(mindate, "%Y"), end = format(maxdate, "%Y"), 
     frequency = 12)
#####

#testing
#####
#sign restriction for identifying shock
constr <- c(+1, -2, -3)

#model method 1 - reject
set.seed(2022)
model1 <- uhlig.reject(Y=vardata, nlags=1, draws=200, subdraws=200, nkeep=1000, KMIN=1,
                       KMAX=6, constrained=constr, constant=TRUE, steps=20)

irfs1 <- model1$IRFS

vl <- c("FFR", "CPI", "Indpro")

irfplot(irfdraws=irfs1, type="median", labels=vl, save=FALSE, bands=c(0.16, 0.84), 
        grid=TRUE, bw=FALSE)


#model method 2 - penalty
set.seed(2022)
model2 <- uhlig.penalty(Y=vardata, nlags=1, draws=200, subdraws=200, nkeep=1000, KMIN=1,
                       KMAX=6, constrained=constr, constant=TRUE, steps=20)

irfs2 <- model2$IRFS

vl <- c("FFR", "CPI", "Indpro")

irfplot(irfdraws=irfs2, type="median", labels=vl, save=FALSE, bands=c(0.16, 0.84), 
        grid=TRUE, bw=FALSE)
#####

#visualize data series
#####
ffrplot <- fredr(series_id = "DFF",
                 frequency = "m") %>%
  dplyr::select(date, value) %>%
  rename(ffr = value) %>% 
  filter(date >= mindate,
         date <= maxdate) %>%
  ggplot(aes(x = date, y = ffr)) +
  geom_line(color = "darkblue") +
  theme_classic() +
  labs(title = "Federal Funds Rate",
       y = "FFR",
       x = element_blank()) +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11)) 

cpiplot <- fredr(series_id = "USACPIALLMINMEI",
                 frequency = "m") %>%
  dplyr::select(date, value) %>%
  rename(cpi = value) %>%
  filter(date >= mindate,
         date <= maxdate) %>%
  ggplot(aes(x = date, y = cpi)) +
  geom_line(color = "darkblue") + 
  theme_classic() +
  labs(title = "Consumer Price Index",
       y = "CPI",
       x = element_blank()) +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11)) 

indproplot <- fredr(series_id = "INDPRO",
                    frequency = "m") %>%
  dplyr::select(date, value) %>%
  rename(indpro = value) %>%
  filter(date >= mindate,
         date <= maxdate) %>%
  ggplot(aes(x = date, y = indpro)) +
  geom_line(color = "darkblue") + theme_classic() +
  labs(title = "Industrial Production Index",
       y = "INDPRO",
       x = element_blank()) +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11)) 


#####

plot1 <- grid.arrange(ffrplot, cpiplot,
                      indproplot,  ncol = 1)
ggsave(paste(folder,"plot1.png",sep=""), plot1, bg = "transparent")

#compare uncertainty measures
#####
vix <- vix %>%
  mutate(id = "VIX index") %>%
  rename(value = vix)
ted <- ted %>%
  mutate(id = "TED spread") %>%
  rename(value = ted)

group.color <- c(A = "darkblue", B = "darkred", C = "darkgreen")

vis <- epu %>%
  mutate(id = "EPU index") %>%
  rename(value = epu) %>%
  bind_rows(vix) %>%
  bind_rows(ted) %>%
  filter(date >= mindate,
         date <= maxdate) %>%
  ggplot(aes(x = date, y = value, color = id, group = id)) + 
  geom_line() +
  theme_classic() +
  theme(legend.title = element_blank(),
        axis.title = element_blank()) +
  scale_y_log10() + 
  geom_smooth(method = "gam", se = FALSE, linetype = "dashed") +
  scale_fill_manual(values = group.color)

#####
ggsave(paste(folder,"plotcompare.png",sep=""), vis, bg = "transparent")

#visualize epu, indicating the high uncertainty areas
#####
epu %>% filter(date >= mindate,
               date <= maxdate) %>%
  ggplot(aes(x = date, y = epu)) +
  geom_line(color = "darkblue") +
  annotate("rect", xmin = as.Date("2008-01-01"), xmax = as.Date("2014-01-01"), ymin = 0.5, ymax = 260, 
           alpha = .2) +
  annotate("rect", xmin = as.Date("2007-01-01"), xmax = as.Date("2015-01-01"), ymin = 0.5, ymax = 260, 
           alpha = .2) +
  annotate("rect", xmin = as.Date("2000-08-01"), xmax = as.Date("2004-01-01"), ymin = 0.5, ymax = 260, 
           alpha = .2) + 
  annotate("rect", xmin = as.Date("1999-08-01"), xmax = as.Date("2005-01-01"), ymin = 0.5, ymax = 260, 
           alpha = .2) + 
  annotate("rect", xmin = as.Date("1990-06-01"), xmax = as.Date("1994-01-01"), ymin = 0.5, ymax = 260, 
           alpha = .2)+
  annotate("rect", xmin = as.Date("1989-06-01"), xmax = as.Date("1995-01-01"), ymin = 0.5, ymax = 260, 
           alpha = .2)+
  theme_classic() +
  labs(title = "Economic policy uncertainty index",
       y = "EPU",
       x = element_blank()) +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11)) 
#####
ggsave(paste(folder,"epu.png",sep=""), bg = "transparent")

##rolling window
#####
#setup
irflist <- list()
years <- 5
window <- years*12
nlag = 1
first <- 1
last <- 6
constr <- c(+1, -2, -3)
lastobs <- nrow(vardata)-window

#create rw estimations
for (i in 1:lastobs){
  subsample <- vardata[i:window, ] %>% ts()
  set.seed(2022)
  try(
  model <- uhlig.penalty(Y=subsample, nlags=nlag, draws=200, subdraws=200, nkeep=1000, KMIN=first,
                          KMAX=last, constrained=constr, constant=TRUE, steps=20)
  )
  
  irfs <- model$IRFS
  irflist[[i]] <- irfs
  
  window <- window+1
}
rm(subsample, irfs, window, model)

irflist2_l <- vector("list", length = length(irflist))
irflist2 <- vector("list", length = length(irflist))
irflist2_u <- vector("list", length = length(irflist))

x <- matrix(nrow = 20, ncol = 3)
y <- matrix(nrow = 20, ncol = 3)
z <- matrix(nrow = 20, ncol = 3)


for (k in 1:length(irflist)){
  for (j in 1:3){
    for (i in 1:20){
      x[i ,j] <- median(irflist[[k]][1:200, i ,j])
      y[i ,j] <- quantile(irflist[[k]][1:200, i ,j], probs = 0.86)
      z[i ,j] <- quantile(irflist[[k]][1:200, i ,j], probs = 0.14)
      irflist2[[k]] <- x
      irflist2_u[[k]] <- y
      irflist2_l[[k]] <- z
    }
  }
}
rm(x, y, z, irflist)
#####

#order of irfs: ffr, cpi, indpro

#subset the time periods in question
#####
t1 <- lubridate::interval(mindate, as.Date("2007-01-01")) %/% months(1)
t2 <- lubridate::interval(mindate, as.Date("2015-01-01")) %/% months(1) - years*12
int1 <- seq(from = t1, to = t2)
rm(t1, t2)

t3 <- lubridate::interval(mindate, as.Date("1999-08-01")) %/% months(1)
t4 <- lubridate::interval(mindate, as.Date("2005-01-01")) %/% months(1) - years*12
int2 <- seq(from = t3, to = t4)
rm(t3, t4)

t5 <- lubridate::interval(mindate, as.Date("1989-06-01")) %/% months(1)
t6 <- lubridate::interval(mindate, as.Date("1995-01-01")) %/% months(1) - years*12
int3 <- seq(from = t5, to = t6)
rm(t5, t6)


p1_l <- vector("list", length = length(int1))
p1 <- vector("list", length = length(int1))
p1_u <- vector("list", length = length(int1))


p2_l <- vector("list", length = length(int2))
p2 <- vector("list", length = length(int2))
p2_u <- vector("list", length = length(int2))


p3_l <- vector("list", length = length(int3))
p3 <- vector("list", length = length(int3))
p3_u <- vector("list", length = length(int3))



for (i in 1:length(int1)){
  p1_l[[i]] <- irflist2_l[[int1[i]]]
  p1[[i]] <- irflist2[[int1[i]]]
  p1_u[[i]] <- irflist2_u[[int1[i]]]
}
for (i in 1:length(int2)){
  p2_l[[i]] <- irflist2_l[[int2[i]]]
  p2[[i]] <- irflist2[[int2[i]]]
  p2_u[[i]] <- irflist2_u[[int2[i]]]
}
for (i in 1:length(int3)){
  p3_l[[i]] <- irflist2_l[[int3[i]]]
  p3[[i]] <- irflist2[[int3[i]]]
  p3_u[[i]] <- irflist2_u[[int3[i]]]
}
rm(int1, int2, int3)

subset_l <- c(p1_l, p2_l, p3_l)
subset <- c(p1, p2, p3)
subset_u <- c(p1_u, p2_u, p3_u)
#####

#aggregate median irf for plot
#####
irfdf1_l <- matrix(nrow = 20, ncol = length(subset))
irfdf1 <- matrix(nrow = 20, ncol = length(subset))
irfdf1_u <- matrix(nrow = 20, ncol = length(subset))

irfdf2_l <- matrix(nrow = 20, ncol = length(subset))
irfdf2 <- matrix(nrow = 20, ncol = length(subset))
irfdf2_u <- matrix(nrow = 20, ncol = length(subset))

irfdf3_l <- matrix(nrow = 20, ncol = length(subset))
irfdf3 <- matrix(nrow = 20, ncol = length(subset))
irfdf3_u <- matrix(nrow = 20, ncol = length(subset))



for (j in 1:20){
  for(i in 1:length(subset)){
    irfdf1_l[j , i] <- subset_l[[i]][j,1]
    irfdf1[j , i] <- subset[[i]][j,1]
    irfdf1_u[j , i] <- subset_u[[i]][j,1]
  }
}

for (j in 1:20){
  for(i in 1:length(subset)){
    irfdf2_l[j , i] <- subset_l[[i]][j,2]
    irfdf2[j , i] <- subset[[i]][j,2]
    irfdf2_u[j , i] <- subset_u[[i]][j,2]
  }
}

for (j in 1:20){
  for(i in 1:length(subset)){
    irfdf3_l[j , i] <- subset_l[[i]][j,3]
    irfdf3[j , i] <- subset[[i]][j,3]
    irfdf3_u[j , i] <- subset_u[[i]][j,3]
  }
}

irfdata_l <- matrix(nrow = 20, ncol = 3)
irfdata <- matrix(nrow = 20, ncol = 3)
irfdata_u <- matrix(nrow = 20, ncol = 3)
for (i in 1:20){
  irfdata_l[i, 1] <- median(irfdf1_l[i,1:ncol(irfdf1)])
  irfdata_l[i, 2] <- median(irfdf2_l[i,1:ncol(irfdf2)])
  irfdata_l[i, 3] <- median(irfdf3_l[i,1:ncol(irfdf3)])
  
  irfdata[i, 1] <- median(irfdf1[i,1:ncol(irfdf1)])
  irfdata[i, 2] <- median(irfdf2[i,1:ncol(irfdf2)])
  irfdata[i, 3] <- median(irfdf3[i,1:ncol(irfdf3)])
  
  irfdata_u[i, 1] <- median(irfdf1_u[i,1:ncol(irfdf1)])
  irfdata_u[i, 2] <- median(irfdf2_u[i,1:ncol(irfdf2)])
  irfdata_u[i, 3] <- median(irfdf3_u[i,1:ncol(irfdf3)])
}
irfdata_l <- irfdata_l %>% data.frame() %>%
  rename(FFR_l = X1,
         CPI_l = X2,
         INDPRO_l = X3)
irfdata <- irfdata %>% data.frame() %>%
  rename(FFR = X1,
         CPI = X2,
         INDPRO = X3)
irfdata_u <- irfdata_u %>% data.frame() %>%
  rename(FFR_u = X1,
         CPI_u = X2,
         INDPRO_u = X3)

irfdata <- cbind(irfdata, irfdata_l, irfdata_u)
rm(irfdata_l, irfdata_u)

rm(subset_l, irfdf1_l, irfdf2_l, irfdf3_l,
   subset, irfdf1, irfdf2, irfdf3,
   subset_u, irfdf1_u, irfdf2_u, irfdf3_u)


wholeperiod_estim <- model2
rm(model1, model2)

IRFS_total_l <- matrix(nrow = 20, ncol = 3)
IRFS_total <- matrix(nrow = 20, ncol = 3)
IRFS_total_u <- matrix(nrow = 20, ncol = 3)

for (j in 1:3){
  for (i in 1:20){
    IRFS_total_l[i ,j] <- quantile(wholeperiod_estim$IRFS[1:200, i ,j], probs = 0.14)
    IRFS_total[i ,j] <- median(wholeperiod_estim$IRFS[1:200, i ,j])
    IRFS_total_u[i ,j] <- quantile(wholeperiod_estim$IRFS[1:200, i ,j], probs = 0.86)
  }
}

IRFS_total_l <- IRFS_total_l %>% data.frame()%>%
  rename(FFR_l = X1,
         CPI_l = X2,
         INDPRO_l = X3)
IRFS_total <- IRFS_total %>% data.frame()%>%
  rename(FFR = X1,
         CPI = X2,
         INDPRO = X3)
IRFS_total_u <- IRFS_total_u %>% data.frame()%>%
  rename(FFR_u = X1,
         CPI_u = X2,
         INDPRO_u = X3)
IRFS_total <- cbind(IRFS_total, IRFS_total_l, IRFS_total_u)
rm(IRFS_total_l, IRFS_total_u)
#####

#create plots
#####
#plots for the non-rw estimation
ffr1 <- ggplot(IRFS_total, aes(x = c(seq(from = 1, to = 20)), y  = FFR, ymin = FFR_l, ymax = FFR_u)) +
  geom_line() +
  geom_hline(yintercept = 0, color = "red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  theme_light() +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11)) +
#  ggtitle("Orthogonal impulse response from FFR shock") +
  xlab("Month") +
  ylab("FFR")

cpi1 <- ggplot(IRFS_total, aes(x = c(seq(from = 1, to = 20)), y  = CPI, ymin = CPI_l, ymax = CPI_u)) +
  geom_line() +
  geom_hline(yintercept = 0, color = "red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  theme_light() +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11)) +
#  ggtitle("Orthogonal impulse response from FFR shock") +
  xlab("Month") +
  ylab("CPI")

indpro1 <- ggplot(IRFS_total, aes(x = c(seq(from = 1, to = 20)), y  = INDPRO, ymin = INDPRO_l, ymax = INDPRO_u)) +
  geom_line() +
  geom_hline(yintercept = 0, color = "red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  theme_light() +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11)) +
#  ggtitle("Orthogonal impulse response from FFR shock") +
  xlab("Month") +
  ylab("INDPRO")

#plots from rw estimation
ffr2 <- ggplot(irfdata, aes(x = c(seq(from = 1, to = 20)), y  = FFR, ymin = FFR_l, ymax = FFR_u)) +
  geom_line() +
  geom_hline(yintercept = 0, color = "red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  theme_light() +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11)) +
#  ggtitle("Orthogonal impulse response from FFR shock") +
  xlab("Month") +
  ylab("FFR")


cpi2 <- ggplot(irfdata, aes(x = c(seq(from = 1, to = 20)), y  = CPI, ymin = CPI_l, ymax = CPI_u)) +
  geom_line() +
  geom_hline(yintercept = 0, color = "red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  theme_light() +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11)) +
 # ggtitle("Orthogonal impulse response from FFR shock") +
  xlab("Month") +
  ylab("CPI")


indpro2 <- ggplot(irfdata, aes(x = c(seq(from = 1, to = 20)), y  = INDPRO, ymin = INDPRO_l, ymax = INDPRO_u)) +
  geom_line() +
  geom_hline(yintercept = 0, color = "red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  theme_light() +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11)) +
 # ggtitle("Orthogonal impulse response from FFR shock") +
  xlab("Month") +
  ylab("INDPRO")
#####

#combine plots

x1 <- grid.arrange(ffr1, cpi1, indpro1, ncol = 3, top = "Orthogonal impulse responses from FFR shock, baseline model")
x2 <- grid.arrange(ffr2, cpi2, indpro2, ncol = 3, top = "Orthogonal impulse responses from FFR shock, rolling window")

irfplot1 <- grid.arrange(x1, x2, ncol = 1)
ggsave(paste(folder,"irfplot1.png",sep=""), irfplot1, bg = "transparent")


#decompose high uncertainty periods with irfplots
#####
#p1
irfdf1_l <- matrix(nrow = 20, ncol = length(p1))
irfdf1 <- matrix(nrow = 20, ncol = length(p1))
irfdf1_u <- matrix(nrow = 20, ncol = length(p1))

irfdf2_l <- matrix(nrow = 20, ncol = length(p1))
irfdf2 <- matrix(nrow = 20, ncol = length(p1))
irfdf2_u <- matrix(nrow = 20, ncol = length(p1))

irfdf3_l <- matrix(nrow = 20, ncol = length(p1))
irfdf3 <- matrix(nrow = 20, ncol = length(p1))
irfdf3_u <- matrix(nrow = 20, ncol = length(p1))



for (j in 1:20){
  for(i in 1:length(p1)){
    irfdf1_l[j , i] <- p1_l[[i]][j,1]
    irfdf1[j , i] <- p1[[i]][j,1]
    irfdf1_u[j , i] <- p1_u[[i]][j,1]
  }
}

for (j in 1:20){
  for(i in 1:length(p1)){
    irfdf2_l[j , i] <- p1_l[[i]][j,2]
    irfdf2[j , i] <- p1[[i]][j,2]
    irfdf2_u[j , i] <- p1_u[[i]][j,2]
  }
}

for (j in 1:20){
  for(i in 1:length(p1)){
    irfdf3_l[j , i] <- p1_l[[i]][j,3]
    irfdf3[j , i] <- p1[[i]][j,3]
    irfdf3_u[j , i] <- p1_u[[i]][j,3]
  }
}

irfdata_l <- matrix(nrow = 20, ncol = 3)
irfdata <- matrix(nrow = 20, ncol = 3)
irfdata_u <- matrix(nrow = 20, ncol = 3)
for (i in 1:20){
  irfdata_l[i, 1] <- median(irfdf1_l[i,1:ncol(irfdf1)])
  irfdata_l[i, 2] <- median(irfdf2_l[i,1:ncol(irfdf2)])
  irfdata_l[i, 3] <- median(irfdf3_l[i,1:ncol(irfdf3)])
  
  irfdata[i, 1] <- median(irfdf1[i,1:ncol(irfdf1)])
  irfdata[i, 2] <- median(irfdf2[i,1:ncol(irfdf2)])
  irfdata[i, 3] <- median(irfdf3[i,1:ncol(irfdf3)])
  
  irfdata_u[i, 1] <- median(irfdf1_u[i,1:ncol(irfdf1)])
  irfdata_u[i, 2] <- median(irfdf2_u[i,1:ncol(irfdf2)])
  irfdata_u[i, 3] <- median(irfdf3_u[i,1:ncol(irfdf3)])
}
irfdata_l <- irfdata_l %>% data.frame() %>%
  rename(FFR_l = X1,
         CPI_l = X2,
         INDPRO_l = X3)
irfdata <- irfdata %>% data.frame() %>%
  rename(FFR = X1,
         CPI = X2,
         INDPRO = X3)
irfdata_u <- irfdata_u %>% data.frame() %>%
  rename(FFR_u = X1,
         CPI_u = X2,
         INDPRO_u = X3)

irfdata <- cbind(irfdata, irfdata_l, irfdata_u)
rm(irfdata_l, irfdata_u)

rm(p1_l, irfdf1_l, irfdf2_l, irfdf3_l,
   p1, irfdf1, irfdf2, irfdf3,
   p1_u, irfdf1_u, irfdf2_u, irfdf3_u)

ffr3 <- ggplot(irfdata, aes(x = c(seq(from = 1, to = 20)), y  = FFR, ymin = FFR_l, ymax = FFR_u)) +
  geom_line() +
  geom_hline(yintercept = 0, color = "red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  theme_light() +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11)) +
 # ggtitle("Orthogonal impulse response from FFR shock") +
  xlab("Month") +
  ylab("FFR")


cpi3 <- ggplot(irfdata, aes(x = c(seq(from = 1, to = 20)), y  = CPI, ymin = CPI_l, ymax = CPI_u)) +
  geom_line() +
  geom_hline(yintercept = 0, color = "red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  theme_light() +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11)) +
#  ggtitle("Orthogonal impulse response from FFR shock") +
  xlab("Month") +
  ylab("CPI")


indpro3 <- ggplot(irfdata, aes(x = c(seq(from = 1, to = 20)), y  = INDPRO, ymin = INDPRO_l, ymax = INDPRO_u)) +
  geom_line() +
  geom_hline(yintercept = 0, color = "red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  theme_light() +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11)) +
#  ggtitle("Orthogonal impulse response from FFR shock") +
  xlab("Month") +
  ylab("INDPRO")

#p2
irfdf1_l <- matrix(nrow = 20, ncol = length(p2))
irfdf1 <- matrix(nrow = 20, ncol = length(p2))
irfdf1_u <- matrix(nrow = 20, ncol = length(p2))

irfdf2_l <- matrix(nrow = 20, ncol = length(p2))
irfdf2 <- matrix(nrow = 20, ncol = length(p2))
irfdf2_u <- matrix(nrow = 20, ncol = length(p2))

irfdf3_l <- matrix(nrow = 20, ncol = length(p2))
irfdf3 <- matrix(nrow = 20, ncol = length(p2))
irfdf3_u <- matrix(nrow = 20, ncol = length(p2))



for (j in 1:20){
  for(i in 1:length(p2)){
    irfdf1_l[j , i] <- p2_l[[i]][j,1]
    irfdf1[j , i] <- p2[[i]][j,1]
    irfdf1_u[j , i] <- p2_u[[i]][j,1]
  }
}

for (j in 1:20){
  for(i in 1:length(p2)){
    irfdf2_l[j , i] <- p2_l[[i]][j,2]
    irfdf2[j , i] <- p2[[i]][j,2]
    irfdf2_u[j , i] <- p2_u[[i]][j,2]
  }
}

for (j in 1:20){
  for(i in 1:length(p2)){
    irfdf3_l[j , i] <- p2_l[[i]][j,3]
    irfdf3[j , i] <- p2[[i]][j,3]
    irfdf3_u[j , i] <- p2_u[[i]][j,3]
  }
}

irfdata_l <- matrix(nrow = 20, ncol = 3)
irfdata <- matrix(nrow = 20, ncol = 3)
irfdata_u <- matrix(nrow = 20, ncol = 3)
for (i in 1:20){
  irfdata_l[i, 1] <- median(irfdf1_l[i,1:ncol(irfdf1)])
  irfdata_l[i, 2] <- median(irfdf2_l[i,1:ncol(irfdf2)])
  irfdata_l[i, 3] <- median(irfdf3_l[i,1:ncol(irfdf3)])
  
  irfdata[i, 1] <- median(irfdf1[i,1:ncol(irfdf1)])
  irfdata[i, 2] <- median(irfdf2[i,1:ncol(irfdf2)])
  irfdata[i, 3] <- median(irfdf3[i,1:ncol(irfdf3)])
  
  irfdata_u[i, 1] <- median(irfdf1_u[i,1:ncol(irfdf1)])
  irfdata_u[i, 2] <- median(irfdf2_u[i,1:ncol(irfdf2)])
  irfdata_u[i, 3] <- median(irfdf3_u[i,1:ncol(irfdf3)])
}
irfdata_l <- irfdata_l %>% data.frame() %>%
  rename(FFR_l = X1,
         CPI_l = X2,
         INDPRO_l = X3)
irfdata <- irfdata %>% data.frame() %>%
  rename(FFR = X1,
         CPI = X2,
         INDPRO = X3)
irfdata_u <- irfdata_u %>% data.frame() %>%
  rename(FFR_u = X1,
         CPI_u = X2,
         INDPRO_u = X3)

irfdata <- cbind(irfdata, irfdata_l, irfdata_u)
rm(irfdata_l, irfdata_u)

rm(p2_l, irfdf1_l, irfdf2_l, irfdf3_l,
   p2, irfdf1, irfdf2, irfdf3,
   p2_u, irfdf1_u, irfdf2_u, irfdf3_u)

ffr4 <- ggplot(irfdata, aes(x = c(seq(from = 1, to = 20)), y  = FFR, ymin = FFR_l, ymax = FFR_u)) +
  geom_line() +
  geom_hline(yintercept = 0, color = "red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  theme_light() +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11)) +
#  ggtitle("Orthogonal impulse response from FFR shock") +
  xlab("Month") +
  ylab("FFR")


cpi4 <- ggplot(irfdata, aes(x = c(seq(from = 1, to = 20)), y  = CPI, ymin = CPI_l, ymax = CPI_u)) +
  geom_line() +
  geom_hline(yintercept = 0, color = "red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  theme_light() +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11)) +
 # ggtitle("Orthogonal impulse response from FFR shock") +
  xlab("Month") +
  ylab("CPI")


indpro4 <- ggplot(irfdata, aes(x = c(seq(from = 1, to = 20)), y  = INDPRO, ymin = INDPRO_l, ymax = INDPRO_u)) +
  geom_line() +
  geom_hline(yintercept = 0, color = "red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  theme_light() +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11)) +
 # ggtitle("Orthogonal impulse response from FFR shock") +
  xlab("Month") +
  ylab("INDPRO")


#p3
irfdf1_l <- matrix(nrow = 20, ncol = length(p3))
irfdf1 <- matrix(nrow = 20, ncol = length(p3))
irfdf1_u <- matrix(nrow = 20, ncol = length(p3))

irfdf2_l <- matrix(nrow = 20, ncol = length(p3))
irfdf2 <- matrix(nrow = 20, ncol = length(p3))
irfdf2_u <- matrix(nrow = 20, ncol = length(p3))

irfdf3_l <- matrix(nrow = 20, ncol = length(p3))
irfdf3 <- matrix(nrow = 20, ncol = length(p3))
irfdf3_u <- matrix(nrow = 20, ncol = length(p3))



for (j in 1:20){
  for(i in 1:length(p3)){
    irfdf1_l[j , i] <- p3_l[[i]][j,1]
    irfdf1[j , i] <- p3[[i]][j,1]
    irfdf1_u[j , i] <- p3_u[[i]][j,1]
  }
}

for (j in 1:20){
  for(i in 1:length(p3)){
    irfdf2_l[j , i] <- p3_l[[i]][j,2]
    irfdf2[j , i] <- p3[[i]][j,2]
    irfdf2_u[j , i] <- p3_u[[i]][j,2]
  }
}

for (j in 1:20){
  for(i in 1:length(p3)){
    irfdf3_l[j , i] <- p3_l[[i]][j,3]
    irfdf3[j , i] <- p3[[i]][j,3]
    irfdf3_u[j , i] <- p3_u[[i]][j,3]
  }
}

irfdata_l <- matrix(nrow = 20, ncol = 3)
irfdata <- matrix(nrow = 20, ncol = 3)
irfdata_u <- matrix(nrow = 20, ncol = 3)
for (i in 1:20){
  irfdata_l[i, 1] <- median(irfdf1_l[i,1:ncol(irfdf1)])
  irfdata_l[i, 2] <- median(irfdf2_l[i,1:ncol(irfdf2)])
  irfdata_l[i, 3] <- median(irfdf3_l[i,1:ncol(irfdf3)])
  
  irfdata[i, 1] <- median(irfdf1[i,1:ncol(irfdf1)])
  irfdata[i, 2] <- median(irfdf2[i,1:ncol(irfdf2)])
  irfdata[i, 3] <- median(irfdf3[i,1:ncol(irfdf3)])
  
  irfdata_u[i, 1] <- median(irfdf1_u[i,1:ncol(irfdf1)])
  irfdata_u[i, 2] <- median(irfdf2_u[i,1:ncol(irfdf2)])
  irfdata_u[i, 3] <- median(irfdf3_u[i,1:ncol(irfdf3)])
}
irfdata_l <- irfdata_l %>% data.frame() %>%
  rename(FFR_l = X1,
         CPI_l = X2,
         INDPRO_l = X3)
irfdata <- irfdata %>% data.frame() %>%
  rename(FFR = X1,
         CPI = X2,
         INDPRO = X3)
irfdata_u <- irfdata_u %>% data.frame() %>%
  rename(FFR_u = X1,
         CPI_u = X2,
         INDPRO_u = X3)

irfdata <- cbind(irfdata, irfdata_l, irfdata_u)
rm(irfdata_l, irfdata_u)

rm(p3_l, irfdf1_l, irfdf2_l, irfdf3_l,
   p3, irfdf1, irfdf2, irfdf3,
   p3_u, irfdf1_u, irfdf2_u, irfdf3_u)

ffr5 <- ggplot(irfdata, aes(x = c(seq(from = 1, to = 20)), y  = FFR, ymin = FFR_l, ymax = FFR_u)) +
  geom_line() +
  geom_hline(yintercept = 0, color = "red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  theme_light() +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11)) +
 # ggtitle("Orthogonal impulse response from FFR shock") +
  xlab("Month") +
  ylab("FFR")


cpi5 <- ggplot(irfdata, aes(x = c(seq(from = 1, to = 20)), y  = CPI, ymin = CPI_l, ymax = CPI_u)) +
  geom_line() +
  geom_hline(yintercept = 0, color = "red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  theme_light() +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11)) +
 # ggtitle("Orthogonal impulse response from FFR shock") +
  xlab("Month") +
  ylab("CPI")


indpro5 <- ggplot(irfdata, aes(x = c(seq(from = 1, to = 20)), y  = INDPRO, ymin = INDPRO_l, ymax = INDPRO_u)) +
  geom_line() +
  geom_hline(yintercept = 0, color = "red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  theme_light() +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11)) +
#  ggtitle("Orthogonal impulse response from FFR shock") +
  xlab("Month") +
  ylab("INDPRO")
#####


#combine 
y1 <- grid.arrange(ffr3, cpi3, indpro3, ncol = 3, top = "Orthogonal impulse responses from FFR shock, Jan. 2007 - Jan. 2015")
y2 <- grid.arrange(ffr4, cpi4, indpro4, ncol = 3, top = "Orthogonal impulse responses from FFR shock, Aug. 1999 - Jan. 2005")
y3 <- grid.arrange(ffr5, cpi5, indpro5, ncol = 3, top = "Orthogonal impulse responses from FFR shock, Jun. 1989 - Jan. 1995")

irfplot2 <- grid.arrange(y1, y2, y3, ncol = 1) 
ggsave(paste(folder,"irfplot2.png",sep=""), irfplot2, bg = "transparent")




b