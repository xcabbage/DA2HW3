---
title: "DA2HA3"
output: html_document
---


```{r}
library(readr)
library(dplyr)
library(dummies)
library(ggplot2)
library(lmtest)  # for Breush-Godfrey
library(stargazer)
library(pander)
library(stats)
library(data.table)
library(readxl)
library(vars)
library(zoo)  # easier work with dates
library(gridExtra)  # to combine graphs
library(urca)  # for unit root and cointegration tests
library(sandwich)  # for robust estimation of covariance matrix


source('C:/Users/KaposztassyG/Documents/_CEU/DA2/HW3/da_helper_functions.R')

turkey2 <- read_excel('C:/Users/KaposztassyG/Documents/_CEU/DA2/HW3/turkey2.xlsx')

turkey2$date <-as.Date(turkey2$date,"%y-%m-%d")

turkey2 <- na.omit(turkey2)

setDT(turkey2)

str(turkey2)

turkey2 <- turkey2 %>% mutate(
  month = as.numeric(format(date, '%m')),
  ln_er = log(exch_rate),
  ln_ind = log(tot_ind),
  ln_tot_ind = log(tot_ind),
  ln_exch_rate = log(exch_rate),
  dpct_er = 100*(ln_er - lag(ln_er)),
  dpct_ind = 100*(ln_ind - lag(ln_ind)),
  ld_er = d(log(exch_rate)),
  ld_ind = d(log(tot_ind))
)

rm(list=ls())
turkey <- read.csv("/Users/Andi/Desktop/turkey2.csv")

#Data transformation
View(turkey)
names(turkey)
class(turkey$date)


#Exploratory plots
turkey2 %>% ggplot(aes(x = date, y = tot_ind)) + geom_line(size = 0.5, col='red') + theme_minimal()
turkey2 %>% ggplot(aes(x = date, y = exch_rate)) + geom_line(size = 0.5, col='green') + theme_minimal()
turkey2 %>% ggplot(aes(x = date, y = ln_tot_ind)) + geom_line(size = 0.5, col='darkred') + theme_minimal()
turkey2 %>% ggplot(aes(x = date, y = ln_exch_rate)) + geom_line(size = 0.5, col='darkgreen') + theme_minimal()

# Phillips-Perron tests
variables_to_test <- c("date", "tot_ind", "exch_rate", "ln_tot_ind", "ln_exch_rate")

for (variable in variables_to_test) {
  print(variable)
  pperron(turkey2[[variable]])
}

myreg1 <- lm(tot_ind ~ exch_rate, turkey2[date > "2001-03-01"])
summary(myreg1)

turkey3 <- turkey2 %>% 
  mutate(
  ld_tot_ind = d(log(tot_ind)),
  ld_exch_rate = d(log(exch_rate))
)
str(turkey3)

myreg2 <- lm(ld_tot_ind ~ ld_exch_rate, turkey3[date > "2001-03-01"])
summary(myreg2)

#Task 2 - Do a VAR with Y and X

Irf <- function(...) {  # to calculate simple IRF instead of orthogonal by default
  irf(..., ortho = FALSE)
}

turkey3 <- as.data.frame(turkey3)
turkey3 <- na.omit(turkey3)
  
exch_rate_var <- data.frame(cbind(turkey3$ld_tot_ind, turkey3$ld_exch_rate))
ind_outp_var <- data.frame(cbind(turkey3$ln_exch_rate, turkey3$ln_tot_ind))

colnames(exch_rate_var) <- c("ld_tot_ind","ld_exch_rate")
colnames(ind_outp_var) <- c("ln_exch_rate","ln_tot_ind")


var1_ger <- VAR(exch_rate_var, p = 1)
summary(var1_ger)
plot(Irf(var1_ger, impulse = 'ld_tot_ind', response = 'ld_exch_rate'))

var2_ger <- VAR(exch_rate_var, p = 2)
summary(var2_ger)
plot(Irf(var2_ger, impulse = 'ld_tot_ind', response = 'ld_exch_rate'))

var3_ger <- VAR(exch_rate_var, p = 3)
summary(var3_ger)
plot(Irf(var3_ger, impulse = 'ld_tot_ind', response = 'ld_exch_rate'))

var5_ger <- VAR(exch_rate_var, p = 5)
summary(var5_ger)
plot(Irf(var5_ger, impulse = 'ld_tot_ind', response = 'ld_exch_rate'))

#VAR with exchange rate

var_er <- VAR(ind_outp_var, p = 3)
summary(var_er)
plot(Irf(var_er, impulse = 'ln_exch_rate', response = 'ln_tot_ind'))


#Granger causality
causality(var_er, cause = 'ln_tot_ind') # reject that no Granger-causality
causality(var_er, cause = 'ln_exch_rate')  # reject that no Granger-causality

causality(var5_ger, cause = 'ld_tot_ind') # reject that no Granger-causality
causality(var5_ger, cause = 'ld_exch_rate')  # reject that no Granger-causality




# table1

myreg1 <- lm(d(dpct_er) ~ 1, turkey3)
summary_r(myreg1, se = 'newey-west', max_lag = 25)

myreg2 <- lm(d(dpct_ind) ~ 1, turkey3)
summary_r(myreg2, se = 'newey-west', max_lag = 25)

myreg3 <- lm(dpct_er ~ factor(month), turkey3)
summary_r(myreg3, se = 'newey-west', max_lag = 25)

myreg4 <- lm(dpct_ind ~ factor(month), turkey3)
summary_r(myreg4, se = 'newey-west', max_lag = 25)

stargazer_r(
  list(myreg1, myreg2, myreg3, myreg4), 
  se = 'newey-west', max_lag = 25, digits = 2,
  omit.stat=c("LL", "aic", "ser", "f", "adj.rsq", "rsq")
)


# table2

myreg5 <- lm(dpct_ind ~ dpct_er, turkey3)
myreg6 <- lm(dpct_ind ~ dpct_er + factor(month), turkey3)

stargazer_r(
  list(myreg5, myreg5, myreg5, myreg5, myreg6), 
  se = list('robust', 'newey-west', 'newey-west', 'newey-west', 'newey-west'),
  max_lag = list(NA, 2, 13, 25, 25),
  digits = 2,
  omit.stat=c("LL", "aic", "ser", "f", "adj.rsq", "rsq"),
  dep.var.caption = 'LHS: % change in price',
  dep.var.labels.include = FALSE,
  column.labels = c('robust SE', 'NW SE, lags 2', 'NW SE, lags 13', 'NW SE, lags 25', 'NW SE, lags 25')
)

# table3

myreg7 <- lm(dpct_ind ~ dpct_er + lag(dpct_er), turkey3)
# with automatic lag inclusion
myreg7 <- lm(
  as.formula(
    paste("dpct_ind ~ ", lags(dpct_er, 0:1))
  ), 
  turkey3
)

# everything together
max_lag <- list(1, 4, 7, 13)
myregs <- lapply(
  max_lag, function(l) {
    lm(
      as.formula(paste("dpct_ind ~", lags(dpct_er, 0:l))),
      turkey3
    )
  }
)
myreg8 <- lm(
  as.formula(
    paste("dpct_ind ~", lags(dpct_er, 0:13), "+ factor(month)")
  ),
  turkey3    
)
myregs[[5]] <- myreg8

stargazer_r(
  myregs, 
  se = 'newey-west',
  max_lag = 25,
  digits = 2,
  omit.stat=c("LL", "aic", "ser", "f", "adj.rsq", "rsq"),
  dep.var.caption = 'LHS: % change in ind.output',
  dep.var.labels.include = FALSE,
  column.labels = c(sapply(max_lag, function(l) paste0("L(", l, ")")), "L(13)")
)

# table4

myreg9 <- lm(
  as.formula(
    paste("dpct_ind ~", lags(dpct_er, c(1, 4, 7, 13)))
  ),
  turkey3
)
myreg9

der_lag <- c(1, 4, 7, 13)
dder_max_lag <- c(0, 3, 6, 12)
myregs2 <- lapply(
  seq_along(der_lag), function(i) {
    lm(
      as.formula(
        paste(
          "dpct_ind ~", 
          lags(dpct_er, der_lag[i]), "+", lags(d(dpct_er), 0:dder_max_lag[i])
        )
      ),
      turkey3
    )
  }
)
myregs2 <- c(list(myreg2), myregs2)
stargazer_r(
  myregs2, 
  se = 'newey-west',
  max_lag = 25,
  digits = 2,
  omit.stat=c("LL", "aic", "ser", "f", "adj.rsq", "rsq"),
  dep.var.caption = 'LHS: % change in ind.output',
  dep.var.labels.include = FALSE,
  column.labels = c(sapply(der_lag, function(l) paste0("L(", l, ")")))
)
```
