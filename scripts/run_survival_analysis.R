#### Corrected-paths copy of DGRP_survival_analysis.R for local reproduction
#### Original: Dr Rita Ibrahim & Dr Adam Dobson
#### Paths fixed via here(); tSNE section skipped (needs XGBoost WNS output not yet generated)
#### Do NOT modify the original script.

rm(list=ls())

library(survival)
library(coxme)
library(car)
library(dplyr)
library(emmeans)
library(ggplot2)
library(effectsize)
library(ggh4x)
library(knitr)
library(ciTools)
library(here)
library(survminer)
library(rms)
library(see)

jmp2r <- function(x){
  for(i in 1:nrow(x)){
    xx <- x[i,]
    if(xx$event[1] > 0){
      xxx <- matrix(rep(as.matrix(xx), xx[1,"event"]), ncol=ncol(xx), byrow=T)
      if(exists("a")){a <- rbind(a, xxx)}else{(a <- xxx)}
    }
  }
  a <- as.data.frame(a)
  for(i in 1:ncol(a)){
    colnames(a)[i] <- colnames(x)[i]
    if(class(a[,i]) != class(x[,i])) {a[,i] <- as.character(a[,i]); class(a[,i]) <- class(x[,i])}
  }
  a <- a[,!colnames(a) %in% c("count", "deathIsZero_censorIsOne", "event")]
  a
}

d_0 <- read.csv(here("Survival-data", "DGRP-starvationresistance.csv"))
d <- jmp2r(d_0)
d$DGRP <- as.factor(d$DGRP)

dd <- datadist(d)
options(datadist = 'dd')

psm1 <- psm(Surv(time, censor) ~ DGRP,
            data = d,
            dist = "logistic")

psm1_emm <- emmeans(psm1, ~ DGRP)
emmeanss <- as.data.frame(psm1_emm)
colnames(emmeanss) <- c("DGRP", "emmean", "SE", "df", "CI05", "CI95")
emmeanss <- emmeanss[order(emmeanss$emmean), ]

emmeans_df <- as.data.frame(emmeanss)
emmeans_df$DGRP <- gsub("-", "", emmeans_df$DGRP)

write.csv(emmeans_df, file = here("Emmeans.csv"), row.names = FALSE)

lower_threshold <- quantile(emmeans_df$emmean, 0.20)
upper_threshold <- quantile(emmeans_df$emmean, 0.80)

sensitive <- c()
resistant  <- c()
for (i in 1:nrow(emmeans_df)) {
  if (emmeans_df$emmean[i] < lower_threshold) {
    sensitive <- append(sensitive, emmeans_df[i, ]$DGRP)
  } else if (emmeans_df$emmean[i] > upper_threshold) {
    resistant <- append(resistant, emmeans_df[i, ]$DGRP)
  }
}

cat("Sensitive lines:", length(sensitive), "\n")
cat("Resistant lines:", length(resistant), "\n")

sensitive_df <- as.data.frame(sensitive)
sensitive_df$sensitive <- gsub("-", "", sensitive_df$sensitive)
resistant_df <- as.data.frame(resistant)
resistant_df$resistant <- gsub("-", "", resistant_df$resistant)

write.csv(sensitive_df, file = here("sensitive_df_20pct_emmean.csv"), row.names = FALSE)
write.csv(resistant_df, file = here("resistant_df_80pct_emmean.csv"), row.names = FALSE)

cat("Done. Files written to", here(), "\n")
