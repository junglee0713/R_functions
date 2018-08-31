###=====
### Collection of R utility functions
###=====

library(tidyverse)
library(reshape2)

###=====
###  pval lmer: LRT based LMM
###=====

pval_lmer <- function(formula, data) {
  
  require(lme4)
  require(dplyr)
  require(reshape2)
  
  formula <- as.formula(formula)
  parsed_formula <- reshape2:::parse_formula(formula)  ### parse the formula
  dep_var <- as.character(parsed_formula[[1]])   ### dependent variable
  pred_vars <- as.character(parsed_formula[[2]]) ### predictor variables
  random_effect_vars <- pred_vars[grepl("(", pred_vars, fixed=TRUE)]
  fixed_effect_vars <- pred_vars[!(pred_vars %in% random_effect_vars)]
  
  ### full model for likelihood comparison
  fit_full <- lmer(formula, data, REML=F)
  
  anova_table_list <- list(NULL)
  
  for (i in 1:length(fixed_effect_vars)) {
    curr_fixed_effect_var <- fixed_effect_vars[i]
    reduced_pred_vars <- setdiff(pred_vars, curr_fixed_effect_var)
    reduced_formula <- as.formula(paste(c(dep_var,  "~", paste(reduced_pred_vars, collapse=" + ")), collapse=" "))
    fit_reduced <- lmer(reduced_formula, data, REML = F)
    anova_table_list[[i]] <- data.frame(Term=curr_fixed_effect_var, 
                                        anova(fit_reduced, fit_full)["fit_full", c("Chisq", "Chi Df", "Pr(>Chisq)")],
                                        stringsAsFactors = F)
  }
  
  anova_table <- bind_rows(anova_table_list) %>%
    rename(`$p$-value`=Pr..Chisq.)
  return(anova_table)
}

###=====
###  is_outlier: detect outlier
###=====

is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}

###=====
###  tpairwise: tidy result from pairwise t-test
###=====

tpairwise <- function(data, response, group) {
  df <- data %>%
    select_(response, group) %>%
    setNames(c("y", "gr")) %>%
    mutate(gr = as.factor(gr))
  res <- pairwise.t.test(df$y, df$gr, p.adjust.method = "none")$p.value %>%
    as.data.frame()
  comp <- character()
  pval <- double()
  for (i in 1:ncol(res)) {
    gr1 <- colnames(res)[i] 
    for (j in i:nrow(res)) {
      gr2 <- rownames(res)[j]
      comp <- c(comp, paste0(gr1, " vs ", gr2))
      pval <- c(pval, res[gr2, gr1])
    }
  }
  out <- data.frame(comp = comp, pval = pval, stringsAsFactors = F)
  return(out)
}


