#!/usr/bin/env Rscript
# Debug test to understand column structure

library(BUGSnet)
library(dplyr)

data(thrombolytic)
thrombo.slr <- data.prep(arm.data = thrombolytic,
                         varname.t = "treatment",
                         varname.s = "study")

thrombo.model <- nma.model(data = thrombo.slr,
                           outcome = "events",
                           N = "sampleSize",
                           reference = "SK",
                           family = "binomial",
                           link = "log",
                           effects = "random",
                           type = "consistency")

jags_result <- nma.run(model = thrombo.model,
                       n.adapt = 500,
                       n.burnin = 500,
                       n.iter = 1000)

cat("trt.key:", paste(jags_result$trt.key, collapse = ", "), "\n\n")

samples_matrix <- do.call(rbind, jags_result$samples) %>% data.frame()
cat("Column names in samples_matrix:\n")
cat(paste(colnames(samples_matrix)[1:15], collapse = "\n"), "\n\n")

cat("Trying to access d.tPA:\n")
print(head(samples_matrix[[paste0("d.", "tPA")]]))

cat("\nTrying with backticks:\n")
print(head(samples_matrix[[paste0("`d.", "tPA", "`")]]))

