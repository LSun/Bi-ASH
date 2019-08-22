#############
## Plot the results from misspecified_sims.R
#############

suppressPackageStartupMessages(library(tidyverse))
library(broom)
library(ggthemes)

## Read in and clean data ------------------------------------------------------
ncmat <- readRDS("~/GitHub/Bi-ASH/output/nc_sims.rds")

subnamevec <- colnames(ncmat)[-(1:10)]
# num_method <- length(subnamevec) / 3
num_method <- length(subnamevec) / 2
subnamevec[1:num_method] <- str_c("MSE_", subnamevec[1:num_method])
subnamevec[(num_method + 1):(2 * num_method)] <- str_c("AUC_", subnamevec[(num_method + 1):(2 * num_method)])
# subnamevec[(2*num_method + 1):(3*num_method)] <- str_c("cov_", subnamevec[(2*num_method + 1):(3*num_method)])

colnames(ncmat)[-(1:10)] <- subnamevec

ncmatdf <- as_tibble(ncmat)

ncmatdf %>%
  gather(-log2foldsd,
         -Ngene,
         -log2foldmean,
         -skip_gene,
         -current_seed,
         -nullpi,
         -Nsamp,
         -ncontrols,
         -prop_cont,
         -poisthin,
         key = "metric_method",
         value = "value") %>%
  select(nullpi, metric_method, value) %>%
  separate(metric_method, into = c("metric", "method"), sep = "_") ->
  ncmatdf_long

ncmatdf_long$nullpi <- paste0("pi[0] == ", ncmatdf_long$nullpi)
ncmatdf_long$metric <- factor(ncmatdf_long$metric, levels = c("MSE", "AUC"))
ncmatdf_long$method <- factor(ncmatdf_long$method, levels = c("ols", "ruv2", "ruv3", "cate", "ashr", "biashr"))

plot.nc.sim <- ggplot(ncmatdf_long[ncmatdf_long$method != "ashr", ], aes(x = method, y = value, group = method)) +
  geom_boxplot(aes(fill = method)) +
  facet_grid(metric ~ nullpi, scale = "free_y", labeller = label_parsed) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 15),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 15),
    strip.text.x = element_text(size = 15),
    strip.text.y = element_text(size = 15)
  )

ggsave("~/GitHub/Bi-ASH/output/nc_sims.pdf", width = 10, height = 8)