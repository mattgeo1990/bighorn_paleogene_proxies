library(ggplot2)

file <- paste0(
  "/Users/allen/Documents/GitHub/",
  "bighorn_paleogene_proxies/data/processed/",
  "BHB_multiproxy_final.csv"
)

d <- read.csv(file, na.strings = c("NA", ""))

# Put the two sets of sample means into long format
plot_data <- rbind(
  data.frame(
    dataset = "U-M",
    sample_id = d$MLA_horizon_id,
    T47_C = d$IPLD47_mean_T47_C
  ),
  data.frame(
    dataset = "Caltech",
    sample_id = d$Snell_sample_id,
    T47_C = d$Snell_mean_T47_C
  )
)

plot_data <- subset(plot_data, is.finite(T47_C))
plot_data$dataset <- factor(
  plot_data$dataset,
  levels = c("U-M", "Caltech")
)

# Equal-weight mean of sample means and t-based 95% CI
summarize_group <- function(x) {
  test <- t.test(x, conf.level = 0.95)
  
  data.frame(
    n = length(x),
    mean_T47_C = mean(x),
    sd_T47_C = sd(x),
    se_T47_C = sd(x) / sqrt(length(x)),
    CI95_lower_C = test$conf.int[1],
    CI95_upper_C = test$conf.int[2]
  )
}

summary_table <- do.call(
  rbind,
  lapply(
    split(plot_data$T47_C, plot_data$dataset),
    summarize_group
  )
)

summary_table$dataset <- rownames(summary_table)
rownames(summary_table) <- NULL

summary_table <- summary_table[
  ,
  c(
    "dataset", "n", "mean_T47_C", "sd_T47_C",
    "se_T47_C", "CI95_lower_C", "CI95_upper_C"
  )
]

print(summary_table, digits = 3)

# Raw sample means + group means and 95% CIs
ggplot(plot_data, aes(x = dataset, y = T47_C, color = dataset)) +
  geom_jitter(
    width = 0.10,
    height = 0,
    alpha = 0.65,
    size = 2.5
  ) +
  geom_errorbar(
    data = summary_table,
    aes(
      x = dataset,
      y = mean_T47_C,
      ymin = CI95_lower_C,
      ymax = CI95_upper_C
    ),
    inherit.aes = FALSE,
    width = 0.13,
    linewidth = 0.9,
    color = "black"
  ) +
  geom_point(
    data = summary_table,
    aes(x = dataset, y = mean_T47_C),
    inherit.aes = FALSE,
    shape = 21,
    size = 4,
    stroke = 1,
    fill = "white",
    color = "black"
  ) +
  scale_color_manual(values = c(`U-M` = "#0072B2", Caltech = "#D55E00")) +
  labs(
    x = NULL,
    y = expression(T[47]~degree*C),
    title = expression("Comparison of U-M and Caltech "*T[47]),
    subtitle = "Points are sample means; open circles and bars are group means ± 95% CI"
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none")
