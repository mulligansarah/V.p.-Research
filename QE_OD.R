library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(purrr)
library(broom)
library(emmeans)

data <- read.csv('QuantumEfficiency_OD.csv')

data <- data |> 
  mutate(Hour = as.factor(Hour)) |> 
  ungroup()
QE <- ggplot(data, aes(Hour, QuantumEfficiency, fill = Treatment)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  theme_classic() +
  scale_fill_manual(values = c("W" = "darkslateblue", "WV" = "aquamarine3")) +
  labs(title = "Quantum Efficiency for each Treatment over 60 Hours", x = 'Hour', y = 'Quantum Efficiency (Fv/Fm)')+
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

anova_model <- aov(QuantumEfficiency ~ Treatment * Hour, data = data)
summary(anova_model)
posthoc <- emmeans(anova_model, pairwise ~ Treatment | Hour)
posthoc$contrasts

sig_table <- as.data.frame(posthoc$contrasts) |> 
  mutate(sig = case_when(p.value < 0.001 ~ "***", p.value < 0.01  ~ "**", p.value < 0.05  ~ "*")) |>
  select(Hour, contrast, estimate, p.value, sig)

star_positions <- data |>
  group_by(Hour) |>
  summarize(y = max(QuantumEfficiency, na.rm = TRUE) * 1.05)
sig_plot_data <- left_join(star_positions, sig_table, by = "Hour")
QE +
  geom_text(data = sig_plot_data,aes(x = Hour, y = y, label = sig),inherit.aes = FALSE,size = 6)