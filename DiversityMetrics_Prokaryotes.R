library(tidyverse)
library(tidyr)
library(vegan)
library(dplyr)
library(broom)
library(purrr)
library(ggplot2)
library(ggpubr)
library(viridis)

data <- read.csv('level-4_altered.csv')
div <- data |> 
  rowwise() |>
  mutate(richness = specnumber(c_across(`d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Vibrionales`:`d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Gammaproteobacteria`)),
         shan     = diversity(c_across(`d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Vibrionales`:`d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Gammaproteobacteria`)),
         simp     = diversity(c_across(`d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Vibrionales`:`d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Gammaproteobacteria`), "simpson"),
         invsimp  = diversity(c_across(`d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Vibrionales`:`d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Gammaproteobacteria`), "inv")) |>
  ungroup()
average_div <- div |> 
  group_by(Treatment, Time) |> 
  summarize(across(richness:invsimp, list(mean = mean, sd = sd)))
div0 <- filter(div, Time == 0)
div60 <- filter(div, Time == 60)
metrics <- c("richness", "shan", "simp", "invsimp")
run_tests <- function(df) {
  lapply(metrics, function(m)
    t.test(df[[m]] ~ df$Treatment)
  ) |> setNames(metrics)
}
tests_time0  <- run_tests(div0)
tests_time60 <- run_tests(div60)

results_time0 <- map_dfr(metrics,
                         ~ tidy(tests_time0[[.x]]),
                         .id = "metric") |>
  mutate(time = 0)

results_time60 <- map_dfr(metrics,
                          ~ tidy(tests_time60[[.x]]),
                          .id = "metric") |>
  mutate(time = 60)
results_table <- bind_rows(results_time0, results_time60)
results_table <- results_table |> select(time, everything())
results_table$metric <- metrics[as.integer(results_table$metric)]
results_table <- results_table |>
  mutate(sig = case_when(p.value < 0.001 ~ "***",
                         p.value < 0.01  ~ "**",
                         p.value < 0.05  ~ "*"))
plot_data <- results_table |>
  select(time, metric, estimate1, estimate2, sig) |>
  pivot_longer(cols = c(estimate1, estimate2),
               names_to = "Treatment",
               values_to = "Mean") |>
  mutate(Treatment = recode(Treatment,
                            estimate1 = "W",
                            estimate2 = "WV"))
star_positions <- plot_data |>
  group_by(time, metric) |>
  summarise(y = max(Mean) * 1.05,   # slightly above tallest bar
            sig = unique(sig),
            .groups = "drop")

div0proplot <- ggplot(subset(plot_data, time == 0), aes(x = metric, y = Mean, fill = Treatment)) +
  geom_col(position = position_dodge(width = 0.7)) +
  geom_text(data = subset(star_positions, time == 0),
            aes(x = metric, y = y, label = sig),
            inherit.aes = FALSE, size = 6) +
  scale_fill_manual(values = c("W" = "darkslateblue", "WV" = "aquamarine3")) +
  theme_classic(base_size = 14)+
  labs(title = "Hour 0", x = 'Metric', y = 'Mean of Metric')+
  theme(plot.title = element_text(face = "italic", hjust = 0.5),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
div60proplot <- ggplot(subset(plot_data, time == 60), aes(x = metric, y = Mean, fill = Treatment)) +
  geom_col(position = position_dodge(width = 0.7)) +
  geom_text(data = subset(star_positions, time == 60),
            aes(x = metric, y = y, label = sig),
            inherit.aes = FALSE, size = 6) +
  scale_fill_manual(values = c("W" = "darkslateblue", "WV" = "aquamarine3")) +
  theme_classic(base_size = 14)+
  labs(title = "Hour 60", x = 'Metric', y = 'Mean of Metric')+
  theme(plot.title = element_text(face = "italic", hjust = 0.5),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
divproplot <- ggarrange(div0proplot, div60proplot, ncol = 2, common.legend = TRUE, legend = "right")
divproplot <- annotate_figure(divproplot,
                                top = text_grob("Diversity Metrics of Prokaryotes by Treatment and Time", face = "bold", size = 16, hjust = 0.5))


#Other diversity metrics
data_l = data |> 
  pivot_longer(cols = 4:409, 
               names_to = "Order", 
               values_to = "NumberofASVs")
gammaDiv = length(unique(data_l$Order))
betaDiv = div |> 
  group_by(Treatment, Time) |> 
  summarise(alpha = mean(richness, na.rm = TRUE),
            gamma = gammaDiv,
            beta_a = gamma - alpha)
a <- ggplot(betaDiv, aes(Time, alpha, color = Treatment))+
  geom_point()+
  geom_line()+
  scale_y_continuous(limits = c(0,400))+
  labs(x = 'Hour', y = expression(alpha ~ 'diversity'))+
  theme_classic() +
  scale_color_manual(values = c("W" = "darkslateblue",
                                "WV" = "aquamarine3"))

b <- ggplot(betaDiv, aes(Time, beta_a, color = Treatment))+
  geom_point()+
  geom_line()+
  scale_y_continuous(limits = c(0,400))+
  labs(x = 'Hour', y = expression(beta ~ 'diversity')) +
  theme_classic() +
  scale_color_manual(values = c("W" = "darkslateblue",
                                "WV" = "aquamarine3"))

alphabetaproplot <- ggarrange(a,b,nrow =1, common.legend = TRUE, legend = 'bottom', align = 'h')
alphabetaproplot <- annotate_figure(alphabetaproplot,
                                 top = text_grob("Alpha and Beta Diversity of Prokaryotes Over Time",
                                                 face = "bold",
                                                 size = 16,
                                                 hjust = 0.5))

#Dissimilarities
pro_comm <- data |> 
  group_by(Treatment, Time) |> 
  summarise(across(`d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Vibrionales`:`d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Gammaproteobacteria`, mean)) |>
  unite("SampleID", Treatment, Time, sep = "_") |>  
  column_to_rownames(var = "SampleID")

jacprodist <- vegdist(pro_comm, method = "jaccard")
jacprodist
pro.nmds.jc <- metaMDS(pro_comm, distance = "jaccard", k = 2, try = 100)
pro.nmds.jc
brayprodist = vegdist(pro_comm, method = "bray")
brayprodist
pro.nmds.bc = metaMDS(pro_comm, distance = "bray", k = 2, try = 100)

pro_spp <- data |> 
  group_by(Treatment, Time) |> 
  summarise(across(c(d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Vibrionales:d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Gammaproteobacteria), mean))
nmdspro_output <- bind_rows(bc = data.frame(pro.nmds.bc[["points"]]),jc = data.frame(pro.nmds.jc[["points"]])) |> 
  mutate(site = rep(unique(paste(pro_spp$Treatment, pro_spp$Time, sep = "_")), times = 2),
         Dissimilarity = rep(c("Bray", "Jaccard"), each = length(unique(paste(pro_spp$Treatment, pro_spp$Time, sep = "_")))))
nmdspro_bc <- subset(nmdspro_output, Dissimilarity == "Bray")
nmdspro_bc_plot <- ggplot(nmdspro_bc, aes(MDS1, MDS2, color = as.factor(site))) +
  geom_point(size = 4) +
  scale_color_viridis_d(option = "turbo") +
  labs(color = "Treatment", title = "Bray-Curtis") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "italic", hjust = 0.5))
nmdspro_jc <- subset(nmdspro_output, Dissimilarity == "Jaccard")
nmdspro_jc_plot <- ggplot(nmdspro_jc, aes(MDS1, MDS2, color = as.factor(site))) +
  geom_point(size = 4) +
  scale_color_viridis_d(option = "turbo") +
  labs(color = "Treatment", title = "Jaccard") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "italic", hjust = 0.5))
nmds_pro <- ggarrange(nmdspro_bc_plot, nmdspro_jc_plot,ncol = 2,common.legend = TRUE,legend = "bottom")
nmds_pro <- annotate_figure(nmds_pro,top = text_grob("NMDS of Prokaryote Communities", face = "bold", size = 16, hjust = 0.5))
