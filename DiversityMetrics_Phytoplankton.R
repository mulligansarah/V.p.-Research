library(tidyverse)
library(tidyr)
library(vegan)
library(dplyr)
library(broom)
library(purrr)
library(ggplot2)
library(ggpubr)
library(viridis)

data <- read.csv('ChemTax_Results_altered.csv')
div <- data |> 
  rowwise() |>
  mutate(richness = specnumber(c_across(`All.Cyano`:`Dinoflagellates`)),
    shan     = diversity(c_across(`All.Cyano`:`Dinoflagellates`)),
    simp     = diversity(c_across(`All.Cyano`:`Dinoflagellates`), "simpson"),
    invsimp  = diversity(c_across(`All.Cyano`:`Dinoflagellates`), "inv")) |>
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

div0phytoplot <- ggplot(subset(plot_data, time == 0), aes(x = metric, y = Mean, fill = Treatment)) +
  geom_col(position = position_dodge(width = 0.7)) +
  geom_text(data = subset(star_positions, time == 0),
            aes(x = metric, y = y, label = sig),
            inherit.aes = FALSE, size = 6) +
  scale_fill_manual(values = c("W" = "darkslateblue", "WV" = "aquamarine3")) +
  theme_classic(base_size = 14)+
  labs(title = "Hour 0", x = 'Metric', y = 'Mean of Metric')+
  theme(plot.title = element_text(face = "italic", hjust = 0.5),
                                  axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
div60phytoplot <- ggplot(subset(plot_data, time == 60), aes(x = metric, y = Mean, fill = Treatment)) +
  geom_col(position = position_dodge(width = 0.7)) +
  geom_text(data = subset(star_positions, time == 60),
            aes(x = metric, y = y, label = sig),
            inherit.aes = FALSE, size = 6) +
  scale_fill_manual(values = c("W" = "darkslateblue", "WV" = "aquamarine3")) +
  theme_classic(base_size = 14)+
  labs(title = "Hour 60", x = 'Metric', y = 'Mean of Metric')+
  theme(plot.title = element_text(face = "italic", hjust = 0.5),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
divphytoplot <- ggarrange(div0phytoplot, div60phytoplot, ncol = 2, common.legend = TRUE, legend = "right")
divphytoplot <- annotate_figure(divphytoplot,
  top = text_grob("Diversity Metrics by Treatment and Time", face = "bold", size = 16, hjust = 0.5)
)

#Other diversity metrics
data_l = data |> 
  pivot_longer(cols = 7:11, 
               names_to = "Taxonomy", 
               values_to = "Chl-Quantity") 
gammaDiv = length(unique(data_l$Taxonomy))
betaDiv = div |> 
  group_by(Treatment, Time) |> 
  summarise(alpha = mean(richness, na.rm = TRUE),
            gamma = gammaDiv,
            beta_a = gamma - alpha)

a <- ggplot(betaDiv, aes(Time, alpha, color = Treatment))+
  geom_point()+
  geom_line()+
  scale_y_continuous(limits = c(0,7))+
  labs(x = 'Hour', y = expression(alpha ~ 'diversity'))+
  theme_classic() +
  scale_color_manual(values = c("W" = "darkslateblue",
                               "WV" = "aquamarine3"))

b <- ggplot(betaDiv, aes(Time, beta_a, color = Treatment))+
  geom_point()+
  geom_line()+
  scale_y_continuous(limits = c(0,7))+
  labs(x = 'month', y = expression(beta ~ 'diversity')) +
  theme_classic() +
  scale_color_manual(values = c("W" = "darkslateblue",
                                "WV" = "aquamarine3"))

alphabetaplot <- ggarrange(a,b,nrow =1, common.legend = TRUE, legend = 'bottom', align = 'h')
alphabetaplot <- annotate_figure(alphabetaplot,
  top = text_grob("Alpha and Beta Diversity of Phytoplankton Over Time",
                  face = "bold",
                  size = 16,
                  hjust = 0.5))

#Dissimilarities
phyto_comm <- data |> 
  group_by(Treatment, Time) |> 
  summarise(across(`All.Cyano`:`Dinoflagellates`, mean)) |>
  unite("SampleID", Treatment, Time, sep = "_") |>  
  column_to_rownames(var = "SampleID")

jacphytodist <- vegdist(phyto_comm, method = "jaccard")
jacphytodist
phyto.nmds.jc <- metaMDS(phyto_comm, distance = "jaccard", k = 2, try = 100)
phyto.nmds.jc
brayphytodist = vegdist(phyto_comm, method = "bray")
brayphytodist
phyto.nmds.bc = metaMDS(phyto_comm, distance = "bray", k = 2, try = 100)

phyto_spp <- data |> 
  group_by(Treatment, Time) |> 
  summarise(across(c(All.Cyano, Dinoflagellates), mean))
nmdsphyto_output <- bind_rows(bc = data.frame(phyto.nmds.bc[["points"]]),jc = data.frame(phyto.nmds.jc[["points"]])) |> 
  mutate(site = rep(unique(paste(phyto_spp$Treatment, phyto_spp$Time, sep = "_")), times = 2),
    Dissimilarity = rep(c("Bray", "Jaccard"), each = length(unique(paste(phyto_spp$Treatment, phyto_spp$Time, sep = "_")))))
nmdsphyto_bc <- subset(nmdsphyto_output, Dissimilarity == "Bray")
nmdsphyto_bc_plot <- ggplot(nmdsphyto_bc, aes(MDS1, MDS2, color = as.factor(site))) +
  geom_point(size = 4) +
  scale_color_viridis_d(option = "turbo") +
  labs(color = "Treatment", title = "Bray-Curtis") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "italic", hjust = 0.5))
nmdsphyto_jc <- subset(nmdsphyto_output, Dissimilarity == "Jaccard")
nmdsphyto_jc_plot <- ggplot(nmdsphyto_jc, aes(MDS1, MDS2, color = as.factor(site))) +
  geom_point(size = 4) +
  scale_color_viridis_d(option = "turbo") +
  labs(color = "Treatment", title = "Jaccard") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "italic", hjust = 0.5))
nmds_phyto <- ggarrange(nmdsphyto_bc_plot, nmdsphyto_jc_plot,ncol = 2,common.legend = TRUE,legend = "bottom")
nmds_phyto <- annotate_figure(nmds_phyto,top = text_grob("NMDS of Phytoplankton Communities", face = "bold", size = 16, hjust = 0.5))

#Partitioning
phyto_spp.2 <- data_l |>
  mutate(PA = if_else(`Chl-Quantity` > 0, 1, 0),
    TreatmentTime = paste(Treatment, Time, sep = "_")) |>
  group_by(TreatmentTime, Taxonomy) |>
  summarise(mean_count = mean(`Chl-Quantity`, na.rm = TRUE),
    mean_PA = mean(PA, na.rm = TRUE),
    .groups = "drop")
phyto_spp.2 <- phyto_spp.2 |>
  mutate(Taxonomy = recode(Taxonomy, "Diatoms...Haptos" = "Diatoms&Haptos"))
Cophyto_dist <- ggplot(phyto_spp.2, aes(x = Taxonomy, y = TreatmentTime, fill = mean_count)) +
  geom_raster() +
  scale_fill_viridis_c(option = 'mako', trans = "log1p") +
  labs(y = 'Treatment & Time', x = 'Taxonomy', fill = 'Mean Chlorophyll Content', title = 'Mean Abundance of Taxonomic Groups with Treatment and Time') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5))