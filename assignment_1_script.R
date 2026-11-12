#### Snipets 
##***************************
## Assignment-1 : Shagufta, Yeshi, Diana, Rachel
##
## 2026-11-02
##
##***************************

## Packages used -------
library(tidyverse)
conflicted::conflicts_prefer(dplyr::filter())
library(viridis)
# + scale_color/fill_viridis_c/d()
theme_set(theme_light())

##Load the data -----

trichinella_spiralis <- read.csv("../data/GCA_001447595.1.csv")
trichinella_zimbabwensis <- read_csv("../data/GCA_001447665.1.csv")

##WRANGLE -----
###make a vector containing strings of non TE names 
ls_non_te = c("simple_repeat", "Simple_repeat",
              "Low_complexity", "Satellite",
              "rRNA", "snRNA", "ARTEFACT",
              "Unknown")

####filter repeat dataset for only TE's 
trichinella_zimbabwensis_te_cleaned <- trichinella_zimbabwensis|> 
  filter(matching_class %notin% ls_non_te)

trichinella_spiralis_te_cleaned <- trichinella_spiralis |> 
  filter(matching_class %notin% ls_non_te)

#####Find TE lengths
trichinella_spiralis_te <- trichinella_spiralis_te_cleaned |> 
  mutate(
    te_length = qry_width,
    te_category = case_when(
      str_detect(matching_class, "^DNA|^RC") ~ "DNA_transposon",
      str_detect(matching_class, "^(LINE|LTR|PLE|SINE)") ~ "Retrotransposon"))

trichinella_zimbabwensis_te <- trichinella_zimbabwensis_te_cleaned |> 
  mutate(
    te_length = qry_width,
    te_category = case_when(
      str_detect(matching_class, "^DNA|^RC") ~ "DNA_transposon",
      str_detect(matching_class, "^(LINE|LTR|PLE|SINE)") ~ "Retrotransposon"))

######Calculate TE category proportions for each species
spiralis_proportions <- trichinella_spiralis_te |>
  count(te_category) |>
  mutate(proportion = n / sum(n),
         species = "T. spiralis")

zimb_proportions <- trichinella_zimbabwensis_te |>
  count(te_category) |>
  mutate(proportion = n / sum(n),
         species = "T. zimbabwensis") 

#######combine both datasets for TE length comparison 
te_fragment_df <- bind_rows(
  trichinella_spiralis_te |>
    mutate(species = "T. spiralis"),
  trichinella_zimbabwensis_te  |> 
    mutate(species = "T. zimbabwensis")) |>
  mutate(
    species = factor(species, levels = c("T. spiralis", "T. zimbabwensis")),
    te_category = factor(te_category,
                         levels = c("DNA_transposon", "Retrotransposon")))

#STATISTICAL ANALYSIS-----
###chi square test for independence of TE category and species 
counts <- matrix(
  c(120, 917,
    207, 1072),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    Species = c("T. spiralis", "T. zimbabwensis"),
    TE_class = c("DNA transposons", "Retrotransposons")))

chisq.test(counts)

####calculate standard-error of proportions
p_spiralis_DNA <- 0.1157184
n_spiralis_DNA <- 120

sqrt(p_spiralis_DNA * (1 - p_spiralis_DNA) / n_spiralis_DNA)

p_spiralis_RT <- 0.8842816
n_spiralis_RT <- 917

sqrt(p_spiralis_RT * (1 - p_spiralis_RT) / n_spiralis_RT)

p_zimb_DNA <- 0.1618452
n_zimb_DNA <- 207

sqrt(p_zimb_DNA * (1 - p_zimb_DNA) / n_zimb_DNA)

p_zimb_RT <- 0.8381548
n_zimb_RT <- 1072

sqrt(p_zimb_RT * (1 - p_zimb_RT) / n_zimb_RT)

#####add st error to proportions dataset
spiralis_proportions_er <- spiralis_proportions |> 
  mutate(st_error = c(0.02920155, 0.0105636))

zimb_proportions_er <- zimb_proportions |> 
  mutate(st_error = c(0.02559923, 0.01124902))

###two way anova for TE length
aov(te_length ~ species * te_category, data = te_fragment_df)


#GRAPH 1 EXPLORATORY ------
###Histogram comparing TE length between both species

figure_1 <- ggplot() +
  geom_freqpoly(
    data = trichinella_zimbabwensis_te,
    aes(x = te_length, colour = "T.zimbabwensis"), 
    alpha = 0.6,
    bins = 200,
    linewidth = 1.0
  ) +
  geom_freqpoly(
    data = trichinella_spiralis_te, 
    aes(x = te_length, colour = "T.spiralis"),
    alpha = 0.6,
    bins = 200,
    linewidth = 1.0 
  ) +
  labs(
    x = "TE Length (bp)",
    y = "TE Count",
    colour = ""
  ) +
  coord_cartesian(xlim= c(0,300)
  ) +
  theme(legend.position = "bottom")

figure_1

#GRAPH 2 FINAL RESULTS----
####barchart for TE class proportion 

figure_2 <- ggplot() +
  geom_col(
    data = spiralis_proportions_er,
    aes(x = te_category, y = proportion, fill = "T.spiralis"),
    position = position_nudge(x = -0.2),
    width = 0.35,
    colour = "black",
    linewidth = 0.5
  ) +
  geom_errorbar(
    data = spiralis_proportions_er,
    aes(x = te_category,
        ymin = proportion - st_error,
        ymax = proportion + st_error),
    position = position_nudge(x = -0.2),
    width = 0.12
  ) +
  geom_col(
    data = zimb_proportions_er,
    aes(x = te_category, y = proportion, fill = "T.zimbabwensis"),
    position = position_nudge(x = 0.2),
    width = 0.35,
    colour = "black",
    linewidth = 0.5
  ) +
  geom_errorbar(
    data = zimb_proportions_er,
    aes(x = te_category,
        ymin = proportion - st_error,
        ymax = proportion + st_error),
    position = position_nudge(x = 0.2),
    width = 0.12
  ) + 
  labs(
    x = "",
    y = "Proportion of Total TE's (%)",
    fill = ""
  ) +
  scale_fill_manual(
    values = c(
      "T.spiralis" = "#E69F00",
      "T.zimbabwensis" = "#1B98E0ff")
  ) +
  scale_y_continuous(labels = scales::percent
  ) +
  theme(legend.position = "bottom")

figure_2

##GRAPH 3 FINAL RESULTS -----
##TE length

figure_3 <- ggplot(te_fragment_df, aes(x = species, y = te_length, fill = species)
) +
  geom_boxplot(outlier.shape = NA, 
               width = 0.7,
               alpha = 0.8, 
               linewidth = 0.4, 
               color = "black"
  ) +
  geom_jitter(width = 0.2,
              size = 0.4, 
              alpha = 0.2
  ) +
  facet_wrap(~ te_category, nrow = 1
  ) +
  scale_y_log10(
    breaks = c(50, 100, 200, 500, 1000, 2000),
  ) + scale_fill_manual(
    values = c("T. spiralis" = "#0072B2",
               "T. zimbabwensis" = "#E69F00"),
    name = NULL
  ) +
  labs(
    x = NULL,
    y = "TE length (bp, log10 scale)"
  ) +
  theme(legend.position = "bottom")

figure_3
