##***************************
## Assignment-1 : Shagufta, Yeshi, Diana, Rachel
## Primary author of this script: Yeshi
## 2026-11-02
##
## Assignment-2 edits by Hasan Jasem (secondary contributor), 2026-11-19:
##   Edit 1. Made the three figures consistent: one shared colour palette for
##           the two species (Fig 2 and Fig 3 previously had the colours
##           swapped), italic species names, cleaner axis/category labels,
##           and simpler figure code that plots combined data instead of
##           adding one layer per species
##   Edit 2. Wrote process_te() to replace the cleaning/classifying code
##           that was repeated for each species
##   Edit 3. Added a data check and summary step (repeats kept after
##           filtering, TE counts and lengths per class) and saved the
##           figures to an output folder
##   Minor.  Defined %notin% in the script so it runs in a fresh R session,
##           printed the ANOVA summary table, consistent read_csv(),
##           removed unused viridis package, shared theme set once at top
##***************************

## Packages used -------
library(tidyverse)
conflicted::conflicts_prefer(dplyr::filter())

## Shared settings -------
# EDIT 1: legend position and italic species names set once here,
# instead of repeating theme(legend.position = "bottom") in every figure
theme_set(theme_light() +
            theme(legend.position = "bottom",
                  legend.text = element_text(face = "italic")))

# EDIT 1: one colour per species, used by every figure (colour-blind friendly)
species_colours <- c("T. spiralis"     = "blue",
                     "T. zimbabwensis" = "orange")

# EDIT 1: readable labels for the TE categories (no underscores in figures)
te_labels <- c(DNA_transposon  = "DNA transposons",
               Retrotransposon = "Retrotransposons")

# MINOR: %notin% is not built into R; defining it here means the script
# runs in a fresh R session, not only when it is already in the saved workspace
`%notin%` <- Negate(`%in%`)

##Load the data -----
# MINOR: both files now use read_csv() so they load the same way
trichinella_spiralis <- read_csv("../data/GCA_001447595.1.csv")
trichinella_zimbabwensis <- read_csv("../data/GCA_001447665.1.csv")

##WRANGLE -----
###make a vector containing strings of non TE names
ls_non_te = c("simple_repeat", "Simple_repeat",
              "Low_complexity", "Satellite",
              "rRNA", "snRNA", "ARTEFACT",
              "Unknown")

# EDIT 2: one function replaces the filter + classify code that was
# written out twice (once per species)
process_te <- function(repeat_df, species_name) {
  repeat_df |>
    filter(matching_class %notin% ls_non_te) |>
    mutate(
      te_length = qry_width,
      te_category = case_when(
        str_detect(matching_class, "^DNA|^RC") ~ "DNA_transposon",
        str_detect(matching_class, "^(LINE|LTR|PLE|SINE)") ~ "Retrotransposon"),
      species = species_name)
}

trichinella_spiralis_te <- process_te(trichinella_spiralis, "T. spiralis")
trichinella_zimbabwensis_te <- process_te(trichinella_zimbabwensis, "T. zimbabwensis")

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
# EDIT 2: species is already added by process_te(), so no mutate() needed here
te_fragment_df <- bind_rows(trichinella_spiralis_te,
                            trichinella_zimbabwensis_te) |>
  mutate(
    species = factor(species, levels = c("T. spiralis", "T. zimbabwensis")),
    te_category = factor(te_category,
                         levels = c("DNA_transposon", "Retrotransposon")))

#DATA CHECK -----
# EDIT 3: check how many repeats were kept after removing non-TEs
filter_check <- tibble(
  species = c("T. spiralis", "T. zimbabwensis"),
  all_repeats = c(nrow(trichinella_spiralis), nrow(trichinella_zimbabwensis)),
  tes_kept = c(nrow(trichinella_spiralis_te), nrow(trichinella_zimbabwensis_te))) |>
  mutate(percent_kept = round(tes_kept / all_repeats * 100, 1))

filter_check

# EDIT 3: summary of TE counts and lengths per species and TE category.
# An NA row here would mean some TEs did not match either category.
te_summary <- te_fragment_df |>
  group_by(species, te_category) |>
  summarise(n_te = n(),
            median_length_bp = median(te_length),
            mean_length_bp = round(mean(te_length), 1),
            .groups = "drop")

te_summary

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

# EDIT 1: combine both species so Figure 2 can be drawn with one layer
te_proportions_er <- bind_rows(spiralis_proportions_er, zimb_proportions_er) |>
  mutate(species = factor(species, levels = c("T. spiralis", "T. zimbabwensis")))

###two way anova for TE length
# MINOR: wrapped in summary() so the ANOVA table (F values, p-values) prints
summary(aov(te_length ~ species * te_category, data = te_fragment_df))


#GRAPH 1 EXPLORATORY ------
###Frequency polygon comparing TE length between both species
# EDIT 1: one layer using the combined data instead of one layer per species;
# colours now match Figures 2 and 3

figure_1 <- ggplot(te_fragment_df, aes(x = te_length, colour = species)) +
  geom_freqpoly(bins = 200, linewidth = 1.0, alpha = 0.8) +
  scale_colour_manual(values = species_colours, name = NULL) +
  coord_cartesian(xlim = c(0, 300)) +
  labs(x = "TE length (bp)", y = "TE count")

figure_1

#GRAPH 2 FINAL RESULTS----
####bar chart for TE class proportion
# EDIT 1: one geom_col + one geom_errorbar with position_dodge replaces the
# four separate layers; same species colours as the other figures

figure_2 <- ggplot(te_proportions_er,
                   aes(x = te_category, y = proportion, fill = species)) +
  geom_col(position = position_dodge(width = 0.8),
           width = 0.7, colour = "black", linewidth = 0.5) +
  geom_errorbar(aes(ymin = proportion - st_error,
                    ymax = proportion + st_error),
                position = position_dodge(width = 0.8),
                width = 0.15) +
  scale_fill_manual(values = species_colours, name = NULL) +
  scale_x_discrete(labels = te_labels) +
  scale_y_continuous(labels = scales::percent,
                     expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Proportion of total TEs")

figure_2

##GRAPH 3 FINAL RESULTS -----
##TE length
# EDIT 1: shared colours (previously opposite to Figure 2), readable facet
# labels, italic species names on the x axis

figure_3 <- ggplot(te_fragment_df |> filter(!is.na(te_category)),
                   aes(x = species, y = te_length, fill = species)) +
  geom_boxplot(outlier.shape = NA,
               width = 0.7,
               alpha = 0.8,
               linewidth = 0.4,
               colour = "black") +
  geom_jitter(width = 0.2, size = 0.4, alpha = 0.2) +
  facet_wrap(~ te_category, nrow = 1, labeller = as_labeller(te_labels)) +
  scale_y_log10(breaks = c(50, 100, 200, 500, 1000, 2000)) +
  scale_fill_manual(values = species_colours, name = NULL) +
  labs(x = NULL, y = "TE length (bp, log10 scale)") +
  theme(axis.text.x = element_text(face = "italic"))

figure_3

#SAVE FIGURES -----
# EDIT 3: save each figure to an output folder so they can be reused
# without re-running the code (creates the folder if it does not exist)
dir.create("../output", showWarnings = FALSE)
ggsave("../output/figure1_te_length_freqpoly.png", plot = figure_1, width = 6, height = 4.5)
ggsave("../output/figure2_te_class_proportion.png", plot = figure_2, width = 6, height = 4.5)
ggsave("../output/figure3_te_length_boxplot.png", plot = figure_3, width = 7, height = 4.5)
