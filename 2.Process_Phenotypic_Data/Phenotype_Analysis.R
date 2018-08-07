library(tidyverse)
library(broom)
library(grid)
library(gridExtra)
library(cowplot)

##Read in the F2 data to summarize the populations
data <- read_delim(file = "All_Phenotypes_F2.txt",
                   delim = "\t")


##Split the reciprocal populations into B73 and Mo17
data_formatted <- data %>%
  mutate(Type = ifelse(Population %in% c("B73/PHW30","PHW30/B73"),"B73",
                       ifelse(Population %in% c("Mo17/PHW30","PHW30/Mo17"),"Mo17",Population)))

##Summarize the reciprocal populations
data_summary_pops <- data_formatted %>%
  filter(!Population %in% c("F1","Parent")) %>%
  group_by(Population) %>%
  summarise(Avg = mean(LA, na.rm = TRUE),
            SD = sqrt(var(LA,na.rm = TRUE)),
            n = n(),
            Min = min(LA,na.rm = TRUE),
            Max = max(LA, na.rm = TRUE))

##Summarize the combined populations
data_summary_combined_pops <- data_formatted %>%
  filter(!Population %in% c("F1","Parent")) %>%
  group_by(Type) %>%
  summarise(Avg = mean(LA, na.rm = TRUE),
            SD = sqrt(var(LA,na.rm = TRUE)),
            n = n(),
            Min = min(LA,na.rm = TRUE),
            Max = max(LA, na.rm = TRUE))

##Summarize the F1 lines
data_summary_F1_Parents <- data_formatted %>%
  filter(Population %in% c("F1","Parent")) %>%
  group_by(Genotype) %>%
  summarise(Avg = mean(LA, na.rm = TRUE),
            SD = sqrt(var(LA,na.rm = TRUE)),
            n = n(),
            Min = min(LA,na.rm = TRUE),
            Max = max(LA, na.rm = TRUE))

# #Write these summaries to text files
# write_delim(data_summary_pops, "F2_Population_Summary.txt",
#             delim = "\t")
# write_delim(data_summary_F1_Parents, "Parents_F1_Summary.txt",
#             delim = "\t")
# write_delim(data_summary_combined_pops, "Combined_Groups_Summary.txt",
#             delim = "\t")

f2_hist <- data_formatted %>%
    filter(Type %in% c("B73","Mo17")) %>%
    mutate(Type = ifelse(Type == "B73","B73~Population","Mo17~Population")) %>%
    select(-Population) %>%
    rename(Population = Type) %>%
    mutate(Generation = 'Generation:F[2]') %>%
    ggplot() +
    geom_histogram(aes(LA, ..density..), binwidth = 5) +
    facet_wrap(~ Population + Generation,
               nrow = 1,
               labeller = label_parsed) +
    theme_bw() +
    xlab("Leaf Angle (Degrees)") +
    scale_x_continuous(limits = c(40,80), breaks = seq(40,80,5))

####Test the hypothesis that the reciprocal crosses are equal
b73_test <- t.test(x = data_formatted %>% filter(Population == "B73/PHW30") %>% pull(LA) %>% na.omit(),
       y = data_formatted %>% filter(Population == "PHW30/B73") %>% pull(LA) %>% na.omit(),
       var.equal = FALSE)
b73_test$estimate[1] - b73_test$estimate[2]
b73_test$p.value

mo17_test <- t.test(x = data_formatted %>% filter(Population == "Mo17/PHW30") %>% pull(LA) %>% na.omit(),
       y = data_formatted %>% filter(Population == "PHW30/Mo17") %>% pull(LA) %>% na.omit(),
       var.equal = FALSE)

mo17_test$estimate[1] - mo17_test$estimate[2]
mo17_test$p.value

##########F3 Summary###########
##Read in the F2 and F3 data file, this file has all the F3 phenotypes and also the F2 in one
##file
f3_data <- read_delim("All_F2_F3_LA_Data.txt",
                       delim = "\t") %>%
  mutate(Type = ifelse(Population %in% c("B73/PHW30","PHW30/B73"),"B73",
                       ifelse(Population %in% c("Mo17/PHW30","PHW30/Mo17"),"Mo17",Population)))

##Summarizes the F3 combined population
f3_pop_summary <- f3_data %>%
  select(Genotype,Type,Population,F3) %>%
  na.omit() %>%
  group_by(Type) %>%
  summarise(Avg = mean(F3,na.rm = TRUE),
            SD = sqrt(var(F3,na.rm = TRUE)),
            n = n(),
            Min = min(F3,na.rm = TRUE),
            Max = max(F3, na.rm = TRUE))

##Summarizes the F2 groups, that made it to next generation combined population
f3_pop_group <- f3_data %>%
  select(Genotype,Type,Population,Group,F2) %>%
  na.omit() %>%
  group_by(Type,Group) %>%
  summarise(Avg = mean(F2,na.rm = TRUE),
            SD = sqrt(var(F2,na.rm = TRUE)),
            n = n(),
            Min = min(F2,na.rm = TRUE),
            Max = max(F2, na.rm = TRUE))

##Test if there are any significant difference between hybrids
data_f1_test <- data_formatted %>%
  filter(Population == "F1") %>%
  select(Genotype, LA) %>%
  do(tidy(aov(LA ~ Genotype, .)))

##Makes sure that the F2 selected groups are significantly different
data_group_test <- data_formatted %>%
  filter(!Population %in% c("F1","Parent")) %>%
  select(Population, LA) %>%
  do(tidy(aov(LA ~ Population, .)))


##Writes the F3 pop summary to a text file
# write_delim(f3_pop_summary, "Phenotype_Analysis/Summary/F3_Pops_Summary.txt",
#             delim = "\t")





######The following are to plot the supplemental figure 1#########

##Distribution of phenotypic data
##
parental_data_B73 <- data_summary_F1_Parents %>%
  select(Genotype,Avg) %>%
  filter(Genotype %in% c("B73","PHW30")) %>%
  mutate(Population = "B73~Population")
parental_data_Mo17 <- data_summary_F1_Parents %>%
  select(Genotype,Avg) %>%
  filter(Genotype %in% c("Mo17","PHW30")) %>%
  mutate(Population = "Mo17~Population")


parental_data_formatted_f2 <- parental_data_B73 %>%
  mutate(Generation = "Generation:F[2]") %>%
  bind_rows(parental_data_Mo17 %>%
              mutate(Generation = "Generation:F[2]"))
            
parental_data_formatted_f3 <- parental_data_B73 %>%
  mutate(Generation = "Generation:F[2:3]") %>%
  bind_rows(parental_data_Mo17 %>%
              mutate(Generation = "Generation:F[2:3]"))           
  



f3_hist <- f3_data %>%
  select(-Group,-F2) %>%
  rename(LA = F3) %>%
  mutate(Type = ifelse(Population %in% c("B73/PHW30","PHW30/B73"),"B73~Population",
                       ifelse(Population %in% c("Mo17/PHW30","PHW30/Mo17"),"Mo17~Population",
                              Population)),
         Generation = 'Generation:F[2:3]') %>%
  select(-Population) %>%
  rename(Population = Type) %>%
  ggplot() +
  geom_histogram(aes(LA, ..density..), binwidth = 5) +
  facet_wrap(~ Population + Generation,
             nrow = 1,
             labeller = label_parsed) +
  theme_bw() +
  xlab("Leaf Angle (Degrees)") +
  scale_x_continuous(limits = c(40,80), breaks = seq(40,80,5)) +
  scale_y_continuous(limits = c(0,0.09)) +
  theme(strip.text.x = element_text(size = 6,
                                    margin = margin(.15,0,.1,0,"cm")),
        strip.background = element_rect(colour = "black", fill = NA))



f2_hist_parents <- f2_hist +
  scale_y_continuous(limits = c(0,0.09)) +
  geom_point(data = parental_data_formatted_f2,
             aes(x = Avg),
             y = 0.075,
             size = 2) +
  geom_text(data = parental_data_formatted_f2,
            aes(x = Avg, y = 0.0875, label = Genotype),
            size = 2,
            fontface = "bold") +
    theme(strip.text.x = element_text(size = 6,
                                      margin = margin(.15,0,.1,0,"cm")),
          strip.background = element_rect(colour = "black", fill = NA))
plot_grid(f2_hist_parents, f3_hist, nrow = 2)

# save_plot("phenotypic_histogram3.png",
#           base_height = 4, base_width = 7,
#           plot_grid(f2_hist_parents, f3_hist, nrow = 2),
#           base_aspect_ratio = .6
#           )

f3_hist <- f3_data %>%
  select(-F3) %>%
  rename(LA = F2) %>%
  mutate(Type = ifelse(Population %in% c("B73/PHW30","PHW30/B73"),"B73~Population",
                       ifelse(Population %in% c("Mo17/PHW30","PHW30/Mo17"),"Mo17~Population",
                              Population)),
         Generation = 'Generation:F[2]') %>%
  select(-Population) %>%
  rename(Population = Type) %>%
  ggplot() +
  geom_histogram(aes(LA, ..density..), binwidth = 5) +
  facet_wrap(~ Population + Generation,
             nrow = 1,
             labeller = label_parsed)
