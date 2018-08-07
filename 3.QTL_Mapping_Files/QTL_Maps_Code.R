library(tidyverse)
library(broom)
library(cowplot)
library(grid)
library(gridExtra)


###Compiled the LOD score outputs from IcIM software into one file to create 
###the QTL maps (for plotting figure 2)
##Input the data
data <- read_tsv(file = "All_Results.txt") %>% 
  mutate(Color = Chromosome %% 2)

##This is pulling out the chromosomes where there is at least one significant QTL
##across all the data.
chrsPull <- data %>%
  select(Chromosome,LOD) %>%
  filter(LOD >= 3.0) %>%
  select(Chromosome) %>%
  unique() %>%
  arrange(Chromosome) %>%
  pull(Chromosome)

##Since we have a consensus map, we can pull the markers from all the populations and 
##use that as the order (for plotting). We want the genetic positions to align across
##the populations as we plot
unique_order <- data %>%
  filter(Chromosome %in% chrsPull) %>%
  select(Chromosome,Position) %>%
  unique() %>%
  arrange(Chromosome,Position) %>%
  mutate(Order = row_number())
  

##This pulls the order of the markers into the data frame so we can plot
data_formatted <- data %>%
  filter(Chromosome %in% chrsPull) %>%
  left_join(unique_order,by = c("Chromosome","Position"))

##This calculates the chromosome midpoints and max. It is for plotting purposes, so we have
##the chromsome name in the middle of the chr and a break there
midpoints <- data_formatted %>%
  select(Chromosome,Order,LOD) %>%
  group_by(Chromosome) %>%
  summarise(Max = max(Order),
            Min = min(Order),
            Max_LOD = max(LOD)) %>%
  #filter(Max_LOD >= 3.0) %>%
  mutate(Break = round(((Max-Min)/2 + Min),0)) %>%
  ungroup()


##This is going to plot each LOD map (4) to a data.frame. We are going to add stuff specifically
##to each map            
maps <- data_formatted %>%
  group_by(Population,TraitName) %>%
  do(
    plots = ggplot(data = .) +
      geom_hline(yintercept = 3, color = "red", size = .3) +
      geom_line(aes(x = Order, y = LOD, group = Chromosome)) +
      theme_bw() + 
      scale_y_continuous(name = " \n",
                         #limits = c(0,max(.$LOD) + 1),
                         limits = c(0,8),
                         breaks = seq(from = 0, to = 7, by = 1))  +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    )

##This is going to add the chromsome breaks and chromsome lines
for (i in 1:dim(maps)[1]) {
  holder <- midpoints
  maps$plots[i][[1]] <- maps$plots[i][[1]] +
    scale_x_continuous(breaks = holder$Break,
                       labels = holder$Chromosome,
                       expand = c(.005,.005),
                       limits = c(0,max(unique_order$Order)),
                       name = element_blank()) +
    geom_vline(data = holder,
               aes(xintercept = Max),color = "lightGray", 
               size = .5, linetype = "dotted")
}

##This is going to actually plot the maps and add some labels
popLabel_Size <- 7
locDiff <- .04
allPlots <- plot_grid(maps$plots[1][[1]],
          maps$plots[2][[1]],
          maps$plots[3][[1]],
          maps$plots[4][[1]],
          ncol = 1) + 
  # draw_plot_label(label = "Selected Chromosomes",
  #                                    x = .55,y = .005, 
  #                                    size = 12,
  #                                    hjust = 0.5,vjust = 0) +
  draw_plot_label(c("A","B","C","D"),
                                  c(0,0,0,0),
                                  c(1,.75,.5,.25),
                  size = 12) +
  draw_label(label = "LOD",
             x = .05, y = .5,
             size = 13,
             hjust = 0,vjust = 0,
             angle = 90,
             fontface = "bold") +
  draw_label(label = expression(B73~Population:~F[2]),
             size = popLabel_Size,
             x = .17,
             y = 1 - 0 - locDiff,
             fontface = "bold") +
  draw_label(label = expression(B73~Population:~F[2:3]),
             size = popLabel_Size,
             x = .17,
             y = 1 - .25 - locDiff,
             fontface = "bold") +
  draw_label(label = expression(Mo17~Population:~F[2]),
             size = popLabel_Size,
             x = .17,
             y = 1 - .50 - locDiff,
             fontface = "bold") +
  draw_label(label = expression(Mo17~Population:~F[2:3]),
             size = popLabel_Size,
             x = .17,
             y = 1 - .75 - locDiff,
             fontface = "bold")

allPlots <- add_sub(allPlots,"Selected Chromosomes", 
                    vpadding = grid::unit(0, "lines"),
                    y = 1.3, x = .4, hjust = 0)
##This plots and saves figure 2 from paper
plot_grid(allPlots)


# ggsave(filename = "QTL_Maps_consensus_v2.pdf",
#        plot = plot_grid(allPlots),
#        dpi = 600,
#        width = 7,
#        height = 5)


##############################Supplemental Figure #######################
##This is going to create supplemental figure 2 which compares the single marker scan
##with the LOD mapping

##Created LOD map with all chromosomes
##Figures out the midpoints now for the LOD data with all the chromosomes
midpoints_all <- data %>%
  select(Population,Chromosome,Position) %>%
  unique() %>%
  arrange(Population,Chromosome,Position) %>%
  group_by(Population) %>%
  mutate(Order = row_number()) %>%
  ungroup() %>%
  group_by(Population, Chromosome) %>%
  summarise(Max = max(Order),
            Min = min(Order),
            Break = (Max + Min) / 2) %>% 
  ungroup()

max_all <- midpoints_all %>% 
  select(Population,Chromosome,Max,Min)
##This creates a plotting function to plot a genetic map given a population name and traitID
LOD_plotter <- function(data,Pop,Trait, midpoints) {
  data <- data
  Pop <- Pop
  Trait <- Trait
  midpoints <- midpoints %>%
    filter(Population == Pop)
  data %>%
    filter(Population == Pop,
           TraitID == Trait) %>%
    arrange(Chromosome,Position) %>%
    mutate(Order = row_number()) %>%
    ggplot() +
    geom_vline(data = midpoints %>%
                 filter(Chromosome != 10),
               aes(xintercept = Max),
               linetype = "dotted",
               size = 1,
               color = "light gray") +
    geom_line(aes(x = Order, 
                  y = LOD, 
                  group = Chromosome, 
                  color = factor(Color)),
              size = .75,
              show.legend = FALSE) +
    xlab("Chromosome") +
    scale_x_continuous(labels = c(1:10),
                       breaks = midpoints %>%
                         pull(Break),
                       expand = c(0.025,0)) +
    scale_color_manual(values = c("black","blue")) +
    theme(text = element_text(size = 8),
          axis.text = element_text(size = 8))
}

##Creating the supplementary comparison file of the single marker scan
##with the lnikage mapping results
##This imports the single marker scan data
b73_single <- read_tsv("Single_Marker_Scan/pvalues_B73.txt")
mo17_single <- read_tsv("Single_Marker_Scan/pvalues_Mo17.txt")

##Combines the data into a single data frame
man_data <- b73_single %>%
  mutate(Population = "B73") %>%
  bind_rows(mo17_single %>%
              mutate(Population = "Mo17")) %>%
  gather(Trait,P,-Population,-SNP,-CHR,-BP) %>%
  rename(TraitID = Trait) %>%
  mutate(TraitID = ifelse(TraitID == "F2",1,2))

##Determines the midpoints for that data again for plotting purposes
##It also pulls in the max position for the genetic maps to align the maps so the
##chromosomes overlap

##This is the data used to recalculate the position along the x-axis for the single 
##marker data
renumber_man <- man_data %>%
  select(Population,CHR,BP) %>%
  unique() %>%
  group_by(Population) %>%
  mutate(Order = row_number()) %>% 
  group_by(Population, CHR) %>%
  summarise(Max2 = max(Order),
            Min2 = min(Order)) %>% 
  ungroup() %>% 
  mutate(Chromosome = CHR) %>% 
  left_join(max_all,by = c("Population","Chromosome")) %>%
  ungroup() %>%
  mutate(Ratio = (Max - Min)/(Max2 - Min2)) %>% 
  select(Population,CHR,Min,Ratio)

##This creates the mid points for this new ordered data
midpoints_man <- renumber_man %>% 
  left_join(man_data %>% 
              select(Population,CHR,BP) %>% 
              unique(),
            by = c("Population","CHR")) %>% 
  arrange(Population,CHR,BP) %>% 
  group_by(Population,CHR) %>%
  mutate(Order = (((row_number() - 1) * Ratio) + Min)) %>% 
  summarise(Max = max(Order),
            Min = min(Order),
            Break = (Max + Min)/2) %>% 
  ungroup()



##Creates a plotting function for the single marker scan data that
##uses pop and tratID information
manhattan_plotter <- function(data,Pop,Trait, midpoints) {
  data <- data
  Pop <- Pop
  Trait <- Trait
  midpoints <- midpoints %>%
    filter(Population == Pop)
  data %>%
    filter(Population == Pop,
           TraitID == Trait) %>%
    arrange(CHR,BP) %>%
    left_join(renumber_man,by = c("Population","CHR")) %>%
    arrange(Population,CHR,BP) %>% 
    group_by(Population,CHR) %>%
    mutate(Order = (((row_number() - 1) * Ratio) + Min),
           Color = CHR %% 2) %>% 
    ungroup() %>% 
    ggplot() +
    geom_vline(data = midpoints %>%
                 filter(CHR != 10),
               aes(xintercept = Max),
               linetype = "dotted",
               size = .75,
               color = "light gray") +
    geom_point(aes(x = Order, y = -log10(P),
                   color = factor(Color)),
               show.legend = FALSE,
               size = .5) +
    xlab("Chromosome") +
    scale_x_continuous(labels = c(1:10),
                       breaks = midpoints %>%
                         pull(Break),
                       expand = c(0.025,0)) +
    scale_color_manual(values = c("black","blue")) +
    theme(text = element_text(size = 8),
          axis.text = element_text(size = 8))
}

##This actually plots each one of the maps and combines with with the single marker scan map
plot1 <- 
  plot_grid(LOD_plotter(data = data,
                        Pop = "B73",
                        Trait = 1,
                        midpoints = midpoints_all),
            manhattan_plotter(data = man_data,
                              Pop = "B73",
                              Trait = 1,
                              midpoints = midpoints_man),
            ncol = 1)
plot2 <- 
  plot_grid(LOD_plotter(data = data,
                        Pop = "B73",
                        Trait = 2,
                        midpoints = midpoints_all),
            manhattan_plotter(data = man_data,
                              Pop = "B73",
                              Trait = 2,
                              midpoints = midpoints_man),
            ncol = 1)
plot3 <- 
  plot_grid(LOD_plotter(data = data,
                        Pop = "Mo17",
                        Trait = 1,
                        midpoints = midpoints_all),
            manhattan_plotter(data = man_data,
                              Pop = "Mo17",
                              Trait = 1,
                              midpoints = midpoints_man),
            ncol = 1)
plot4 <- 
  plot_grid(LOD_plotter(data = data,
                        Pop = "Mo17",
                        Trait = 2,
                        midpoints = midpoints_all),
            manhattan_plotter(data = man_data,
                              Pop = "Mo17",
                              Trait = 2,
                              midpoints = midpoints_man),
            ncol = 1)



##This combines the individual plot and makes one big plot seen in supplemental figure 2
popLabel_Size <- 5
locDiff <- .01
allPlots_man_B73 <- plot_grid(plot1,
                      plot2,
                      ncol = 1,
                      labels = c("A","B"),
                      label_size = 10) +
  draw_label(label = expression(B73~Population:~F[2]),
             size = popLabel_Size,
             x = .14,
             y = 1 - locDiff) +
  draw_label(label = expression(B73~Population:~F[2:3]),
             size = popLabel_Size,
             x = .14,
             y = .5 - locDiff)
  
allPlots_man_Mo17 <- plot_grid(plot3,
                               plot4,
                               ncol = 1,
                               labels = c("C","D"),
                               label_size = 10) +
  draw_label(label = expression(Mo17~Population:~F[2]),
             size = popLabel_Size,
             x = .14,
             y = 1 - locDiff) +
  draw_label(label = expression(Mo17~Population:~F[2:3]),
             size = popLabel_Size,
             x = .14,
             y = .5 - locDiff)

plot_grid(allPlots_man_B73,
          allPlots_man_Mo17,
          ncol = 1)
ggsave(filename = "Combined_ComparisonFigure_Supplemental_color.pdf",
       plot = plot_grid(allPlots_man_B73,
                        allPlots_man_Mo17,
                        ncol = 1),
       dpi = 600,
       width = 7,
       height = 9)
