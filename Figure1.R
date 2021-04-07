## FIGURE 1 SCRIPT
# Emily Bean
# 1/2021


## multi-panel with plots of quantity, 260/280, 260/230

## ---- get data ----

sessionInfo()
source("ColorScript.R")
require(tidyverse)
require(ggpubr)

# read in data and clean
dat <- read.table("DNA_Quantities.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>% 
  # replace "too low" with "0"
  mutate(QuBit = as.numeric(str_replace_all(QuBit, "too low", "0")),
         # get Kit information from the SampleID
         Kit = as.factor(sapply(str_split(SampleID, "_"), '[', 1)),
         # get Chicken information from the SampleID
         Replicate = sapply(str_split(SampleID, "_"), '[', 2),
         # do some massaging to make all controls uniformly named
         Replicate = case_when(
           Replicate %in% "PC" ~ "PosControl",
           Replicate %in% "NC" ~ "NegControl",
           Replicate != "PC" && Replicate != "NC" ~ Replicate
         ),
         # replace negative numbers with 0
         Nanodrop_cleanedPCR = replace(Nanodrop_cleanedPCR, which(Nanodrop_cleanedPCR < 0), 0),
         Nanodrop_260280_cleanedPCR = replace(Nanodrop_260280_cleanedPCR, which(Nanodrop_260280_cleanedPCR < 0), 0),
         Nanodrop_260230_cleanedPCR = replace(Nanodrop_260230_cleanedPCR, which(Nanodrop_260230_cleanedPCR < 0), 0),
         Nanodrop_extrcDNA = replace(Nanodrop_extrcDNA, which(Nanodrop_extrcDNA < 0), 0),
         Nanodrop_260280_extrcDNA = replace(Nanodrop_260280_extrcDNA, which(Nanodrop_260280_extrcDNA < 0), 0),
         Nanodrop_260230_extrcDNA = replace(Nanodrop_260230_extrcDNA, which(Nanodrop_260230_extrcDNA < 0), 0)) %>% 
  # rename kits to add a space between
  ### APRIL 2021 UPDATE: RENAME "KIT" TO "METHOD"
  mutate(Method = case_when(
    Kit %in% "Kit1" ~ "Method 1",
    Kit %in% "Kit2" ~ "Method 2",
    Kit %in% "Kit3" ~ "Method 3",
    Kit %in% "Kit4" ~ "Method 4"
  ))

# create dataframe without controls
datnocon <- dat %>% filter(!Replicate == "PosControl" & !Replicate == "NegControl") 
  

### ---- panel a; quantity ----

# dot plot
a <- ggdotplot(data = datnocon, x = "Method", y = "QuBit",
          fill = "Method", binwidth = 1, size = 1.5,
          add = c("median_iqr"), error.plot = "pointrange", add.params = list(size = 1),
          xlab = FALSE, ylab = "DNA quantity (ng/\u03BCL)", legend = "none") +
  # add manual colors
  scale_fill_manual(values = methodcolswithspaces) +
  # add extra white space for asterisks
  scale_y_continuous(expand = expansion(mult = 0, add = c(0.5, 3.5))) +
  #stat_pvalue_manual(statdf, label = "p.adj.signif", label.size = 7)
  geom_text(label = "a", x = 1, y = 15, size = 6) +
  geom_text(label = "a", x = 2, y = 15, size = 6) +
  geom_text(label = "b", x = 3, y = 25, size = 6) +
  geom_text(label = "b", x = 4, y = 25, size = 6)


## ---- panel b; 260/280 ----

# create scale column for difference from ideal
datnocon <- datnocon %>% 
  mutate(scale260280 = Nanodrop_260280_extrcDNA - 1.80)

# dot plot
b <- ggdotplot(data = datnocon, x = "Method", y = "scale260280",
          fill = "Method", binwidth = 0.1, size = 1.4,
          add = c("median_iqr"), error.plot = "pointrange", #add.params = list(size = 1.5),
          xlab = FALSE, ylab = "Difference from 1.80 ideal", legend = "none") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", size = 1) +
  # add manual colors
  scale_fill_manual(values = methodcolswithspaces)+
  # add significance labels
  geom_text(label = "*", x = 2, y = 0.5, size = 8) +
  geom_text(label = "*", x = 4, y = 0.5, size = 8)


## ---- panel c; 260/230 ----

# create scale column for difference from ideal
datnocon <- datnocon %>% 
  mutate(scale260230 = Nanodrop_260230_extrcDNA - 2.00)

# dot plot
c <- ggdotplot(data = datnocon, x = "Method", y = "scale260230",
               fill = "Method", binwidth = 0.8, size = 0.13, 
               add = c("median_iqr"), error.plot = "pointrange", #add.params = list(size = 1.5),
               xlab = FALSE, ylab = "Difference from 2.00 ideal", legend = "none") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", size = 1) +
  # add manual colors
  scale_fill_manual(values = methodcolswithspaces) +
  scale_x_discrete(expand = expansion(mult = -0.8, add = c(0, 3)) ) +
  # add significance labels
  geom_text(label = "*", x = 1, y = -01.5, size = 8) +
  geom_text(label = "*", x = 2, y = -0.7, size = 8) +
  geom_text(label = "*", x = 3, y = -1.5, size = 8) +
  geom_text(label = "*", x = 4, y = -1.5, size = 8)


## ---- add together ----

g1 <- ggarrange(a, b, c, ncol = 3, nrow = 1,
                widths = c(1, 1.3, 1.5),
          labels = c("a", "b", "c"))

ggsave(filename = "./Manuscript/Figure1_Renamed.jpeg", plot = g1,  device = "jpeg", dpi = 600, height = 4, width = 20, units = "in")


## ---- Supplementary figure -----

# A260/230 ratio is recovered after PCR and AmpPure cleaning

# create scale column
datnocon$scaleClean <- datnocon$Nanodrop_260230_cleanedPCR - 2.00

# plot
s <- ggdotplot(data = datnocon, x = "Method", y = "scaleClean",
          fill = "Method", binwidth = 0.8, size = 0.4, 
          add = c("median_iqr"), error.plot = "pointrange", #add.params = list(size = 1.5),
          xlab = FALSE, ylab = "Difference from 2.00 ideal", legend = "none") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", size = 1) +
  # add manual colors
  scale_fill_manual(values = methodcolswithspaces) +
  # add significance asterisk
  geom_text(label = "*", x = 2, y = 2.5, size = 7)

# save
ggsave(filename = "./Manuscript/SuppFigure1_Renamed.jpeg", plot = s,  device = "jpeg", dpi = 600, height = 4, width = 5.5, units = "in")

