Nanodrop & QuBit Analysis
================
Emily Bean
1/11/2021

``` r
# load tidyverse package
if(!require(tidyverse)) {
  install.packages(tidyverse)
  require(tidyverse)
}

# load ggplot2
if(!require(ggplot2)) {
  install.packages(ggplot2)
  require(ggplot2)
}

# load ggpubr
if(!require(ggpubr)) {
  install.packages(ggpubr)
  require(ggpubr)
}

# load FSA for dunn test
if(!require(FSA)) {
  install.packages(FSA)
  require(FSA)
}

# set working directory
PATH <- "~/The Pennsylvania State University/Ganda, Erika - Shared-Trello-Projects/Kit Comparison/"
setwd(PATH)

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
         Nanodrop_260230_extrcDNA = replace(Nanodrop_260230_extrcDNA, which(Nanodrop_260230_extrcDNA < 0), 0))

# create dataframe without controls
datnocon <- dat %>% filter(!Replicate == "PosControl" & !Replicate == "NegControl")
```

## QuBit and Nanodrop Analysis

QuBit is used for DNA quantity and Nanodrop ratios are used for
extraction quality. All variables were found to have a non-normal
distribution (zero-skewed in all cases) that was not improved by
transformation, so non-parametric hypothesis tests were used. Multiple
comparisons are not necessary. The Dunn test was used as a post-hoc for
Kruksal Wallis.

In all variables, negative values are first replaced with 0.

*NOTE: As of 1/11/2021, the “Kit” numbers match Natalia’s “Method”
numbers*

### QuBit - Extracted DNA Quantity

There are differences between Kit 1 and 3, Kit 2 and 2, Kit 1 and 4, and
Kit 2 and 4.

``` r
# check for normality with a histogram
hist(dat$QuBit, main = "QuBit")
```

![](Nanodrop_QuBit_Analysis_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
# shapiro wilks test
shapiro.test(dat$QuBit) # not looking good
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  dat$QuBit
    ## W = 0.93188, p-value = 0.008006

``` r
# do Kruskal-Wallis
(mod <- kruskal.test(QuBit ~ Kit, data = dat))
```

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  QuBit by Kit
    ## Kruskal-Wallis chi-squared = 9.7807, df = 3, p-value = 0.02053

``` r
# KW post hoc is Dunn test
(dunnTest(QuBit ~ Kit, data = dat))
```

    ##    Comparison           Z    P.unadj     P.adj
    ## 1 Kit1 - Kit2 -0.10937300 0.91290665 1.0000000
    ## 2 Kit1 - Kit3 -2.30412454 0.02121565 0.1272939
    ## 3 Kit2 - Kit3 -2.19475154 0.02818143 0.1127257
    ## 4 Kit1 - Kit4 -2.22391767 0.02615398 0.1307699
    ## 5 Kit2 - Kit4 -2.11454467 0.03446877 0.1034063
    ## 6 Kit3 - Kit4  0.08020687 0.93607273 0.9360727

``` r
# visualize with a dot plot
ggdotplot(data = dat, x = "Kit", y = "QuBit",
          ylab = "ng/uL", 
          # add error bars
          add = "mean_se",
          error.plot = "errorbar", add.params = list(size = 2),
          # add title 
          title = "QuBit with mean + standard error")
```

![](Nanodrop_QuBit_Analysis_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

### Are there differences between individual chickens?

As a quality check to ensure that there are not certain chickens skewing
results that are otherwise attributed to kit differences, we will
analyze differences between chickens (called “Replicates” here).

There are no differences in QuBit DNA quantity between chickens
(Kruskal-Wallis p \> 0.05).

``` r
# do Kruskal-Wallis
(mod <- kruskal.test(QuBit ~ Replicate, data = datnocon)) # no differences
```

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  QuBit by Replicate
    ## Kruskal-Wallis chi-squared = 12.748, df = 9, p-value = 0.1743

``` r
# visualize with bar plot
ggbarplot(data = datnocon, x = "Replicate", y = "QuBit", group = "Kit", fill = "Kit",
          orientation = "horiz", title = "No differences in QuBit between chickens")
```

![](Nanodrop_QuBit_Analysis_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

## 260/280 ratio of extracted DNA

The 260/280 ratio is a measure of purity, as nucleic acids absorb at 260
nm and proteins absorb at 280 nm. The ideal value for this measurement
is 1.80 and anything less than \~1.60 is considered problematic for
downstream analysis.

Nanodrop was performed on extracted DNA before AmpPure cleaning.

To analyze, data was transformed to “difference from ideal” as `value
- 1.80`. To determine if each kit had higher or lower quality than
ideal, a non-parametric two-tailed Wilcoxon Rank Sum test was performed
for each kit to assess significant differences from 0.

This test reveals that Kit 2 and Kit 4 both differ from 0.

``` r
# make a scaled column
datnocon <- datnocon %>% 
  mutate(scale260280 = Nanodrop_260280_extrcDNA - 1.80)
# make a scaled column
dat <- dat %>% 
  mutate(scale260280 = Nanodrop_260280_extrcDNA - 1.80)

# set up loop for each kit
kit <- unique(datnocon$Kit)
output <- data.frame()

for(i in 1:length(kit)) {
  
  # two-tailed Wilcox rank sum test to test if the mean of each kit is different than 0
  w <- wilcox.test(datnocon$scale260280[datnocon$Kit == kit[[i]]], mu = 0, alternative = "two.sided",
                   exact = FALSE, conf.int = TRUE)
  # collect output
  output <- rbind(output,
                  data.frame(kit = kit[[i]],
                             Vstat = w$statistic,
                             pval = round(w$p.value, 3),
                             row.names = NULL))
  
}

output
```

    ##    kit Vstat  pval
    ## 1 Kit1  35.0 0.155
    ## 2 Kit2   6.5 0.037
    ## 3 Kit3  24.0 0.906
    ## 4 Kit4   7.0 0.041

### Differences between chickens

As above, we do a test to ensure that there are no differences between
chickens (replicates) that could skew the results. There are no
differences between chickens.

``` r
# do Kruskal-Wallis
(mod <- kruskal.test(scale260280 ~ Replicate, data = datnocon)) # no differences
```

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  scale260280 by Replicate
    ## Kruskal-Wallis chi-squared = 15.244, df = 9, p-value = 0.08445

``` r
# visualize with bar plot
ggbarplot(data = datnocon, x = "Replicate", y = "scale260280", group = "Kit", fill = "Kit",
          orientation = "horiz", xlab = "Difference from ideal 1.80")
```

![](Nanodrop_QuBit_Analysis_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

## 260/230 ratio of extracted DNA

The 260/230 ratio is a measure of contaminants other than proteins;
salts and phenols in particular absorb at 230 nm. The ideal ratio is
2.00 and low ratios are problematic for downstream analysis.

This data is interesting because almost all of the ratios are very low,
as shown by this plot:

``` r
# how well do the samples cluster around the ideal 2.00?
ggscatter(data = dat, x = "LabID", y = "Nanodrop_260230_extrcDNA") +
  geom_hline(yintercept = 2, color = "red") 
```

![](Nanodrop_QuBit_Analysis_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

The data was scaled similarly to 260/280 ratio except 2.00 was used as
the ‘ideal’ value.

The Wilcox rank rum test reveals that all 4 kits have significant
differences from 0.

``` r
# make a scaled column
datnocon <- datnocon %>% 
  mutate(scale260230 = Nanodrop_260230_extrcDNA - 2.0)
# make a scaled column
dat <- dat %>% 
  mutate(scale260230 = Nanodrop_260230_extrcDNA - 2.0)
# set up loop for each kit
kit <- unique(datnocon$Kit)
output <- data.frame()

for(i in 1:length(kit)) {
  
  # two-tailed Wilcox rank sum test to test if the mean of each kit is different than 0
  w <- wilcox.test(datnocon$scale260230[datnocon$Kit == kit[[i]]], mu = 0, alternative = "two.sided",
                   exact = FALSE, conf.int = TRUE)
  # collect output
  output <- rbind(output,
                  data.frame(kit = kit[[i]],
                             stat = w$statistic,
                             pval = round(w$p.value, 5),
                             row.names = NULL))
  
}

output
```

    ##    kit stat    pval
    ## 1 Kit1    0 0.00589
    ## 2 Kit2    0 0.00592
    ## 3 Kit3    0 0.00592
    ## 4 Kit4    0 0.00583

### Are there differences between individual chickens?

There are no differences in 260/230 ratio between chickens.

``` r
# another quality check; are there differences between our replicates?
# do Kruskal-Wallis
(mod <- kruskal.test(scale260230 ~ Replicate, data = datnocon)) # no differences
```

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  scale260230 by Replicate
    ## Kruskal-Wallis chi-squared = 3.7194, df = 9, p-value = 0.9289

``` r
# visualize with bar plot
ggbarplot(data = datnocon, x = "Replicate", y = "scale260230", group = "Kit", fill = "Kit",
          orientation = "horiz", xlab = "Difference from ideal")
```

![](Nanodrop_QuBit_Analysis_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

### Supplementary material: A260/230 is recovered after PCR and cleaning

Nanodrop was performed on the samples after PCR and AmpPure cleaning
procedures. The same scaling and statistical test was applied to this
data.

The purity of samples in Kit 1, 2, and 4 is recovered (p \> 0.05)
however, Kit 2 still has a significant difference from 0 after cleaning
(p = 0.03).

``` r
# make a scaled column
datnocon <- datnocon %>% 
  mutate(cleanscale = Nanodrop_260230_cleanedPCR - 2.0)
# make a scaled column
dat <- dat %>% 
  mutate(cleanscale = Nanodrop_260230_cleanedPCR - 2.0)
# set up loop for each kit
kit <- unique(datnocon$Kit)
output <- data.frame()

for(i in 1:length(kit)) {
  
  # two-tailed Wilcox rank sum test to test if the mean of each kit is different than 0
  w <- wilcox.test(datnocon$cleanscale[datnocon$Kit == kit[[i]]], mu = 0, alternative = "two.sided",
                   exact = FALSE, conf.int = TRUE)
  # collect output
  output <- rbind(output,
                  data.frame(kit = kit[[i]],
                             stat = w$statistic,
                             pval = round(w$p.value, 5),
                             row.names = NULL))
  
}

output
```

    ##    kit stat    pval
    ## 1 Kit1   44 0.10292
    ## 2 Kit2    6 0.03120
    ## 3 Kit3   43 0.12603
    ## 4 Kit4   44 0.10292
