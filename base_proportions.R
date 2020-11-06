#######################################################################################
###						           ><(((º> 	base_proportions.R  <º)))><  	    			      		###
###	Script to calculate and plot the base proportions from fq.gz files' base counts ###
#######################################################################################

### Contact: Eric Garcia, e1garcia@odu.edu  <º)))><
### This script is optimized to calculate base pair proportions, position by position, of DNA fragments coming from single digest RADseq, paired-End sequencing data.

# First, prepare an input tab-delimited file(s) with individuals/files in rows and data as columns. 
# First column being the file names, then counts of each of the bases (in this order: A,C,G,T,N), position, which end of read, and finally, read direction. 
# You can use the "base_calculator.sh" script in GitHub to do this for you https://github.com/philippinespire/scripts_and_code

# Thus, your input file(s) should have 9 columns. Head example file:
# ```
# file	            A	    C	    G   	T   	N positions frag_end    readDir
# PIRE*.F.fq.gz	    38347	40248	40488	35540	3 1-1       restrct_end F
# ```

# This example analyzes untrimmed forward (F) and reverse (R) reads and trimmed forward (R1) and reverse (R2) reads from Siganus spinus (Ssp) baited (b) and unbaited (unb) data.

# Place the current R script file in the same directory with your input file(s)
# This R script can further modify this file(s), calculate means, graph variation in base proportions, etc

# Clean your global environment
rm(list = ls())

# Set your working directory to the same where your script is
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

# Load packages (install.packages() if necessary)
library(broom)
library(glue)
library(janitor)
library(lubridate)
library(magrittr)
library(readxl)
library(tidyverse)
library(zoo)

# Import tsv datasets (input files). This example:
Ssp_b_F    <- read_tsv("allbase_prop_Ssp_baited_78.F_win1_st1_en140_theONE.tsv")
Ssp_b_R    <- read_tsv("allbase_prop_Ssp_baited_78.R_win1_st1_en150_theONE.tsv")
Ssp_b_R1   <- read_tsv("allbase_prop_Ssp_baited_R1_win1_st1_en140_theONE.tsv")
Ssp_b_R2   <- read_tsv("allbase_prop_Ssp_baited_R2_win1_st1_en150_theONE.tsv")
Ssp_unb_F  <- read_tsv("allbase_prop_Ssp_unbaited_F_win1_st1_en140_theONE.tsv")
Ssp_unb_R  <- read_tsv("allbase_prop_Ssp_unbaited_R_win1_st1_en150_theONE.tsv")
Ssp_unb_R2 <- read_tsv("allbase_prop_Ssp_unbaited_R1_win1_st1_en140_theONE.tsv")
Ssp_unb_R1 <- read_tsv("allbase_prop_Ssp_unbaited_R2_win1_st1_en150_theONE.tsv")

# Check that your files were imported correctly
# In this occasion the "F" from the readDir column was getting converted into <logical> "FALSE" values
# Thus, I will simply remove readDir column and recovered the read type from the file names later
Ssp_unb_F  <- Ssp_unb_F [,c(1:8)]
Ssp_unb_R  <- Ssp_unb_R [,c(1:8)]
Ssp_unb_R1 <- Ssp_unb_R1[,c(1:8)]
Ssp_unb_R2 <- Ssp_unb_R2[,c(1:8)]
Ssp_b_F    <- Ssp_b_F   [,c(1:8)]
Ssp_b_R    <- Ssp_b_R   [,c(1:8)]
Ssp_b_R1   <- Ssp_b_R1  [,c(1:8)]
Ssp_b_R2   <- Ssp_b_R2  [,c(1:8)]

### All forward read are supposed to start with the restriction motif and all reverse reads are supposed to exclude it. Yet, this often is not the case.
# Depending on your analysis, you can use this code to limit your data to only forward reads that do have the restriction site and reverse reads that do not have it:
#Ssp_b_R1 <- Ssp_unb_R1 %>%
#  filter(frag_end=="restrct_end")
#Ssp_b_R2 <- Ssp_b_R2 %>%
#  filter(frag_end=="sheared_end")

### The "base_calculator.sh" script provides important information in the file's name
# We will use the names to split this information into variable (columns) but first:
# check if all files have the same name schemes across datasets, i.e. names provide the same info in the same order.
Ssp_unb_F [1,1]
Ssp_unb_R [1,1]
Ssp_unb_R1[1,1]
Ssp_unb_R2[1,1]
Ssp_b_F   [1,1]
Ssp_b_R   [1,1]
Ssp_b_R1  [1,1]
Ssp_b_R2  [1,1]

# In this example, there are 2 different name schemes,
# one for the baited "20200331_PIRE-Ssp-AAtu_005-Capture_CKDL200150313-1a-DY0088-AK1681_H7T2LBBXX_L78.F.fq.gz"
# and another for the unbaited "PIRE2019-Ssp-A-Atu_010-Plate1Pool5Seq1-1E-L4.R1.fq.gz"
# Thus, create one table for all baited and one for all unbaited files
tbl_unb <- bind_rows(Ssp_unb_F, Ssp_unb_R, Ssp_unb_R1, Ssp_unb_R2)
tbl_b <- bind_rows( Ssp_b_F, Ssp_b_R, Ssp_b_R1, Ssp_b_R2)

# Now, we are ready clean names, omit missing data if any,
# and then, use the file names to "separate" valuable info into variables

tbl_b <- tbl_b %>% 
  clean_names() %>%       
  na.omit()  %>%                                                                                                                   
  separate(col=file,into=c("XProject","sp","time_period_individual","xlib_method","x1","x2","extension"),sep="-", remove=FALSE) %>%
  separate(col=time_period_individual,into=c("time_period","individual"),sep="_", remove=TRUE) %>%
  separate(col=extension,into=c("x3","read","x4","x5"),sep="\\.", remove=TRUE) %>%
  separate(col=positions,into=c("st_pos","en_pos"),sep="-", remove=TRUE) %>%
  select(-starts_with("x")) 

tbl_unb <- tbl_unb %>% 
  clean_names() %>%       
  na.omit()  %>%                                                                                                                   
  separate(col=file,into=c("XProject","sp","time_period","indiv","xlib_method","x1","extension"),sep="-", remove=FALSE) %>%
  separate(col=indiv,into=c("popu","individual"),sep="_", remove=TRUE) %>%
  separate(col=extension,into=c("x3","read","x5","x6"),sep="\\.", remove=TRUE) %>%
  separate(col=positions,into=c("st_pos","en_pos"),sep="-", remove=TRUE) %>%
  select(-starts_with("x")) 


# For this analysis, I only want "Albatross" or "Contemporary" in "time_period" columns 
# (instead of AAtu, CGub, A, and C)
tbl_b$time_period <- gsub("AAtu", "Albatross", tbl_b$time_period)
tbl_b$time_period <- gsub("CGub", "Contemporary", tbl_b$time_period)
tbl_unb$time_period <- gsub("A", "Albatross", tbl_unb$time_period)
tbl_unb$time_period <- gsub("C", "Contemporary", tbl_unb$time_period)

# Digest data into a long table and calculate the base proportions per position for each read type and time period
tbl_b_long <- tbl_b %>% 
  group_by(time_period, read, st_pos) %>%
  summarise(a = sum(a), c = sum(c), g = sum(g), t = sum(t), n = sum(n)) %>%
  pivot_longer(col=a:n, names_to="base", values_to="base_count") %>%
  group_by(time_period, read, st_pos) %>%
  mutate(proportion = base_count*100/sum(base_count))

tbl_unb_long <- tbl_unb %>% 
  group_by(time_period, read, st_pos) %>%
  summarise(a = sum(a), c = sum(c), g = sum(g), t = sum(t), n = sum(n)) %>%
  pivot_longer(col=a:n, names_to="base", values_to="base_count") %>%
  group_by(time_period, read, st_pos) %>%
  mutate(proportion = base_count*100/sum(base_count))



### Now you have two tables ready to be plotted.
# Repeat the following code for each table
# Instead of duplicating the code, 
#you can simply switch back and forward between "b" and "unb" to specify the table to be analyzed


# For cleaner looking graphs, I like capitalizing the base names
tbl_unb_long$base <- gsub("a", "A", tbl_unb_long$base)
tbl_unb_long$base <- gsub("c", "C", tbl_unb_long$base)
tbl_unb_long$base <- gsub("g", "G", tbl_unb_long$base)
tbl_unb_long$base <- gsub("t", "T", tbl_unb_long$base)
tbl_unb_long$base <- gsub("n", "N", tbl_unb_long$base)

# In order to plot from position 1 to 150, set st_pos as a factor with a sequence equal to the max read length(150) as the levels.
# Then, to fit a model to your data, coerce st_pos as numeric since st_pos is ploted in the x-axis
tbl_unb_long$st_pos <- factor(tbl_unb_long$st_pos, levels = c(1:150))
tbl_unb_long$st_pos <- as.numeric(tbl_unb_long$st_pos)

## You can set the order in which you wish to plot the bases with:
tbl_b_long$base <- factor(tbl_b_long$base, levels = c("A", "C", "G", "T", "N"))

# Plot the base proportions of unbaited data by type of read
tbl_b_long %>% 
#  filter(read == c("R1","R2")) %>%          # if you only want to filter by read (trimmed R1 and R2 in this example)
  filter(base == c("A","C","G","T")) %>%     # if you want to exclude Ns
  ggplot(aes(x=st_pos, y=proportion, color=time_period)) +
  labs(title ="Base proportion per position in Ssp unbaited untrimmed (F-R) and trimmed (R1-R2) files" , x = "Position", y = "Base Proportion") +
#  geom_point() +                            # if you would like to plot data points
  geom_smooth(method = "lm", se = TRUE, fullrange= TRUE) +
  theme(axis.text.x =element_text(size=8, angle=90, hjust=1)) +
  facet_grid(base~read) 

# You can also apply conditional filtering, or filtering by a combination of diff. values from diff. variables, with "filter(xor"
# Example, If I only care for the restriction end in the forward trimmed reads and sheared end in reverse,  I can limit my data this way:
# tbl_b_long %>% 
#  filter(xor((read == "R1" & frag_end == "restrct_end"), (read == "R2" & frag_end =="sheared_end"))) %>% 

### Repeat plot with baited data.
# Plots are highly customizable analytically and aesthetically  
# (you can limit to any variable or value. Fewer read types, bases, or positions, etc)  
# (you can show data points, fit a line, or both, etc)


###### Calculate mean and stdev of base-proportions in any specified group of files

# Ex 1: time_period, read, and base
means_b <- tbl_b_long %>% 
  filter(base == c("A","C","G","T")) %>%
  group_by(time_period, read, base) %>%
  summarise(mean_base_prop = mean(proportion)) 

sd_b <- tbl_b_long %>% 
  filter(base == c("A","C","G","T")) %>%
  group_by(time_period, read, base) %>%
  summarise(sd_means = sd(proportion)) 

stats_b <- cbind(means_b, sd_means=sd_b$sd_means)

# In only a subset of the data. Untrimmed reverse reads in this case "R"
means_R <- tbl_b_long %>% 
  filter(base == c("A","C","G","T")) %>%
  filter(read == "R") %>%
  group_by(time_period, read, base) %>%
  summarise(mean_base_prop = mean(proportion))

sd_R <- tbl_b_long %>% 
  filter(base == c("A","C","G","T")) %>%
  filter(read == "R") %>%
  group_by(time_period, read, base) %>%
  summarise(sd_means = sd(proportion))

stats_R <- cbind(means_R, sd_means=sd_R$sd_means)


# Plot the means with error bars
mean_b_plot <- ggplot(stats_b, aes(x=base, y=mean_base_prop, color=time_period)) + 
  geom_bar(stat="identity", fill="white", position=position_dodge()) +
  geom_errorbar(aes(ymin=mean_base_prop-sd_means, ymax=mean_base_prop+sd_means), width=.2,
                position=position_dodge(.9)) +
  ylim(0,50) +
  labs(title ="Mean base proportion per Read type and Time period", x = "Base", y = "Mean Proportion + stdev") +
  scale_color_manual(values=c("black", "blue")) +
  facet_grid(.~read)

mean_b_plot

