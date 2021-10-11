##########################################################################
## Make figures for the analysis of spread amd intensity
##
## The result file has extea columns Run and Experiment
##
## Jerome Boulanger
##########################################################################

rm(list = ls())
options(crayon.enabled = FALSE) # to work in emacs ESS
library(tidyverse)
library(rstatix)
library(ggpubr) # for stat_pvalue_manual
library(gridExtra)
library(grid)

#######################
### Load result file
#######################
filename = "./results.csv"
df <- read.csv(filename)
## Make sure that the levels are in the specified order
df$Condition <- factor(df$Condition, level=c("Control","KO"))

#######################
### Spread
#######################
# pairwise t test
pwc1 <- df %>%
    group_by(Run, Experiment) %>%
    t_test(Object.Spread ~ Condition) %>%
    add_significance()  %>%
    add_xy_position(x = "Condition")

## reset the xmin and xmax to this works with the scales="free" option in facet_grid
pwc1$xmin=1
pwc1$xmax=2

## tweak the run labels
run.labs = c("X","Y","Z")
names(run.labs) <- c(1,2,3)

## make the graphs
ggplot(df,aes(x=Condition, y=Object.Spread)) + geom_boxplot() + expand_limits(y = 0) + scale_y_continuous(name = paste("Spread [","\u03bc","m]",sep="")) + scale_x_discrete(guide = guide_axis(angle = 90))  + facet_grid(Run~Experiment,scales="free",labeller=labeller(Experiment=label_both,Run=run.labs)) + stat_pvalue_manual(pwc1) + theme_classic()

ggsave("spread.pdf")

########################
### Intensity
########################

## Intensities channel name
channels = c("ROI.Mean.DAPI","ROI.Mean.GFP","ROI.Mean.AF555")
channels.labs = lapply(channels, function (x){ gsub("\\."," ",gsub("ROI\\.Mean\\.","",x))})

# to label the facets
mylabeller <- function (k) {
    print(k)
    gsub("\\."," ",gsub("ROI\\.Mean\\.","",channels[strtoi(k)]))
}

# Combine all Mean Intensity graphs in a grid Experiment x Channel
intensity <- data.frame()
for (k in 1:length(channels)) {
    tmp <- df %>% filter(Run==1) %>%
        filter(get(channels[k])>0) %>%
        mutate(Intensity=get(channels[k])/max(get(channels[k])),Channel=k)
    intensity <- rbind(intensity,tmp)
}

pwc2 <- intensity %>%
    group_by(Experiment,Channel) %>%
    t_test(Intensity ~ Condition) %>%
    add_significance()  %>%
    add_xy_position(x = "Condition")
# we have only 2 conditions per experiment
pwc2$xmin=1
pwc2$xmax=2

ggplot(intensity, aes(x=Condition, y=Intensity)) + geom_boxplot() + expand_limits(y = 0) + facet_grid(Channel~Experiment,scales="free",space="free",labeller=labeller(Channel=mylabeller,Experiment=label_both)) + stat_pvalue_manual(pwc2) + theme_classic()

ggsave("intensity.pdf")
