### Illustrate the exploring of the Colocalization.csv table using ggplot

library(ggplot2)
library(dplyr)

# get a filename
csvname <- file.choose(new = FALSE)

# load the file
df <- read.csv(csvname)

# Plot red/green colocalization
df %>% filter(Ch1=="red")  %>% filter(Ch2=="green") %>% ggplot(aes(x=Image,y=Pearson)) + geom_boxplot() + ggtitle("Ch1 / Ch2")

# Plot red/blue colocalization
df %>% filter(Ch1=="red")  %>% filter(Ch2=="blue") %>% ggplot(aes(x=Image,y=Pearson)) + geom_boxplot() + ggtitle("Ch1 / Ch3")
