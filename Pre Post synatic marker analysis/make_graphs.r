library(ggplot2)
library(tidyverse)
library(infer)

folder = file.choose(new = FALSE)

# load the information about the series extracted from the lif file
# extract the condition by parsing the serie name
info <- read_csv(file.path(folder, "info.csv"),show_col_types = FALSE) %>%
  mutate(Condition = ifelse(str_detect(`Serie name`,".*wt.*"),"WT","dPDZ"))  %>%
  rename(`Image Name`=`Output  file`)

# reorder the groups
info$Condition <- ordered(info$Condition, levels=c("WT","dPDZ"))

# load the overlap table and make a join to get the conditions
overlap <- read_csv(file.path(folder, "result","overlaps.csv"),show_col_types = FALSE) %>%
  inner_join(info,by="Image Name")

# number of chained cluster by condition
overlap %>% ggplot(aes(x=`Condition`,y=`# chained`)) + geom_boxplot()

# ratio number of chained cluster vs number of bassoon
overlap %>% ggplot(aes(x=`Condition`,y=`# chained`/`bassoon`))+ geom_boxplot()

overlap %>% t_test(formula=`# chained`~`bassoon`)

# load the chained clusters table and add the condition
cl <- read_csv(file.path(folder, "result","clusters.csv"),show_col_types = FALSE)

# count the number of clusters per image from the cluster.csv file
cl %>%
  group_by(`Image Name`) %>%
  summarise(N = n()) %>%
  inner_join(info,by='Image Name') %>%
  ggplot(aes(x=`Condition`,y=`N`)) +
  geom_boxplot() +
  ggtitle("Number of 3 colours clusters per image")

