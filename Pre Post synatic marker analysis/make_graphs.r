library(ggplot2)
library(tidyverse)
library(infer)

# load the information about the series extracted from the lif file
# extract the condition by parsing the serie name
info <- read_csv("/media/micro-nas-1/Greger_lab/imogen/202622/info.csv",show_col_types = FALSE) %>%
  mutate(Condition = ifelse(str_detect(`Serie name`,".*wt.*"),"WT","dPDZ"))  %>%
  rename(`Image Name`=`Output  file`)



# load the overlap table and make a join to get the conditions
overlap <- read_csv("/media/micro-nas-1/Greger_lab/imogen/202622/result/overlaps.csv",show_col_types = FALSE) %>%
  na_if("NaN") %>%
  inner_join(info,by="Image Name")

# number of chained cluster by condition
overlap %>% ggplot(aes(x=`Condition`,y=`# chained`))+ geom_point()

# ratio number of chained cluster vs number of bassoon
overlap %>% ggplot(aes(x=`Condition`,y=`# chained`/`bassoon`))+ geom_boxplot()



# load the chained clusters table and add the condition
cl <- read_csv("/media/micro-nas-1/Greger_lab/imogen/202622/result/clusters.csv",show_col_types = FALSE)

%>%


  na_if("NaN") %>%
  inner_join(info,by="Image Name")

# count the number of clusters per image
cl %>%
  group_by(`Image Name`) %>%
  filter(`Distance [glua2:bassoon]` < 0.25) %>%
  summarise(N = n()) %>%
  inner_join(info,by='Image Name')

%>%
  ggplot(aes(x=`Condition`,y=`n`)) +
  geom_boxplot() +
  ggtitle("Number of 3 colours clusters per image")

cl %>%
  ggplot(aes(x=`Condition`,y=`Distance [psd95:glua2]`)) +
  geom_violin()
