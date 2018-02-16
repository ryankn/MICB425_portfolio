## @knitr metadata26jan18
#Libraries
#install.packages("tidyverse")
library("tidyverse")

#Data Import
metadata <- read.table(file="DS_Friday/26Jan18/Saanich.metadata.txt", header=TRUE, row.names=1, sep="\t", na.strings="NAN")

#Exercise 1
OTU <- read.table(file="DS_Friday/26Jan18/Saanich.OTU.txt", header=TRUE, row.names=1, sep="\t", na.strings="NAN")

#Exercise 2
metadata %>% rownames_to_column('sample') %>% 
  filter(CH4_nM >= 100 & Temperature_C <= 10) %>% 
  column_to_rownames('sample') %>% 
  select(Depth_m,CH4_nM,Temperature_C)

newtable <- 
  metadata %>% rownames_to_column('sample') %>% 
  select(matches("nM|sample")) %>% 
  mutate(N2O_uM = N2O_nM/1000, Std_N2O_uM = Std_N2O_nM/1000, CH4_uM = CH4_nM/1000, Std_CH4_uM = Std_CH4_nM/1000) %>% 
  column_to_rownames('sample')

#Exercise 3
metadata %>% rownames_to_column('sample') %>% 
  select(matches("nM|sample")) %>% 
  mutate(N2O_uM = N2O_nM/1000, Std_N2O_uM = Std_N2O_nM/1000, CH4_uM = CH4_nM/1000, Std_CH4_uM = Std_CH4_nM/1000) %>% 
  column_to_rownames('sample') 