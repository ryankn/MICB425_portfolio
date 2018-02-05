#Libraries
#install.packages("tidyverse")
library("tidyverse")

#Data Import
metadata <- read.table(file="DS_Friday/26Jan18/Saanich.metadata.txt", header=TRUE, row.names=1, sep="\t", na.strings="NAN")

OTU <- read.table(file="DS_Friday/26Jan18/Saanich.OTU.txt", header=TRUE, row.names=1, sep="\t", na.strings="NAN")

#Doin a thing
metadata %>% select(O2_uM)

metadata %>% select(matches("O2|oxygen"))

# select also works with starts_with, ends_with, contains

metadata %>% filter(O2_uM == 0) %>% select(Depth_m)

metadata %>% rownames_to_column('sample') %>% 
  filter(CH4_nM >= 100 & Temperature_C <= 10) %>% 
  column_to_rownames('sample') %>% 
  select(Depth_m,CH4_nM,Temperature_C)

newtable <- 
  metadata %>% rownames_to_column('sample') %>% 
  select(matches("nM|sample")) %>% 
  mutate(N2O_uM = N2O_nM/1000, Std_N2O_uM = Std_N2O_nM/1000, CH4_uM = CH4_nM/1000, Std_CH4_uM = Std_CH4_nM/1000) %>% 
  column_to_rownames('sample')

metadata %>% rownames_to_column('sample') %>% 
  select(matches("nM|sample")) %>% 
  mutate(N2O_uM = N2O_nM/1000, Std_N2O_uM = Std_N2O_nM/1000, CH4_uM = CH4_nM/1000, Std_CH4_uM = Std_CH4_nM/1000) %>% 
  column_to_rownames('sample') %>% 
  View

