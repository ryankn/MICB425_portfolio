## @knitr DS_Assignment20180208
#Package Installation
#install.packages("tidyverse")
#source("https://bioconductor.org/biocLite.R")
#biocLite("phyloseq")

#Libraries
library("tidyverse")
library("phyloseq")

#Data Import
new_OTUs <- 
  read.table("DS_Friday/Assignment20180208/Saanich.OTU.new.txt",
             header = TRUE, sep = "\t", row.names = 1, na.strings = "NAN")
new_metadata <- 
  read.table("DS_Friday/Assignment20180208/Saanich.metadata.new.txt",
             header = TRUE, sep = "\t", row.names = 1, na.strings = "NAN")
load("DS_Friday/Assignment20180208/phyloseq_object.RData") 

#Exercise 1
ggplot(new_metadata, aes(x = CH4_nM, y = Depth_m)) +
  geom_point(color = "purple", shape = 17)

#Exercise 2
new_metadata %>%
  mutate(Temperature_F = Temperature_C * 9 / 5 + 32) %>%
  ggplot(aes(x = Temperature_F, y = Depth_m)) +
  geom_point()

#Exercise 3
physeq_percent = transform_sample_counts(physeq, function(x) 100 * x/sum(x))
plot_bar(physeq_percent, fill="Domain") + 
  geom_bar(aes(fill=Domain), stat="identity") +
  labs(x = "Sample depth", y = "Relative abundance (%)", title = "Domains from 10 to 200 m in Saanich Inlet")

#Exercise 4
new_metadata %>%
  select(matches("uM|depth"),-matches("Std"),-H2S_uM) %>%
  gather(key = "Nutrient", value = "Concentration", -Depth_m) %>%
  ggplot(., aes(x = Depth_m, y = Concentration)) +
  geom_point() +
  geom_line() +
  facet_wrap( ~ Nutrient, scales = "free") +
  theme(legend.position = "none") +
  labs(x = "Depth (m)", y = expression(paste("Concentration (", mu, "M)")))
