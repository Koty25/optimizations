#!/usr/bin/env Rscript

library(tidyverse)
library(dplyr)
library(magrittr)
library(RColorBrewer)
library(colorRamps)

argas <- commandArgs(TRUE)
csv = paste ("test.",argas[1], ".csv",sep="")
pdf = paste ("test.",argas[1],".", argas[2], ".pdf",sep="")

df <- read_csv(csv, col_names=FALSE, col_types=cols()) %>%
        rename(Benchmark = X1, region = X2, Duration = X3, Percentage = X4) %>%
        select(-contains("X")) %>%
        #	as.numeric(Duration)
        print
options(digits=22)
df$Duration <- as.numeric(df$Duration)
if (argas[2] != 'all'){
	df <- df %>% filter(Benchmark == argas[2])
}


print(df)

sapply(df, mode)
sapply(df, class)

colourCount <- length(unique(df$region)) # number of levels
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
        
df %>%
        ggplot(aes(x = Benchmark, y = Duration, fill = region)) +
        theme_bw(base_size=16) +
        geom_bar(position="dodge", stat="identity") +
        xlab("Benchmark") +
        ylab("Tempo [s]")
        
ggsave(pdf, width=16, height=9)