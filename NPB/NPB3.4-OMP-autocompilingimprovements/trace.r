library(tidyverse)
library(dplyr)
library(magrittr)

#!/usr/bin/env Rscript
argas <- commandArgs(TRUE)
csv1 = paste ("times.",argas[1], ".csv",sep="")
pdf1 = paste ("times.",argas[1], ".pdf",sep="")
csv2 = paste ("cache.flops.",argas[1], ".old.csv",sep="")
pdf2 = paste ("cache.",argas[1], ".pdf",sep="")
pdf3 = paste ("flops.",argas[1], ".pdf",sep="")


df <- read_csv(csv1, col_names=FALSE, col_types=cols()) %>%
        rename(Benchmark = X1, Duration = X2) %>%
        select(-contains("X")) %>%
        print
        
df %>%
        ggplot(aes(x = Benchmark, y = Duration, Age = "original")) +
        theme_bw(base_size=16) +
        geom_bar(position="dodge", stat="identity") +
        scale_fill_brewer(palette = "Set1") +
        xlab("Benchmark") +
        ylab("Tempo [s]")
        
ggsave(pdf1, width=16, height=9)

df2 <- read_csv(csv2, col_names=FALSE, col_types=cols()) %>%
        rename(Benchmark = X1, Cache = X2) %>%
        select(-contains("X")) %>%
        print
        
df2 %>%
        ggplot(aes(x = Benchmark, y = Cache)) +
        theme_bw(base_size=16) +
        geom_bar(position="dodge", stat="identity") +
        scale_fill_brewer(palette = "Set1") +
        xlab("Benchmark") +
        ylab("Taxa de acerto da cache [%]")
        
ggsave(pdf2, width=16, height=9)

df3 <- read_csv(csv2, col_names=FALSE, col_types=cols()) %>%
        rename(Benchmark = X1, Flops = X3) %>%
        select(-contains("X")) %>%
        print
        
df3 %>%
        ggplot(aes(x = Benchmark, y = Flops)) +
        theme_bw(base_size=16) +
        geom_bar(position="dodge", stat="identity") +
        scale_fill_brewer(palette = "Set1") +
        xlab("Benchmark") +
        ylab("MFlops")
        
ggsave(pdf3, width=16, height=9)
