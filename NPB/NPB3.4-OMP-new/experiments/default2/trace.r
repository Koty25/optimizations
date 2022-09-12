library(tidyverse)
library(dplyr)
library(magrittr)

#!/usr/bin/env Rscript
argas <- commandArgs(TRUE)
csv = paste ("times.",argas[1], ".csv",sep="")
pdf = paste ("times.",argas[1], ".pdf",sep="")

df <- read_csv(csv, col_names=FALSE, col_types=cols()) %>%
        rename(Benchmark = X1, Duration = X2) %>%
        select(-contains("X")) %>%
        print
        
df %>%
        ggplot(aes(x = Benchmark, y = Duration)) +
        theme_bw(base_size=16) +
        geom_bar(position="dodge", stat="identity") +
        scale_fill_brewer(palette = "Set1") +
        xlab("Benchmark") +
        ylab("Tempo [s]")
        
ggsave(pdf, width=16, height=9)
