library(tidyverse)
library(dplyr)
library(magrittr)

#!/usr/bin/env Rscript
argas <- commandArgs(TRUE)
csv = paste ("times.",argas[1], ".csv",sep="")
pdf = paste ("times.",argas[1], ".pdf",sep="")

df <- read_csv("times.ncc.csv", col_names=FALSE, col_types=cols()) %>%
        rename(Benchmark = X1, Duration = X2) %>%
        select(-contains("X")) %>%
        print

ncol1 <- character(lenght(df$Benchmark))
for(i in 1:lenght(df$Benchmark)) {
ncol1[i] = nec
}
df$Arch <- ncol1

print(df)

df2 <- read_csv("times.intel.csv", col_names=FALSE, col_types=cols()) %>%
        rename(Benchmark = X1, Duration = X2) %>%
        select(-contains("X")) %>%
        print

ncol2 <- character(lenght(df2$Benchmark))
for(i in 1:lenght(df2$Benchmark)) {
ncol2[i] = intel
}
df2$Arch <- ncol2

print(df2)

total <- merge(df,df2,by="Benchmark")

print(total)

        
total %>%
        ggplot(aes(x = Benchmark, y = Duration, fill=Arch)) +
        theme_bw(base_size=16) +
        geom_bar(position="dodge", stat="identity") +
        scale_fill_brewer(palette = "Set1") +
        xlab("Benchmark") +
        ylab("Tempo [s]")
        
ggsave("necvshost.pdf", width=16, height=9)
