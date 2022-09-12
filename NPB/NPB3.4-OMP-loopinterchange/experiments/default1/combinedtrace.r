library(tidyverse)
library(dplyr)
library(magrittr)

#!/usr/bin/env Rscript


df <- read_csv("times.ncc.csv", col_names=FALSE, col_types=cols()) %>%
        rename(Benchmark = X1, Duration = X2) %>%
        select(-contains("X")) %>%
        print
        
ncol1 <- character(length(df$Benchmark))
for(i in 1:length(df$Benchmark)) {
ncol1[i] <- 'nec'
}
df$Arch <- ncol1

print(df)


df2 <- read_csv("times.intel.csv", col_names=FALSE, col_types=cols()) %>%
        rename(Benchmark = X1, Duration = X2) %>%
        select(-contains("X")) %>%
        print
        
ncol2 <- character(length(df2$Benchmark))
for(i in 1:length(df2$Benchmark)) {
ncol2[i] <- 'intel'
}
df2$Arch <- ncol2

print(df2)


total <- merge(df,df2,by=c("Benchmark","Duration", "Arch"), all= TRUE)

print(total)

makespan = 5

vp <- total %>%
	group_by(Benchmark, Arch) %>%
	summarize(Duration.Percentage = sum(Duration)/makespan)
	
print(vp)
        
vp %>%
        ggplot(aes(x = Benchmark, y = Duration.Percentage, fill=Arch)) +
        theme_bw(base_size=16) +
        geom_bar(position="dodge", stat="identity") +
        scale_fill_brewer(palette = "Set1") +
        xlab("Benchmark") +
        ylab("Tempo [s]")
        
ggsave("necvshost.pdf", width=16, height=9)
