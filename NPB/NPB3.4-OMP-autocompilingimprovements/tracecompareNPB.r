library(tidyverse)
library(dplyr)
library(magrittr)

#!/usr/bin/env Rscript
argas <- commandArgs(TRUE)
csv1 = paste ("times.ncc.old.csv",sep="")
csv2 = paste ("times.ncc.new.csv",sep="")
csv3 = paste ("cache.flops.ncc.old.csv",sep="")
csv4 = paste ("cache.flops.ncc.new.csv",sep="")
pdf2 = paste ("cache.newvsold.pdf",sep="")
pdf3 = paste ("flops.newvsold.pdf",sep="")
pdf4 = paste ("taxavet.newvsold.pdf",sep="")
pdf = paste ("comparisonnewvsog.png",sep="")

df <- read_csv(csv1, col_names=FALSE, col_types=cols()) %>%
        rename(Benchmark = X1, Duration = X2) %>%
        select(-contains("X")) %>%
        print

ncol1 <- character(length(df$Benchmark))
for(i in 1:length(df$Benchmark)) {
ncol1[i] <- 'original'
}
df$Age <- ncol1

df2 <- read_csv(csv2, col_names=FALSE, col_types=cols()) %>%
        rename(Benchmark = X1, Duration = X2) %>%
        select(-contains("X")) %>%
        print

ncol2 <- character(length(df2$Benchmark))
for(i in 1:length(df2$Benchmark)) {
ncol2[i] <- 'otimizada'
}
df2$Age <- ncol2

df3 <- read_csv(csv3, col_names=FALSE, col_types=cols()) %>%
        rename(Benchmark = X1, Cache = X2, Flops = X3, Vet = X4, Secache = X5, Seflops = X6, Sevet = X7) %>%
        select(-contains("X")) %>%
        print

ncol3 <- character(length(df3$Benchmark))
for(i in 1:length(df3$Benchmark)) {
ncol3[i] <- 'original'
}
df3$Age <- ncol3

df4 <- read_csv(csv4, col_names=FALSE, col_types=cols()) %>%
        rename(Benchmark = X1, Cache = X2, Flops = X3, Vet = X4, Secache = X5, Seflops = X6, Sevet = X7) %>%
        select(-contains("X")) %>%
        print

ncol4 <- character(length(df4$Benchmark))
for(i in 1:length(df4$Benchmark)) {
ncol4[i] <- 'otimizada'
}
df4$Age <- ncol4

total <- merge(df,df2,by=c("Benchmark","Duration", "Age"), all= TRUE)
total3 <- merge(df3,df4,by=c("Benchmark","Cache", "Secache", "Age"), all= TRUE)
total4 <- merge(df3,df4,by=c("Benchmark","Flops", "Seflops", "Age"), all= TRUE)
total5 <- merge(df3,df4,by=c("Benchmark","Vet", "Sevet", "Age"), all= TRUE)

print(total)
print(total3)
print(total4)
        

total4 %>%
        ggplot(aes(x = Benchmark, y = Flops, fill = Age)) +
        theme_bw(base_size=16) +
        labs(fill = "Versão") +
        geom_bar(position="dodge", stat="identity") +
        geom_errorbar(aes(x=Benchmark, ymin=Flops-Seflops, ymax=Flops+Seflops), width=0.5, colour="black", alpha=0.8, size=0.5, position = position_dodge(0.9)) +
        scale_fill_brewer(palette = "Set1") +
        xlab("Benchmark") +
        ylab("MFLOPS") +
        theme_bw() +
        theme(axis.text=element_text(size=20), axis.title=element_text(size=32), legend.title = element_text(size = 28), legend.text = element_text(size = 24))
        
ggsave(pdf3, width=16, height=9)

total5 %>%
        ggplot(aes(x = Benchmark, y = Vet, fill = Age)) +
        theme_bw(base_size=16) +
        labs(fill = "Versão") +
        geom_bar(position="dodge", stat="identity") +
        geom_errorbar(aes(x=Benchmark, ymin=Vet-Sevet, ymax=Vet+Sevet), width=0.5, colour="black", alpha=0.8, size=0.5, position = position_dodge(0.9)) +
        scale_fill_brewer(palette = "Set1") +
        xlab("Benchmark") +
        ylab("Taxa de operações vetoriais [%]") +
        theme_bw() +
        theme(axis.text=element_text(size=20), axis.title=element_text(size=32), legend.title = element_text(size = 28), legend.text = element_text(size = 24))
        
ggsave(pdf4, width=16, height=9)