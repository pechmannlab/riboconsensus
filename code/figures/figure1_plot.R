# Figure 1
# Extracted from notebook (16/7/19)
# Updated SP (14/9/19)
# Updated SP (27/9/19)

# Libraries:
library(jsonlite)
library(reshape2)
library(RColorBrewer)
library(tidyverse)
library(ggplot2)



# A ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
YKL143W_profiles <- read.table("data/figures/figure1/YKL143W_profiles.txt")
YKL143W_profiles <- as.data.frame( t(YKL143W_profiles) )
YKL143W_profiles$position <- 1:nrow(YKL143W_profiles)
profiles_df = melt(YKL143W_profiles, id=c("position"))


svg(file = "figures/figure1/A_YKL143W_profiles.svg", height =18, width = 12)

p0 <- ggplot(data=profiles_df, aes(x=position, y=value)) +
  geom_area(fill='#1F78B4') +
  theme_classic() +
  labs(x = "Codon position", y = "Ribosome density") +
  theme(
    axis.text.y = element_blank(),
    text = element_text(size=40)
  )

p0 + facet_grid(rows = vars(variable), scale = "free_y") +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        axis.line.y = element_line() )

dev.off()


# B ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
YKL143W_consensus <- read.table('data/figures/figure1/YKL143W_consensus.txt')
consensus_df <- as.data.frame(YKL143W_consensus)
colnames(consensus_df) <- 'consensus'
consensus_df$position <- 1:nrow(consensus_df)


svg(file = "figures/figure1/B_YKL143W_consensus.svg", width=12, height=3)

ggplot(data=consensus_df, aes(x=position, y=consensus)) +
  geom_area(fill='#E31A1C') +
  #scale_color_grey(start=0.8, end=0.2) +
  theme_classic() +
  labs(y = "Consensus", x = "Codon position") +
  theme( text = element_text(size=32),
         axis.text.x = element_text(size=32),
         axis.text.y = element_text(size=32)
         )

dev.off()




# C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
corr_real <- read.table("data/figures/figure1/corr2mc.txt", header=T, sep='\t')
corr_rand <- read.table("data/figures/figure1/corr2randmc.txt")
corr_df <- data.frame(corr_real$Correlation, corr_rand[1:nrow(corr_real), 1])
colnames(corr_df) <- c('real', 'random')
corr_df <- melt(corr_df)


svg(file = "figures/figure1/C_correlation_consensus.svg", height = 5, width = 6)

ggplot(data=corr_df, aes(x=value,  fill=variable)) +
  geom_density(alpha=0.5) +
  labs(x = "Correlation", y = "Frequency") + 
  theme_classic() + 
  scale_fill_manual(values = c("#E31A1C", 'grey60'),
                    name = "Data",
                    labels = c("consensus", "random")) +
  theme(
      axis.text = element_text(size=32),
      text = element_text(size=32),
      legend.justification = c(1.1, 1.1),
      legend.position = c(1, 1)
      )

dev.off()



# D ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
all_consensus <- as.data.frame( read.table("data/figures/figure1/all_consensus.txt") )
colnames(all_consensus) <- 'consensus'


svg(file = "figures/figure1/D_histo_consensus.svg", height = 5, width = 6)

ggplot(data = all_consensus) +
  geom_histogram(aes(x = consensus),
                 binwidth = 0.025,
                 fill="#E31A1C",
                 color='black',
                 alpha=0.5
  ) +
  labs(x = "Consensus score",
       y = "Frequency"
  ) +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size=32),
    text = element_text(size=32)
  )

dev.off()


# E ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dataqc <- read.table("data/figures/figure1/data_qc.txt", header=T, sep='\t')

svg(file = "figures/figure1/E_coverage.svg", height = 5, width = 6)

ggplot(dataqc, aes(x=log(EL), y=corr)) + 
  geom_point(aes(size=cons90), show.legend = FALSE, color="#E31A1C", alpha=0.5) + 
  theme_classic()  +
    theme(
    panel.grid = element_blank(),
    axis.text = element_text(size=32),
    text = element_text(size=32)
  ) + 
  labs(x="Average coverage (log)", y="Correlation") 

dev.off()