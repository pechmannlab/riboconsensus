# Figure 1
# Extracted from notebook (16/7/19)
# Updated SP (14/9/19)
# Updated SP (27/9/19)

# Libraries:
library(jsonlite)
library(reshape2)
library(RColorBrewer)
#library(tidyverse)
library(ggplot2)
library(reshape2)


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



# SUPPL A ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
xv_corr_real <- read.table("data/figures/figure1/suppl_corr2mc_xv.txt", header=T, sep='\t')
xv_corr_rand <- read.table("data/figures/figure1/suppl_corr2randmc_xv.txt")
xv_corr_df <- data.frame(xv_corr_real$Correlation, xv_corr_rand[1:nrow(corr_real), 1])
colnames(xv_corr_df) <- c('real', 'random')
xv_corr_df <- melt(xv_corr_df)


svg(file = "figures/figure1/SUPPL_correlation_consensus.svg", height = 5, width = 6)

ggplot(data=xv_corr_df, aes(x=value,  fill=variable)) +
  geom_density(alpha=0.5) +
  labs(x = "Correlation", y = "Frequency") + 
  theme_classic() + 
  #scale_x_continuous(limits=c(-0.2, 0.6)) +
  scale_fill_manual(values = c("#E31A1C", 'grey60'),
                    name = "Data",
                    labels = c("consensus", "random")) +
  theme(
    axis.text = element_text(size=32),
    text = element_text(size=32),
    legend.justification = c(1.1, 1.1),
    legend.position = c(1, 1),
    legend.title = element_blank()
  )

dev.off()



# SUPPL B ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
peaks <- read.table("data/figures/figure1/suppl_peak_stats.txt", header=T)

svg(file = "figures/figure1/SUPPL_peaks.svg", height = 4, width = 4)

ggplot(peaks) + 
  geom_abline(aes(intercept = 0, slope=1), size=0.5, color="red", linetype="twodash" ) +
  geom_point(aes(x=peaks_mean/length, y=peaks_cons/length)) + 
  labs(x="Fraction individual peaks", y="Fraction consensus peaks") + 
  theme_classic() + 
  theme(
    axis.text = element_text(size=20),
    text = element_text(size=20)
  )

dev.off()



# SUPPL C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
biswas <- read.table("data/figures/figure1/suppl_biswas.txt", header=T)
biswas_m <- melt(biswas)

svg(file = "figures/figure1/SUPPL_biswas.svg", height = 4, width = 5)

ggplot(biswas_m) + 
  geom_boxplot(aes(x=variable, y=value)) + 
  scale_x_discrete(labels=c('th=0.001', 'th=0.01', 'th=0.00001', 'sd=1.5', 'sd=0.5', 'ca=0.9', 'ca=0.1')) + 
  labs(x="", y="Correlation") + 
  theme_classic() + 
  theme(
    text = element_text(size=20),
    axis.text.x = element_text(size=20, angle=90, hjust=1), 
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank()
  )

dev.off()


# SUPPL 1bA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
datasets <- read.table("data/figures/figure1/suppl_datasets.txt", header=T)
datasets_m <- melt(datasets, id='ORF')

svg(file = "figures/figure1/SUPPL_datasets.svg", height = 6, width = 8)

ggplot(datasets_m) + 
  geom_boxplot(aes(x=variable, y=value, fill=variable)) + 
  labs(x="", y="Correlation") + 
  scale_fill_viridis_d() + 
  theme_classic() + 
  theme(
    text = element_text(size=20),
    axis.text.x = element_text(size=20, angle=90, hjust=1), 
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(), 
    legend.position = 'none'
  )

dev.off()



# SUPPL 1bA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
lecanda <- read.table("data/figures/figure1/suppl_lecanda.txt", header=T)

svg(file = "figures/figure1/SUPPL_lecanda.svg", height = 6, width = 6)

ggplot(lecanda, aes(x=Correlation)) + 
  geom_histogram(bins=50, color="white", fill="turquoise") + 
  labs(x="Correlation", y="Count") + 
  theme_classic() + 
  theme(
    text = element_text(size=20),
    axis.text.x = element_text(size=20, angle=90, hjust=1), 
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(), 
    legend.position = 'none'
  )

dev.off()
