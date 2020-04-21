library(ggplot2)
library(reshape2)




# A ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
theta = as.data.frame(read.csv("data/figures/figure3/theta.txt", header=T, sep='\t'))
theta = data.frame(theta$p3_5, theta$p3_10, theta$p3_20, theta$p3_80, theta$p3_90, theta$p3_95)
theta_m = melt(theta)


svg(file = "figures/figure3/A_percentiles.svg", height = 5, width = 6)

ggplot(data = theta_m, aes(x=value, fill=variable)) +
  geom_density(alpha=0.5) +
  labs(
    x = "Expected consensus score",
    y = "Frequency"
  ) +
  theme_classic() +
  scale_fill_manual(values = c("darkblue", "blue", "#A6CEE3", "orange", "gold", "yellow"),
                    name = "Percentile",
                    labels = c("5%", "10%", "20%", "80%", "90%", "95%")) +
  scale_x_continuous(breaks=seq(0, 0.06, 0.02), lim=c(0, 0.06)) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size=32),
    text = element_text(size=32),
    legend.justification = c(1.1, 1.1),
    legend.position = c(1, 1)
  )

dev.off() 




# B ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
kmer_counts <- read.table("data/figures/figure3/kmer_counts.txt", header=T, sep='\t')
order.kmer <- sort(kmer_counts$tot, index.return=T, decreasing=T)$ix
kmer_counts$kmer <- factor(kmer_counts$kmer, levels=kmer_counts$kmer[order.kmer])

kmer_counts2 <- kmer_counts[,c(1,2,3,5)]
kmer_counts2 <- melt(kmer_counts2)
kmer_counts2$category <- rep(NA, nrow(kmer_counts2))
for (i in 1:nrow(kmer_counts2)){kmer_counts2$category <- paste(kmer_counts$class, kmer_counts2$variable, sep='')}


svg(file = "figures/figure3/B_counts.svg", height = 5, width = 8)

ggplot(kmer_counts2, aes(x=kmer, y=value, fill=category)) +
  geom_bar(stat="identity", show.legend = FALSE) + 
  theme_classic() + 
  geom_hline(aes(yintercept = 20), size=0.5, linetype="twodash" ) +
  scale_fill_manual(values=c("grey60", "grey40", "gold", "blue", "gold", "blue", "gold", "blue")) + 
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size=32),
    text = element_text(size=32),
    legend.justification = c(1.1, 1.1),
    legend.position = c(1, 1)
  ) + 
  labs(x="Rank", y="Counts") + 
  scale_x_discrete(breaks=c(as.character(levels(kmer_counts$kmer)[1]),
                            as.character(levels(kmer_counts$kmer)[2500]),
                            as.character(levels(kmer_counts$kmer)[5000])), 
                            labels=c(1,2500, 5000) )
  
dev.off()



# C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
kmer_summary <- data.frame(counts=c(nrow(kmer_counts),
                                    sum( kmer_counts$class == "DT" ), 
                                    sum( kmer_counts$class == "nonDT-" ),
                                    sum( kmer_counts$class == "nonDT+")), 
                           class=c("All", "DT", "nDT_lo", "nDT_hi"))
kmer_summary$class <- factor(kmer_summary$class, levels=kmer_summary$class)


svg(file = "figures/figure3/C_summary.svg", height = 4, width = 4)

ggplot(kmer_summary, aes(x=class, y=counts, fill=class))+
  geom_col(show.legend = FALSE) + 
  theme_classic() +
  scale_y_continuous( breaks=c(0, 2500, 5000) ) + 
  scale_fill_manual(values = c("grey40", "orange", "#2ABE29", "#29ABE2") ) + 
  theme(
    axis.text.x = element_text(size=32, angle = 45, hjust = 1), 
    axis.text.y = element_text(size=32),
    axis.line.x = element_blank(), 
    axis.ticks.x = element_blank(),
    text = element_text(size=24)
  ) + 
  labs(x="", y="Counts")

dev.off()



# D ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
kmer_freqs <- read.table("data/figures/figure3/kmer_frequencies.txt", header=T, sep='\t')
kmer_freqs <- kmer_freqs[,c(1,2,5,6,3,4)]
colnames(kmer_freqs) <- c("codon", "DT", "nDT_hi", "nDT_lo", "redundant", "background")
kmer_m <- melt(kmer_freqs)


svg(file = "figures/figure3/D_codons.svg", height = 3, width = 16)

ggplot(data=kmer_m, aes(x=codon, y=value, fill=variable)) +
    geom_bar(stat='identity', position=position_dodge() ) + 
    theme_classic() +
    labs(x = "", y = "Frequency") +
    scale_fill_manual(values=c("orange", "#29ABE2", "#2ABE29", "grey40", "black")) + 
    theme(axis.text.x = element_text(size=18, family='Courier', angle = 90, hjust = 1), 
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank(), 
          axis.text.y = element_text(size=24),
          text = element_text(size=32),
          legend.text = element_text(size=24)
          )

dev.off()



# SUPPL ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

redundancy <- read.table("data/figures/figure3/kmer_redundancy.txt", header=T)
red_m <- melt(redundancy, id='length')


svg(file = "figures/figure3/SUPPL_redundancy.svg", height = 4, width = 5)

ggplot(red_m, aes(x=length, y=value+1, fill=variable)) + 
  geom_col( position="dodge2" ) + 
  scale_y_log10() + 
  scale_fill_manual(values=c("grey20", "grey50", "grey80")) + 
  labs(x="Length", y="Counts") + 
  theme_classic()  +
  theme(
    text = element_text(size=20),
    legend.title = element_blank()
  )


dev.off()


redundancy_mm <- read.table("data/figures/figure3/kmer_redundancy2.txt", header=T)
red2_m <- melt(redundancy_mm, id='length')


svg(file = "figures/figure3/SUPPL_redundancy_MM.svg", height = 4, width = 5)

ggplot(red2_m, aes(x=length, y=value+1, fill=variable)) + 
  geom_col( position="dodge2" ) + 
  scale_y_log10() + 
  scale_fill_manual(values=c("grey20", "grey50", "grey80")) + 
  labs(x="Length", y="Counts") + 
  theme_classic() + 
  theme(
    text = element_text(size=20),
    legend.title = element_blank()
  )

dev.off()
