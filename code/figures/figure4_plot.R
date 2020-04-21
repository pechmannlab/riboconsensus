# figure 4 panels
# oct 2019

library(ggplot2)
library(reshape2)
library(cowplot)


# Figure 4A ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clusters <- read.table("data/figures/figure4/clusters.txt", header=T, sep='\t')
colnames(clusters) <- c("category", "cluster", "not cluster")
clusters_m <- melt(clusters)


svg(file = "figures/figure4/A_clusters.svg", height = 5, width = 4)

ggplot(clusters_m, aes(x=category, y=value, fill=variable)) + 
  geom_col(position="dodge2") + 
  theme_classic() + 
  scale_fill_manual(values=c("grey40", "grey70")) + 
  labs(x="", y="Counts") +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size=24),
    text = element_text(size=24), 
    legend.position = "bottom", 
    legend.title = element_blank(),
    legend.text = element_text(size=24)
  )

dev.off()

#fisher.test(clusters[,c(2,3)])$p.value

# Figure 4B ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DTbias <- read.table("data/figures/figure4/cluster_kmer_position.txt", header=T, sep='\t')
BGbias <- read.table("data/figures/figure4/cluster_bg_position.txt", header=T, sep='\t')
positionbias <- rbind(BGbias, DTbias)


svg(file = "figures/figure4/B_positions.svg", height = 5, width = 3)

ggplot(positionbias, aes(x=class, y=position, fill=factor(class))) + 
  geom_boxplot(show.legend = F) + 
  lims(y=c(0, 1000)) + 
  labs(x="", y="Position") + 
  scale_fill_manual(values=c("grey40", "grey70", "grey40", "grey70")) + 
  theme_classic() +
  theme(
    axis.text.x = element_text(size=24, angle = 90, hjust = 1), 
    axis.text = element_text(size=24), 
    text = element_text(size=24), 
    axis.ticks = element_blank(),
    legend.title = element_blank(),
    legend.position = "bottom"
  )

dev.off()



# Figure 4C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DT       <- read.table("data/figures/figure4/associations_kmers_DT.txt", header=T, sep='\t')
nDThi    <- read.table("data/figures/figure4/associations_kmers_nDT+.txt", header=T, sep='\t')
nDTlo    <- read.table("data/figures/figure4/associations_kmers_nDT-.txt", header=T, sep='\t')


svg(file = "figures/figure4/C_DT_associations.svg", height = 5, width = 6)

dataDT <- data.frame(feature=factor(c("nOpt", "RNA", "+chrg"), levels=c("nOpt", "RNA", "+chrg")), DT=log2(DT$OR[1:3]), nDThi=log2(nDThi$OR[1:3]), nDTlo=log2(nDTlo$OR[1:3]) )
dataDT_m <- melt(dataDT)

data_sig <- as.vector(t(as.matrix( cbind( nDTlo$pval[1:3], nDThi$pval[1:3], DT$pval[1:3]) )) ) # weird re-ordering in ggplot, order of colour and alpha vectors don't correspond to data (!?!)

sel.sig <- data_sig < 0.05 
data_sig[sel.sig] <- 1
data_sig[!sel.sig] <- 0

data_alpha <- data_sig
sel.alpha <- data_alpha == 1
data_alpha[sel.alpha] <- 1
data_alpha[!sel.alpha] <- 0.4


ggplot(data=dataDT_m, aes(x=feature, y=value, fill=variable)) +
  geom_hline(aes(yintercept = 0), size=0.2) + 
  geom_bar(stat='identity', alpha=data_alpha, colour=data_sig,  position=position_dodge() ) + 
  theme_classic() + 
  coord_flip() + 
  labs(y = "log2 OR", x = "") +
  scale_fill_manual(values=c("orange", "#29ABE2", "#2ABE29")) + 
  theme(axis.text.x = element_text(size=24), 
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(), 
        axis.text.y = element_text(size=24),
        text = element_text(size=24),
        legend.text = element_text(size=24),
        legend.title = element_blank(),
        legend.position = c(0.8,0.8)
  )

dev.off()




# SUPPL A ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cluster_lenbias <- read.table("data/figures/figure4/suppl_cluster_kmer_length.txt", header=T)

plot.hi <- ggplot(cluster_lenbias) + 
  geom_point(aes(x=length, y=nclust_hi)) + 
  labs(x="Length", y="Number of high clusters") + 
  theme_classic() + 
  theme(
    text = element_text(size=16)
  )

plot.lo <- ggplot(cluster_lenbias) + 
  geom_point(aes(x=length, y=nclust_low)) + 
  labs(x="Length", y="Number of low clusters") + 
  theme_classic() + 
  theme(
    text = element_text(size=16)
  )

plot.cor <- ggplot(cluster_lenbias) + 
  geom_point(aes(x=nclust_hi, y=nclust_low)) + 
  labs(x="Number of high clusters", y="Number of low clusters") + 
  theme_classic() + 
  theme(
    text = element_text(size=16)
  )


svg(file = "figures/figure4/SUPPL_length_clusters.svg", height = 3, width = 8)

plot_grid(plot.hi, plot.lo, plot.cor, labels ="", ncol = 3, align = 'h')

dev.off()




svg(file = "figures/figure4/SUPPL_position_clusters.svg", height = 4, width =6)

cluster_posbias <- read.table("data/figures/figure4/suppl_DT_clusters.txt", header=T)

ggplot(cluster_posbias) + 
  geom_density(aes(x=rel_position, color=class), size=1.5) + 
  scale_color_manual(values=c("blue", "gold")) + 
  theme_classic() + 
  labs(x="Relative position", y="Density") + 
  theme(
    text = element_text(size=18), 
    legend.title = element_blank()
  )

dev.off()