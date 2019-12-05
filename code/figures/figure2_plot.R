# Figure2
# Extracted from figures/figure_2/R_plotting_figure_2.ipynb notebook (02/8/19)

# Libraries:
library(jsonlite)
library(reshape2)
library(RColorBrewer)
library(tidyverse)
library(pROC)


# A ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TR_consensus <- as.data.frame( read.table("data/figures/figure2/codon_rates_consensus.txt", header=T, sep='\t') )
TR_consensus$nopt <- as.factor( as.numeric( TR_consensus$tAI <= quantile(TR_consensus$tAI, 0.25 ) ) )
order.TR <- sort(TR_consensus$median_mc, index.return=T)$ix
TR_consensus$Codon <- factor(TR_consensus$Codon, levels=TR_consensus$Codon[order.TR])
TR_consensus$rank <- 1:nrow(TR_consensus)


svg(file = "figures/figure2/A_codonTR.svg", height = 5, width = 8)

ggplot(TR_consensus, aes(x=Codon, y=median_mc)) + 
  #geom_hline(aes(yintercept = mean(median_mc)), size=0.5, linetype="twodash" ) + 
  geom_errorbar(aes(ymin=p_25_mc, ymax=p_75_mc), size = 1, width=0, color="grey30") +
  geom_point(aes(color=nopt), shape=15, size=3) + #, show.legend = FALSE) + 
  scale_color_manual(values = c('grey30', "#E31A1C")) + 
  
  theme_classic() +
  theme(axis.text.x = element_blank(),
        #axis.text.x=element_text(size=18, angle = 90, hjust = 1), 
        axis.ticks.x=element_blank(), 
        axis.text.y = element_text(size=32), 
        axis.line.x = element_blank(), 
        text = element_text(size=32), 
        legend.justification = c(1.1, 1.1),
        legend.position = c(1, 1)) + 
  labs(y = "Codon consensus RD", x="Codons")
  
dev.off()




# B ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CTR <- read.table("data/figures/figure2/codonrates_all.txt", header=T, sep='\t')
CTR$nopt <- as.numeric( CTR$tAI <= quantile(CTR$tAI, 0.25 ) )

roc_mc <- roc(CTR$nopt, CTR$median_mc, direction="<")
roc_mn <- roc(CTR$nopt, CTR$median_mn, direction="<")


reverse_rocdf <- function(roc){
  data <- data.frame(spec=roc$specificities, sens=roc$sensitivities)
  newspec <- c()
  newsens <- c()
  unique_spec <- unique(data$spec)
  for (i in 1:length(unique_spec)){
    sel <- data$spec == unique_spec[i]
    newspec <- c(newspec, rep(1-unique_spec[i], sum(sel)))
    newsens <- c(newsens, rev(data$sens[sel]) )
  }
  newdata <- data.frame(specificities=newspec, sensitivities=newsens)
}


roc_df <- data.frame()
result <- rep(0, 20)
for (i in 1:20){
  myroc <- roc(CTR$nopt, CTR[,i+3], direction="<")
  myroc_df <- reverse_rocdf(myroc)
  myroc_df$class <- rep( paste("R",i,sep=""), length(myroc_df$specificities) )
  myroc_df$weight <- rep(0.5, length(myroc_df$specificities)) 
  roc_df <- rbind(roc_df, myroc_df )
  result[i] <- myroc$auc
}

roc_mc_rev <- reverse_rocdf(roc_mc)
roc_mc_rev$class <- rep("MC",length(roc_mc_rev$specificities))
roc_mc_rev$weight <- rep(0.8, length(roc_mc_rev$specificities))

roc_mn_rev <- reverse_rocdf(roc_mn)
roc_mn_rev$class <- rep("MN",length(roc_mn_rev$specificities))
roc_mn_rev$weight <- rep(0.6, length(roc_mn_rev$specificities))

roc_df <- rbind(roc_df, roc_mc_rev, roc_mn_rev)


svg(file = "figures/figure2/B_roc.svg", height = 5, width = 6)

ggplot(roc_df, aes(x=specificities, y=sensitivities, color=class)) + 
  geom_abline(aes(intercept = 0, slope=1), size=0.5, linetype="twodash" ) +
  geom_line(aes(size=weight), show.legend = FALSE) + 
  theme_classic() + 
  scale_color_manual(values = c( "red", "orange",  rep("grey60", 20) ) ) + 
  scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1), labels=c(1, 0.75, 0.5, 0.25, 0)) + 
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1), labels=c(0, 0.25, 0.5, 0.75, 1)) + 
  theme(
    axis.text = element_text(size=32),
    text = element_text(size=32)
  ) + 
  labs(x="Specificity", y="Sensitivity")
  
dev.off()



svg(file = "figures/figure2/B2_roc.svg", height = 5, width = 4.5)

auc <- data.frame(class=factor(c("Consensus", "Mean", "Avg ind."), levels=c("Consensus", "Mean", "Avg ind.")), auc=c(roc_mc$auc, roc_mn$auc, mean(result)) )

ggplot(auc, aes(x=class, y=auc, fill=class)) + 
  geom_col() + 
  theme_classic() +
  scale_y_continuous( lim=c(0, 0.8) ) + 
  scale_fill_manual(values = c("red", "orange", "grey60") ) + 
  theme(
    axis.text.x = element_blank(), #element_text(size=32, angle = 45, hjust = 1), 
    axis.text.y = element_text(size=32),
    axis.line.x = element_blank(), 
    axis.ticks.x = element_blank(),
    text = element_text(size=24)
  ) + 
  labs(x="", y="AUC")

dev.off()







# C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
inhibitory <- c('CGACCG', 'CGAGCG', 'CGAATA', 'CTCCCG', 'CGACGA', 'CGACGG', 'CGACTG', 'GTACCG', 'GTGCGA', 'GTACGA', 'CTGCCG', 'CTGCGA', 'ATACGA', 'ATACGG', 'AGGCGG', 'CTGATA', 'AGGCGA')

cdr_cons <- read.table("data/figures/figure2/codonpair_rates_consensus.txt", header=T, sep='\t')
cdr_cons <- cdr_cons[rowSums(is.na(cdr_cons))==0,]

rankorder <- c()
inhib_idx <- c()
good_inhib <- c() # these are the ones that are present in the set of genes analyzed
for (i in 1:length(inhibitory)){o <- which(as.character(cdr_cons$Codonpair)==inhibitory[i]); if (length(o) > 0) {inhib_idx <- c(inhib_idx, o); good_inhib <- c(good_inhib, i)} }
res <- c()
for (i in 1:length(inhib_idx)){o <- sum(cdr_cons$median > cdr_cons$median[inhib_idx[i]]) + 1; res <- c(res, o)}
rankorder <- data.frame(mc=res)

inhib_dat <- data.frame( x=cdr_cons$median[inhib_idx], y=rep(120, length(inhib_idx)))

svg(file = "figures/figure2/C_codonpairs_histo.svg", height = 5, width = 6)

ggplot(cdr_cons, aes(x=median)) + 
  geom_histogram(aes(x = median),
                 binwidth = 0.005,
                 fill="#E31A1C",
                 color='black',
                 alpha=0.5
  ) +
  labs(x = "Codon duplet consensus RD", y = "Frequency") +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size=32),
    text = element_text(size=32)
  ) + 
  geom_segment(data=inhib_dat, inherit.aes=F, aes(x=x, xend=x, y=y, yend=y-100), arrow=arrow(length = unit(0.02, "npc"), ends="last", type="closed") )
  
dev.off()






# D ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# rankorder naive mean
rankorder_mn <- c()
inhib_idx_mn <- c()
for (i in 1:length(good_inhibitory)){o <- which(as.character(cdr_cons$Codonpair)==good_inhibitory[i]); if (length(o) > 0) {inhib_idx_mn <- c(inhib_idx_mn, o)} }
res_mn <- c()
for (i in 1:length(inhib_idx_mn)){o_mn <- sum(cdr_cons$mn_median > cdr_cons$mn_median[inhib_idx_mn[i]]) + 1; res_mn <- c(res_mn, o_mn)}
rankorder_mn <- data.frame(mn=res_mn)

# rankorder indv. 
cdr <- read.table("data/figures/figure2/codon_duplets_rates.txt", header=T, sep='\t')

cdc <- as.matrix(read.table("codon_duplet_counts.txt"))
cdc.sel <- colMeans(cdc) > 10
cdr <- cdr[,cdc.sel]

good_inhibitory <- inhibitory[good_inhib]

rankorder_ind <- c()
tmp_idx <- c()
for (j in 1:length(good_inhibitory) ){
  o <- which(colnames(cdr) == good_inhibitory[j])
  if ( length(o) > 0) {tmp_idx <- c(tmp_idx, o)} 
}
for (i in 1:20){
  res <- c()
  for (j in 1:length(tmp_idx) ){
    o <- sum( cdr[i,] > cdr[i, tmp_idx[j] ] , na.rm=T) + 1
    res <- c(res, o)
  }
  rankorder_ind <- cbind(rankorder_ind, res)
}
colnames(rankorder_ind) <- paste("R",1:20,sep='')
rankorder_ind <- as.data.frame(rankorder_ind)

ranko_df <- data.frame( rbind(  data.frame( class=rep("mc", nrow(rankorder)),  rank=as.numeric(as.matrix(rankorder))), 
                                data.frame( class=rep("mn", nrow(rankorder_mn)), rank=as.numeric(as.matrix(rankorder_mn))), 
                                data.frame( class=rep("ind", prod(dim(rankorder_ind))), rank=as.vector(as.matrix(rankorder_ind)))
                      ))
  
  

svg(file = "figures/figure2/D_codonpairs_rank.svg", height = 5, width = 4)

ggplot(ranko_df, aes(x=class, y=rank, fill=class)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width=0.2) + 
  theme_classic() +
  scale_fill_manual(values = c("red", "orange", "grey60") ) + 
  theme(
    axis.text.x = element_blank(), #element_text(size=32, angle = 45, hjust = 1), 
    axis.text.y = element_text(size=32),
    axis.line.x = element_blank(), 
    axis.ticks.x = element_blank(),
    text = element_text(size=32), 
    legend.justification = c(1.1, 1.1),
    legend.position = c(0.6, 1)
  ) + 
  labs(x="", y="Rank")

dev.off()





# E ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
RD <- read.table("data/figures/figure2/codonpair_rates_raw10_cons.txt", header=T, sep='\t')
sel.mc <- RD$class == 'mc'
sel.mn <- RD$class == 'mn'

RD.norm <- RD
# add count of 1 so that mean is 1 and negative values (which don't make sense) are limited
RD.norm[sel.mc, 3] <- ( RD.norm[sel.mc, 3] - mean(RD.norm[sel.mc, 3])) / sd(RD.norm[sel.mc, 3])  + 1
RD.norm[sel.mn, 3] <- ( RD.norm[sel.mn, 3] - mean(RD.norm[sel.mn, 3])) / sd(RD.norm[sel.mn, 3]) + 1


svg(file = "figures/figure2/E_RD_normdist.svg", height = 5, width = 5)

ggplot(RD.norm, aes(x=RD, color=class)) + 
  geom_density(size=1.2) + 
  theme_classic() + 
  scale_x_continuous(lim=c(0, 4)) + 
  scale_color_manual(values=c("red", "orange")) + 
  labs(x="Normalized RD", y="Density") + 
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size=32),
    text = element_text(size=32),
    legend.justification = c(1.1, 1.1),
    legend.position = c(0.75, 1)
  )
 
dev.off()