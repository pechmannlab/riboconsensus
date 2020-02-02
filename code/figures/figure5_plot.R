# figure 5 panels

library(ggplot2)
library(reshape2)



# A --------------------------------------------------------------------------------------------------------
ssb <- as.data.frame(read.table("data/figures/figure5/ssb_enrichment.txt", header=T, sep='\t'))
ssb$position <- ssb$position - 30
ssb.DTkmers <- as.data.frame( read.table('data/figures/figure5/ssb.DTkmers', header=T, sep='\t') )
ssb.DTkmers$position <- ssb.DTkmers$position - 30
ssbpoly <- data.frame(id = rep(c("site"), each=2), x = c(0, 7, 7, 0), y = c(rep(c(0, 0.002), each=2) ) )


svg(file = "figures/figure5/A_ssbbinding.svg", height = 5, width = 7)

ggplot(ssb, aes(x=position)) + 
  geom_ribbon(aes(ymax=rand_mean+rand_sd, ymin=rand_mean-rand_sd), fill="#eeeeee", alpha=1) +
  geom_line(aes(x=position, y=firstssb), color="red", size=1.5) + 
  geom_line(aes(x=position, y=allssb), color="orange", size=1) + 
  theme_classic() + 
  lims(x=c(-30,40), y=c(0, 0.05)) + 
  labs(x="Position relative to binding", y="Avg. consensus RD") +
  geom_freqpoly(inherit.aes=F, data=ssb.DTkmers, aes(x=position, y=..density.., color=factor(class) ), 
                linetype=c("twodash"), binwidth=3, show.legend = FALSE) +
  scale_color_manual(values=c("grey20", "grey60")) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size=32),
    text = element_text(size=32)
  ) + 
  geom_segment(aes(x=-30, xend=40, y=0.001, yend=0.001), arrow=NULL)  + 
  geom_polygon(data=ssbpoly, inherit.aes=F, aes(x=x, y=y, fill=id), show.legend= FALSE) + 
  scale_fill_manual(values=(c("gold", "#555555")))

dev.off()




# B --------------------------------------------------------------------------------------------------------
domain <- as.data.frame(read.table("data/figures/figure5/domain_RD.txt", header=T, sep='\t'))
domain.DTkmers <- as.data.frame( read.table('data/figures/figure5/domain.DTkmers', header=T, sep='\t') )
domainpoly <- data.frame(id = rep(c("db", "tunnel"), each=4), x = c(0, 2, 2, 0, 2, 37, 37, 2), y = c(rep(c(0, 0.002), each=2), rep(c(0,0.002),each=2)) )


svg(file = "figures/figure5/B_domainboundaries.svg", height = 5, width = 7)

ggplot(domain, aes(x=position)) + 
  geom_ribbon(aes(ymax=rand_mean+rand_sd, ymin=rand_mean-rand_sd), fill="#eeeeee", alpha=1) +
  geom_line(aes(x=position, y=first_db), color="red", size=1.5) + 
  geom_line(aes(x=position, y=nonssb), color="#29ABE2", size=1) + 
  theme_classic() + 
  lims(x=c(0, 70), y=c(0, 0.035)) + 
  labs(x="Distance from domain boundary", y="Avg. consensus RD") +
  geom_freqpoly(inherit.aes=F, data=domain.DTkmers, aes(x=position, y=..density.., color=factor(class) ), 
                linetype=c("twodash"), binwidth=1, show.legend = FALSE) +
  scale_color_manual(values=c("grey20", "grey60")) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size=32),
    text = element_text(size=32)
  ) + 
  geom_segment(aes(x=0, xend=70, y=0.001, yend=0.001), arrow=NULL)  + 
  geom_polygon(data=domainpoly, inherit.aes=F, aes(x=x, y=y, fill=id), show.legend= FALSE) + 
  scale_fill_manual(values=(c("#2ABE29", "#555555")))

dev.off()

 
