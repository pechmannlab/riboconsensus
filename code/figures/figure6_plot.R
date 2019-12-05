# figure 6 plots

library(ggplot2)
library(reshape2)


GO <- read.table("data/figures/figure6/GO_combined.txt", header=T, sep='\t', na.strings="NA")
GO$GO <- factor(GO$GO, levels=GO$GO)
GOm <- melt(GO)
GOm <- GOm[,c(2,3,4,5)]
colnames(GOm) <- c("GO", "description", 'variable', 'pval')



svg(file = "figures/figure6/GO_analysis.svg", height = 9.5, width = 7)

ggplot(GOm, aes(x = variable, y = GO)) + 
  geom_tile(aes(fill=log10(pval)), color="white", size=0.8) + 
  scale_y_discrete(labels=GO$name, position="left") +
  theme_classic() + 
  labs(x=NULL, y=NULL) + 
  theme(axis.text.x = element_text(size=16, angle = 45, hjust = 1), 
        axis.text.y = element_text(size=12),
        axis.ticks = element_blank(), 
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        text = element_text(size=16),
        legend.justification = c(1.1, 1.1),
        legend.position = c(1, 1)
  ) +
  scale_fill_gradient(low = "forestgreen", high = "white", space = "Lab", 
                      na.value = "white", guide = "colourbar")

dev.off()




