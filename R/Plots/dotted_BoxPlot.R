r2 <- read.csv('../TestforR_Average.txt', sep = '\t')

posn.d <- position_dodge(width=0.4)
my_color = rep('dodgerblue1', times = 6)
ggplot(r2, aes(x=Treatment, y=Average.per.cell, color = Embryo)) + 
  geom_point(alpha = 1, position = posn.d) + scale_color_manual(values = my_color) +
  geom_boxplot(alpha = 0, colour = "black") + 
  theme_classic() + ylim(c(0,2)) +
  theme(legend.position = 'none')
