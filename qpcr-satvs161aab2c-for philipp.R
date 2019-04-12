library(ggplot2)
library(dplyr)

pl1.n <- read.csv("gpr161aab-2c.csv", header=T)

pl1.plot <- pl1.n %>% filter(target %in% c("gpr161a", "gpr161b")) %>% filter(!sample == "gpr161aab.2")




# ##### plotting and statistics
# labels=c(expression(atop(italic('smo')^'+/?', 'control')),
#          expression(atop(italic('smo')^'-/-', 'control')),
#          expression(atop(italic('smo')^'-/-', +italic('NvSmo')))
# )
plotlabels <- c("sat", "aabb")

rm(expr_boxplot)
expr_boxplot <- ggplot(pl1.plot, aes(x = group,
                                    y = d0.norm,
                                    fill = group)) +
  
  geom_boxplot(aes(x = group,
                   y = d0.norm,
                   fill = group),
               position = position_dodge(width = 0.9),
               #size = 0.75,
               alpha = 1,
               outlier.alpha = 0,
               show.legend = F) +
  
  geom_point(aes(x = group,
                 y = d0.norm,
                 fill = group),
             position = position_dodge(width = 0.9),
             size = 1.5,
             alpha = 0.5,
             show.legend = F) +
  
  #stat_compare_means(comparisons = list(c("gpr161aabb", "wt"))) +
  
  facet_wrap( ~ target
              #, scales = "free"
              , strip.position = "bottom", nrow = 1) +
 
  
  
  theme_classic(base_size = 10) +
  
  theme(strip.background = element_rect(fill="#e0e0e0",size = 0),
        strip.text = element_text(face = "bold.italic", size = "10"),
        strip.placement = "outside"
  ) +
  
  theme(axis.title.x = element_blank(),
        axis.ticks = element_line(size=0.75),
        axis.text = element_text(size="8"),
        #axis.line = element_line(size=1)
        #panel.spacing = unit(1.5, "lines")
        ) +
  
  scale_y_continuous(trans=log2_trans(),
                     breaks = trans_breaks("log2", function(x) 2^x, n=4),
                     labels = trans_format("log2", math_format(2^.x))) +
  
  scale_fill_brewer(type = "seq", palette = "Greys") +
  
  #scale_x_discrete(labels=plotlabels) +
  labs(y="log2(normalized expression)") +
  expand_limits(x=2^-4)
  

expr_boxplot

summa <- summarize.ex(pl1.n)
summa <- mutate(summa, cv = sd.ex/mean.ex*100)
summa
