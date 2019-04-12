source("qpcr_functions.R")


pl1 <- process.run("dominik_2018-09-04 sat vs gpr161aab+ 2cell.rdml", opt=T)
pl1 <- pl1 %>% group_by(sample, target) %>% mutate(group = strsplit(sample, "\\.")[[1]][1])

pl1.n <- normalize(pl1, c("ef1a", "rpl13"))

write.csv(pl1.n, "2018-09-04 sat vs gpr161aab+ 2cell.csv")

pl1.plot <- pl1.n %>% filter(target %in% c("gpr161b", "gpr161a")) %>%
  filter(!sample == "gpr161aab.2") # removed outlier sample


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
  
  stat_compare_means(comparisons = list(c("gpr161aab", "wt"))) +
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
        axis.text = element_text(size="8")
        ) +
  scale_y_continuous(trans=log2_trans(),
                     breaks = trans_breaks("log2", function(x) 2^x, n=4),
                     labels = trans_format("log2", math_format(2^.x))) +
  
  scale_fill_brewer(type = "seq", palette = "Greys") +
  labs(y="log2(normalized expression)")

expr_boxplot

egb <- ggplot_build(expr_boxplot)
pvals <- p.adjust(c(as.character(egb$data[[3]]$annotation)), method = "BH")
psig <-  ifelse( pvals <= 0.001, "***",
                 ifelse( pvals <= 0.01, "**",
                         ifelse( pvals <= 0.05, "*",
                                 "")
                 )
)

egb$data[[3]]$annotation <- psig
#egb$data[[3]]$vjust <- egb$data[[3]]$vjust+0.5
plot(ggplot_gtable(egb))

ggsave("gpr161aab vs wt-2cell.pdf", ggplot_gtable(egb), device="pdf", width = 16, height = 9, units = "cm")