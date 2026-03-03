# =====================================================================
#  Hierarchical Division in R
#  Goal: Hierarchical structure division of the co-expression network.
#
#  Input:
#    - summary_overall_cultivar.txt : Features of each module in the co-expression network.
#
#  Output:
#    - Sunburst plot with customized theme and manual colors
#
#  Dependencies:
#    - MEGENA ggrepel dplyr
# =====================================================================

module.table=read.table("summary_overall_cultivar.txt",header=T)
htbl = module.table[,c("module.parent","module.id")]
htbl$category="A"

##Supplement information for the top-level parent modules
root1=data.frame(module.parent="c1_1",module.id="c1_1",category="B")
root2=data.frame(module.parent="c2_1",module.id="c2_1",category="B")
df <- rbind(htbl,root1,root2)
df$category=as.factor(df$category)

sbobj3 = draw_sunburst_wt_fill(module.df = df,
                               feat.col = "category",
                               fill.type = "discrete",border.col="white",border.width=0.5,
                               fill.scale = scale_fill_manual(values = c("A" = "gray","B" = "gray")), 
                               id.col = "module.id",parent.col = "module.parent",min.angle=400,
                               theme.adjust=theme_void())