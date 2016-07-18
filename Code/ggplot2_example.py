library(ggplot2)

library(grid)

GenPhen_plot <-function(name_of_file){
 
 
 sim_scatter <- read.delim(name_of_file, header = FALSE)
 # ex. "/Users/bianca/Documents/graduatevwork/thesis/paper_plots_scripts/data_files/bveQTL_0.3_1000.0_0.0_0.0_1.0_-0.1_sim_no_45.txt"

 colnames(sim_scatter) <- c("genotype","expression")
 sim_scatter$genotype[sim_scatter$genotype == 0] <- "AA"
 sim_scatter$genotype[sim_scatter$genotype == 1] <- "Aa"
 sim_scatter$genotype[sim_scatter$genotype == 2] <- "aa"

 g<- ggplot(sim_scatter, aes(genotype, expression, color=factor(genotype))) + scale_color_manual(values=c("lightseagreen", "darkorange1","lightslateblue"))
 #g <- g + geom_violin(alpha=0.2, color="black")
 g <- g + xlab("genotype") + ylab("phenotype") 
 g<- g + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
               panel.background = element_blank(), axis.line = element_line(colour = "black"))
 g <- g + theme(text = element_text(size=18), axis.text.x = element_text(size=18),  axis.text.y = element_text(size=18),
                legend.title=element_blank(),legend.position=c(0.85,0.95),legend.background=element_blank(),
                legend.key=element_rect(fill=NA))                                                                                                                    
 #g <- g + guides(colour = guide_legend(override.aes = list(size=2))) 
 g <- g + theme(legend.position="none")
 #g <- g+geom_violin(alpha=0.5, color="gray") + geom_jitter(alpha=0.5,
 #                                                          position = position_jitter(width = 0.1))
 g <- g+geom_boxplot(alpha=0.5, color="gray")+geom_jitter(alpha=0.5,
                                                          position = position_jitter(width = 0.1)) 
 g <- g + theme(aspect.ratio = 1 . )
 
 print(g)
 return(g)
}