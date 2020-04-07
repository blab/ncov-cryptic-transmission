# command line script to export pretty ggptree plots

library(ape,quietly = TRUE)        # Analyses of Phylogenetics and Evolution
library(ggtree,quietly = TRUE)     # pretty tree plots
library(cowplot,quietly = TRUE)    # nicer plot defaults
library(lubridate,quietly = TRUE)  # tools to help with date formatting
library(optparse,quietly = TRUE)   # command line argument parsing

option_list <- list(
  make_option(c("-tf", "--treeFile"), type="character", default='example/connectedTree.tre',help="tree file name", metavar="character"),
  make_option(c("-df", "--dataFile"), type="character", default="example/connectedTree.csv",help="data file name", metavar="character"),
  make_option(c("-cf", "--colorField"), type="character", default="timeInfected",help="color field", metavar="character"),
  make_option(c("-lp", "--legendPosition"), type="character", default="right",help="legend position (\'none\' == no legend)", metavar="character"),
  make_option(c("-of", "--outputFile"), type="character", default='derive',help="output file name", metavar="character"),
  make_option(c("-sc", "--scale"), type="character", default="time", help="x-axis scale (time, divergence)", metavar="character")
); 

opt <- parse_args(OptionParser(option_list=option_list))


phyTree <- ape::read.tree(opt$treeFile)
dat<-read.csv(opt$dataFile,stringsAsFactors = FALSE)
dat<-dat[!is.na(dat$timeInfected),]  # remove ancestral nodes with no metadata (from artificial joins)

dat$timeInfected <- date_decimal(dat$timeInfected)

if("PUMA5CE" %in% colnames(dat)){
  dat$PUMA5CE<-as.character(dat$PUMA5CE)
}
if("GEOID" %in% colnames(dat)){
  dat$GEOID<-as.character(dat$GEOID)
}

######################
## Plot with ggtree ##
######################

if( opt$scale=='time'){
  p <- ggtree(phyTree,mrsd = max(dat$timeInfected)) %<+% dat +                # creates visualization and binds metadata
      geom_tippoint(aes(color=get(opt$colorField))) +     # displays information about the sequences
      geom_nodepoint(aes(color=get(opt$colorField))) +    # displays information about internal nodes
      theme_tree2(legend.position=opt$legendPosition) + labs(color=opt$colorField) 
} else {
  p <- ggtree(phyTree) %<+% dat +                # creates visualization and binds metadata
    geom_tippoint(aes(color=get(opt$colorField))) +     # displays information about the sequences
    geom_nodepoint(aes(color=get(opt$colorField))) +    # displays information about internal nodes
    theme_tree2(legend.position=opt$legendPosition) + labs(color=opt$colorField) 
}

if (opt$outputFile=='derive'){
    opt$outputFile <- paste(sub('.tre','',opt$treeFile),opt$colorField,'png',sep='.')
}

ggsave(opt$outputFile, width=6, height = 5, units='in', dpi=300)

file.remove('Rplots.pdf')  # ggsave creates an empty device window from command line


