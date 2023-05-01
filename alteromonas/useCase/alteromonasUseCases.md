# *Alteromonas* Use Cases

**Date:** April 2023
**Primary Author:** Michelle DeMers

This reproducible workflow details the work done in analyzing the CAZyme availability in *Alteromonas*. The first use case looks at the distribution of CAZymes across the phylogeny of *Alteromonas*. The second use case looks at the phylogeny created from the CAZymes, determining how they align with the accepted phylogeny of the genus.

## Use Case 1: CAZyme availability in the genus

### Create Split Pangenome

To analyze the CAZymes directly, it is necessary to create a split pangenome. Within the interactive view of the pangenome, the `Search` tab was navigated to, the `Field` was set to `CAZyme`, `Operator` was set to `==`, and `Term` was set to `KNOWN`. These splits were highlighted on the tree. Then, a bin was made from this selection under the default name, `Bin_1`, in the `default_2` collection. This bin was then turned into its own pangenome with the command:

	anvi-split -p Alteromonas_Pangenome_2/Alteromonas_Pangenome_2-PAN.db -g Alteromonas-v2-GENOMES.db -C default_2_ --bin-id Bin_1 -o CAZYMES
	
This resulting pangenome has 680 gene clusters, reflective of the 680 gene clusters with known CAZyme functions.

### Count CAZymes

A summary was made from the split pangenome by making a defult collection:

	# make default collection in new pangenome
	anvi-script-add-default-collection -p CAZYMES/Bin\_1/PAN.db -g Alteromonas-NEW-GENOMES.db -C DEFAULT -b EVERYTHING
	#make summary
	anvi-summarize -p CAZYMES/Bin\_1/PAN.db -g Alteromonas-NEW-GENOMES.db -C DEFAULT
	
From there, the resulting summary file, `The_Bin_1_split_from_"Alteromonas_Pangenome_2"_gene_clusters_summary.txt`, was turned into an Excel file. Then, a pivot table was made from the genomes v. CAZyme data, with each CAZyme converted into its category. The Pivot Table was built with the Category as columns, the genome names as rows, and the values as the Count of the individual CAZyme types. This pivot table was saved as a .csv file, `Cazymes_per_genome_mags.csv`.

### Processing

To process these results, the CAZyme count file `Cazymes_per_genome_mags.csv` was pulled into R. Processing also needed the pangenome storage file `external-genomes.txt`, the phylogeny file `midpoint_mags_bayesian`, the genome names file `view.txt`, and the genome metadata file `alteromonas_metadata.csv`.

```r
## Get libraries
library("plyr")
library("BiocManager")
library("ggtree")
library("ggplot2")
library("dplyr")
library(stringr)
library("aplot")
library("ggtreeExtra")
library("ggstance")
library("tidyr")
library("ggnewscale")
library("nationalparkcolors")

## Formatting data for the tree
genomesnames <- read.table("path/to/pangenome/files/external-genomes.txt", header=T)
#here we have to go through all of the additional MAGs and format them correctly
test <- grep("contigs", genomesnames$contigs_db_path)
test2 <- genomesnames[test,]
for (p in 1:nrow(test2)){
  num <- str_split(test2$contigs_db_path[p], "_")[[1]][4]
  num <- str_split(num,"-")[[1]][1]
  genomesnames[test[p],2] <- num
}
genomesnames$contigs_db_path<- gsub(".db", "", genomesnames$contigs_db_path)
genomesnames$contigs_db_path<- gsub("GCF_", "", genomesnames$contigs_db_path)

## add the completeness map
completeness <- read.csv("path/to/metdata/alteromonas_metadata.csv")[,c(1,5,38)]
test <- completeness$taxon_oid %in% genomesnames$contigs_db_path
completeness <- completeness[test,]
completeness<-completeness[complete.cases(completeness),]
row.names(completeness)<-completeness$taxon_oid
new1<-match(completeness$taxon_oid,genomesnames$contigs_db_path)
names1 <- genomesnames[new1,1]
row.names(completeness) <- names1
completeness <- data.frame("Completeness"=completeness[,3])
row.names(completeness) <- names1

## add the mags/isolates map
mags <- read.csv("path/to/metdata/alteromonas_metadata.csv")[,c(1,5,42)]
test <- mags$taxon_oid %in% genomesnames$contigs_db_path
mags <- mags[test,]
row.names(mags)<-mags$taxon_oid
new1<-match(mags$taxon_oid,genomesnames$contigs_db_path)
names1 <- genomesnames[new1,1]
row.names(mags) <- names1
mags$Genome.Name...Sample.Name <- names1
mags <- data.frame("Type of Genome"=mags[,3])
row.names(mags) <- names1
colnames(mags)<-c("Type of Genome")

## Import the tree
tree <- read.tree("path/to/tree/midpoint_mags_bayesian")
# get the labels on the tree properly
clade_position <- data.frame(tree$tip.label,
                             c(1:336)
)
colnames(clade_position)<-c("label", "order")

## Import tree tip labels
genomes <- read.table("path/to/metadata/view.txt", sep ="\t", header=T
)
new_order<-match(clade_position$label,genomes$genome_id)
label2_names <- genomes[new_order,2]
d <- data.frame(label = tree$tip.label, label2=label2_names)
tree2 <- full_join(tree, d, by = "label")

## To format the tree with labels
p<-ggtree(tree2) + theme_tree2() + geom_tippoint(size=2) + ggtitle("All Genomes Tree") +theme(legend.position = "none", axis.title.y = element_blank(), plot.title = element_text(size = 12,face = "bold",hjust = 0.5,vjust = 1)) + hexpand(.15) +
  geom_cladelabel(node=523, label="A. macleodii", color=(brewer.pal(10, "Paired"))[2], offset=.0, angle=90, hjust='center', offset.text=.03, barsize=1.5, fontsize=8) + 
  geom_cladelabel(node=494, label="A. mediterranea", color=(brewer.pal(10, "Paired"))[4], offset=.0, angle=90, hjust='center', offset.text=.03, barsize=1.5, fontsize=8) + 
  geom_cladelabel(node=390, label="A. australica", color=(brewer.pal(10, "Paired"))[6], offset=.0, angle=90, hjust='center', offset.text=.04, barsize=1.5, fontsize=8)  + 
  geom_cladelabel(node=463, label="A. naphthalenivorans", color=(brewer.pal(10, "Paired"))[8], offset=.0, angle=90, hjust='center', offset.text=.02, barsize=1.5, fontsize=3)  + 
  geom_cladelabel(node=457, label="A. stellipolaris", color=(brewer.pal(10, "Paired"))[10], offset=.0, hjust='center', angle=90, offset.text=.02, barsize=1.5, fontsize=3) + 
  geom_hilight(node=523, fill=(brewer.pal(10, "Paired"))[1], alpha=.6, extendto = 1.36) + 
  geom_hilight(node=494, fill=(brewer.pal(10, "Paired"))[3], alpha=.6, extendto = 1.4) + 
  geom_hilight(node=390, fill=(brewer.pal(10, "Paired"))[5], alpha=.6, extendto = 1.31)+ 
  geom_hilight(node=463, fill=(brewer.pal(10, "Paired"))[7], alpha=.6, extendto = 1.34) + 
  geom_hilight(node=457, fill=(brewer.pal(10, "Paired"))[9], alpha=.6, extendto = 1.315)

## Adding Metadata
h2 <- gheatmap(p, completeness, offset = .02,                               # offset shifts the heatmap to the right,
               width = 0.07,                              # width defines the width of the heatmap column,
               # color defines the boarder of the heatmap columns
               colnames = T, colnames_angle=0, font.size=2, color=F, hjust=.5) + theme(legend.position = "bottom", legend.title = element_text(size = 10), legend.text = element_text(size = 5), legend.box = "horizontal", legend.margin = margin()) + scale_fill_viridis_c( direction=-1, name = "Completeness")
h2 <- h2 + new_scale_fill()
h2 <- gheatmap(h2, mags, offset = .17,width = 0.07, colnames = T, colnames_angle=0, font.size=2, color=F, hjust=.5) + scale_fill_manual(values = c(brewer.pal(7, name = "Greys" )[2], brewer.pal(7, name = "Greys" )[7]), name = "Type of Genome")
h2 <- h2 + new_scale_fill()

## To format the CAZyme data:
d2 <- read.csv("path/to/metadata/Cazymes_per_genome_mags.csv")
d2 <- d2[,1:7]
tree$tip.label
order<-match( tree$tip.label, d2$Genome)
d2 <- d2[order,]
d2$Genome == tree$tip.label # check
d2 <- gather(d2, Cazyme, Count, AA:PL)
aa <- d2[d2$Cazyme == "AA",]
ce <- d2[d2$Cazyme == "CE",]
pl <- d2[d2$Cazyme == "PL",]
cbm <-  d2[d2$Cazyme == "CBM",]
gt <-  d2[d2$Cazyme == "GT",]
gh <- d2[d2$Cazyme == "GH",]

## To make proportion values, where each denominator is the max from that column:
aa$Count <- aa$Count/19
ce$Count <- ce$Count/26
pl$Count <- pl$Count/22
cbm$Count<-  cbm$Count/5
gt$Count <-  gt$Count/48
gh$Count <- gh$Count/76
barpal <- park_palette("ArcticGates")

## To make the barplot on the side of the figure:
p4 <- facet_plot(h2, panel = 'Carbohydrate-Binding Modules', data = cbm, 
                 geom = geom_barh, 
                 mapping = aes(x = Count, fill = as.factor(Cazyme)), 
                 stat='identity' , fill = barpal[1], color = "#3F3939") + theme(legend.position = "bottom", legend.text = element_text(size = 7)) + scale_x_continuous(limits=c(0,1))
p5 <- facet_plot(p4, panel = 'Auxiliary Activites', data = aa, 
                 geom = geom_barh, 
                 mapping = aes(x = Count, fill = as.factor(Cazyme)), 
                 stat='identity' , fill = barpal[2], color = "#3F3939") + theme(legend.position = "bottom", legend.text = element_text(size = 7))
p6 <- facet_plot(p5, panel = 'Polysaccharide Lyases', data = pl, 
                 geom = geom_barh, 
                 mapping = aes(x = Count, fill = as.factor(Cazyme)), 
                 stat='identity' , fill = barpal[3], color = "#3F3939") + theme(legend.position = "bottom", legend.text = element_text(size = 7))
p7 <- facet_plot(p6, panel = 'Carbohydrate Esterase', data = ce, 
                 geom = geom_barh, 
                 mapping = aes(x = Count, fill = as.factor(Cazyme)), 
                 stat='identity' , fill = barpal[4], color = "#3F3939") + theme(legend.position = "bottom", legend.text = element_text(size = 7))
p8 <- facet_plot(p7, panel = 'GlycosylTransferases', data = gt, 
                 geom = geom_barh, 
                 mapping = aes(x = Count, fill = as.factor(Cazyme)), 
                 stat='identity' , fill = barpal[5], color = "#3F3939") + theme(legend.position = "bottom", legend.text = element_text(size = 7))

p9 <- facet_plot(p8, panel = 'Glycoside Hydrolases', data = gh, 
                 geom = geom_barh, 
                 mapping = aes(x = Count, fill = as.factor(Cazyme)), 
                 stat='identity' , fill = barpal[6], color = "#3F3939") + theme(legend.position = "bottom", legend.text = element_text(size = 8))
p9 <- facet_widths(p9, widths = c(7, 1,1,1,1,1,1))
ggsave(paste0("path/to/final/files/Cazymes_All_Heatmap_and_tree1.pdf"), p9, width=20, height = 35, units = "in")

## The figure above will produce a cut-off version of the tree. The code was re-run without the `scale_x_continuous(limits=c(0,1))` command, and saved as a new file. Then the left (tree) and right (bar plots) were combined into one figure:

p<-ggtree(tree2) + theme_tree2() + geom_tippoint(size=2) + ggtitle("All Genomes Tree") +theme(legend.position = "none", axis.title.y = element_blank(), plot.title = element_text(size = 12,face = "bold",hjust = 0.5,vjust = 1)) + hexpand(.15) +
  geom_cladelabel(node=523, label="A. macleodii", color=(brewer.pal(10, "Paired"))[2], offset=.0, angle=90, hjust='center', offset.text=.03, barsize=1.5, fontsize=8) + 
  geom_cladelabel(node=494, label="A. mediterranea", color=(brewer.pal(10, "Paired"))[4], offset=.0, angle=90, hjust='center', offset.text=.03, barsize=1.5, fontsize=8) + 
  geom_cladelabel(node=390, label="A. australica", color=(brewer.pal(10, "Paired"))[6], offset=.0, angle=90, hjust='center', offset.text=.04, barsize=1.5, fontsize=8)  + 
  geom_cladelabel(node=463, label="A. naphthalenivorans", color=(brewer.pal(10, "Paired"))[8], offset=.0, angle=90, hjust='center', offset.text=.02, barsize=1.5, fontsize=3)  + 
  geom_cladelabel(node=457, label="A. stellipolaris", color=(brewer.pal(10, "Paired"))[10], offset=.0, hjust='center', angle=90, offset.text=.02, barsize=1.5, fontsize=3) + 
  geom_hilight(node=523, fill=(brewer.pal(10, "Paired"))[1], alpha=.6, extendto = 1.36) + 
  geom_hilight(node=494, fill=(brewer.pal(10, "Paired"))[3], alpha=.6, extendto = 1.4) + 
  geom_hilight(node=390, fill=(brewer.pal(10, "Paired"))[5], alpha=.6, extendto = 1.31)+ 
  geom_hilight(node=463, fill=(brewer.pal(10, "Paired"))[7], alpha=.6, extendto = 1.34) + 
  geom_hilight(node=457, fill=(brewer.pal(10, "Paired"))[9], alpha=.6, extendto = 1.315)

## Adding Metadata
h2 <- gheatmap(p, completeness, offset = .02,                               # offset shifts the heatmap to the right,
               width = 0.07,                              # width defines the width of the heatmap column,
               # color defines the boarder of the heatmap columns
               colnames = T, colnames_angle=0, font.size=2, color=F, hjust=.5) + theme(legend.position = "bottom", legend.title = element_text(size = 10), legend.text = element_text(size = 5), legend.box = "horizontal", legend.margin = margin()) + scale_fill_viridis_c( direction=-1, name = "Completeness")
h2 <- h2 + new_scale_fill()
h2 <- gheatmap(h2, mags, offset = .17,width = 0.07, colnames = T, colnames_angle=0, font.size=2, color=F, hjust=.5) + scale_fill_manual(values = c(brewer.pal(7, name = "Greys" )[2], brewer.pal(7, name = "Greys" )[7]), name = "Type of Genome")
h2 <- h2 + new_scale_fill()

## To format the CAZyme data:
d2 <- read.csv("path/to/metadata/Cazymes_per_genome_mags.csv")
d2 <- d2[,1:7]
tree$tip.label
order<-match( tree$tip.label, d2$Genome)
d2 <- d2[order,]
d2$Genome == tree$tip.label # check
d2 <- gather(d2, Cazyme, Count, AA:PL)
aa <- d2[d2$Cazyme == "AA",]
ce <- d2[d2$Cazyme == "CE",]
pl <- d2[d2$Cazyme == "PL",]
cbm <-  d2[d2$Cazyme == "CBM",]
gt <-  d2[d2$Cazyme == "GT",]
gh <- d2[d2$Cazyme == "GH",]

## To make proportion values, where each denominator is the max from that column:
aa$Count <- aa$Count/19
ce$Count <- ce$Count/26
pl$Count <- pl$Count/22
cbm$Count<-  cbm$Count/5
gt$Count <-  gt$Count/48
gh$Count <- gh$Count/76
barpal <- park_palette("ArcticGates")

## To make the barplot on the side of the figure:
p4 <- facet_plot(h2, panel = 'Carbohydrate-Binding Modules', data = cbm, 
                 geom = geom_barh, 
                 mapping = aes(x = Count, fill = as.factor(Cazyme)), 
                 stat='identity' , fill = barpal[1], color = "#3F3939") + theme(legend.position = "bottom", legend.text = element_text(size = 7))
p5 <- facet_plot(p4, panel = 'Auxiliary Activites', data = aa, 
                 geom = geom_barh, 
                 mapping = aes(x = Count, fill = as.factor(Cazyme)), 
                 stat='identity' , fill = barpal[2], color = "#3F3939") + theme(legend.position = "bottom", legend.text = element_text(size = 7))
p6 <- facet_plot(p5, panel = 'Polysaccharide Lyases', data = pl, 
                 geom = geom_barh, 
                 mapping = aes(x = Count, fill = as.factor(Cazyme)), 
                 stat='identity' , fill = barpal[3], color = "#3F3939") + theme(legend.position = "bottom", legend.text = element_text(size = 7))
p7 <- facet_plot(p6, panel = 'Carbohydrate Esterase', data = ce, 
                 geom = geom_barh, 
                 mapping = aes(x = Count, fill = as.factor(Cazyme)), 
                 stat='identity' , fill = barpal[4], color = "#3F3939") + theme(legend.position = "bottom", legend.text = element_text(size = 7))
p8 <- facet_plot(p7, panel = 'GlycosylTransferases', data = gt, 
                 geom = geom_barh, 
                 mapping = aes(x = Count, fill = as.factor(Cazyme)), 
                 stat='identity' , fill = barpal[5], color = "#3F3939") + theme(legend.position = "bottom", legend.text = element_text(size = 7))

p9 <- facet_plot(p8, panel = 'Glycoside Hydrolases', data = gh, 
                 geom = geom_barh, 
                 mapping = aes(x = Count, fill = as.factor(Cazyme)), 
                 stat='identity' , fill = barpal[6], color = "#3F3939") + theme(legend.position = "bottom", legend.text = element_text(size = 8))
p9 <- facet_widths(p9, widths = c(7, 1,1,1,1,1,1))
ggsave(paste0("path/to/final/files/Cazymes_All_Heatmap_and_tree2.pdf"), p9, width=20, height = 35, units = "in")

## For the bar plot on the top of the figure, which dictates the values that the proportions are relative to, the following was run:
data <- data.frame(
  name=c("Carbohydrate-Binding Modules","Auxiliary Activites","Polysaccharide Lyases","Carbohydrate Esterase","GlycosylTransferases","Glycoside Hydrolases") ,  
  Maximum=c(5,19,22,26,48,76)
)
plot <- ggplot(data, aes(x=name, y=Maximum)) + 
  geom_bar(stat = "identity", fill = barpal[1:6], color = "#3F3939")
ggsave(paste0("path/to/final/files/Cazymes_plot.pdf"), plot, width=8, height = 6, units = "in")
```

This is the processing that was completed in R for the first figure. All final figure edits were completed in Adobe Illustrator.

## Use Case 2: CAZyme Phylogenies

For this use case, we are using the 78 isolate genomes. For these genomes, a new pangenome was built and new phylogeny was determined. Conveniently, the pangenome work that we have done started with this dataset as version 1, so most of this analysis is completed. First, the genomes storage and pangenome were built:

	#generate genomes storage
	anvi-gen-genomes-storage -e external-genomes-v1.txt -o Alteromonas-v1-GENOMES.db
	#make pangenome
	anvi-pan-genome -g Alteromonas-v1-GENOMES.db -n Alteromonas_Pan --use-ncbi-blast
	#for version naming
	mv Alteromonas_Pan Alteromonas_Pan_v1
	
Then, core gene clusters were identified:

	anvi-script-compute-bayesian-pan-core -p Alteromonas_Pan_v1/Alteromonas_Pan-PAN.db -g Alteromonas-v1-GENOMES.db --store-in-db
	
Then, we can follow the anvi'o's method for finding a set of 'better' genes for phylogenomic analysis. This is explained more in depth [here](https://merenlab.org/2017/06/07/phylogenomics/), but briefly, we want to find a set of genes that has sufficiently evolved within all of the genomes. Moreover, we want these genes to be free of any large insertions or deletions such that they can be easily aligned. For that, we use the `--max-functional-homogeneity-index` and `--min-geometric-homogeneity-index` filters, respectively. We identify these genes and concatenate and align in one step:
	
	anvi-get-sequences-for-gene-clusters -p Alteromonas_Pangenome_v1/Alteromonas_Pan-PAN.db -g Alteromonas-GENOMES-v1.db --max-num-genes-from-each-genome 1 --min-num-genomes-gene-cluster-occurs 78 --max-functional-homogeneity-index .9 --min-geometric-homogeneity-index .925 --concatenate-gene-clusters -o bc_concatenated-proteins.fa

This should result in 111 core gene clusters, or core genes in our case. After producing these concatenated genes, the trees were built in RAxML version 8.2.12:

	./standard-RAxML/raxmlHPC-SSE3 -s pathway/to/bc_concatenated-proteins.fa -m PROTGAMMAAUTO -p 12345 -n TREE1
	
The filter `-p` specifies a random number for the parsimony inferences. The `-m PROTGAMMAAUTO` is used to choose the model for the dataset based on the model that scores the highest likelihood score on the initial parsimony tree. The `-n` filter specifies the name of the output files. Using the output file with the `bestTree` assignment (calculated as the best tree through maximum likelihood analysis), we display the tree in FigTree version 1.4.4, which can be downloaded [here](https://github.com/rambaut/figtree/releases). In FigTree, we can midpoint root the tree because this dataset does not contain an outgroup (common outgroups for *Alteromonas* are *Pseudoalteromonas atlantica* strains) but does require rooting. Once the tree is midpoint rooted, it is exported in newick file format as currently displayed. 

### Getting CAZyme Gene Phylogenies

From the pangenome made for the isolate genomes, it is possible to undergo the same process as completed above for making a split pangenome that only contains gene clusters that have known functions:

	anvi-split -p Alteromonas_Pangenome_v1/Alteromonas\_Pan-PAN.db -g Alteromonas-v1-GENOMES.db -C default --bin-id Bin_1 -o CAZymes
	
To make phylogeny from any one of these genes, it is necessary to identify single-copy core genes, and thus we searched this split pangenome for nay genes that were identified as `SCG` by anvi'o and have the `CORE` distinction assigned by Bayesian statistics. This resulted in 8 unique CAZyme genes overall:

	AA6
	CE1
	CE3
	CE4
	GH103
	GH77
	GT51
	GT81
	
Then, the gene cluster IDs were retrieved for these gene clusters utilizing the summary table from `anvi-summarize`. The IDs were pasted into .txt files for each of the gene clusters. Alignments of each of these genes were made for the genomes:

	anvi-get-sequences-for-gene-clusters p Alteromonas_Pangenome_v1/Alteromonas_Pan-PAN.db -g Alteromonas-v1-GENOMES.db --gene-cluster-ids-file GENE\_CLUSTER\_ID\_FILE\_.txt -o CAZYME-IDENTIFIER-concatenated-proteins.fa

For each of these alignment files, a tree was made using RAxML version 8.2.12:

	./standard-RAxML/raxmlHPC-SSE3 -s pathway/to/CAZYME-IDENTIFIER-concatenated-proteins.fa -m PROTGAMMAAUTO -p 12345 -n TREE-CAZYME-ID
	
In FigTree, we midpoint rooted the trees because this dataset does not contain an outgroup (common outgroups for *Alteromonas* are *Pseudoalteromonas atlantica* strains) but does require rooting. Once the trees are midpoint rooted, they are exported in newick file format as currently displayed. 

### Processing

All of the trees were viewed and two were chosen that roughly followed the phylogeny using the 'better core,' one was chosen that strays from the accepted phylogeny significantly. These are trees made from the CAZymes GT51 and CE3, and GH103 respectively. To visualize these four trees (one accepted and 3 from CAZymes), R was used. This processing also requires a `view.txt` file for processing the proper names of the genomes on the tree.

```r
## load libraries
library("plyr")
library("BiocManager")
library("ggtree")
library("ggplot2")
library("dplyr")
library(stringr)
library("aplot")

## import newick file for the CE3 CAZyme tree:
tree <- read.tree("path/to/raw/newick/files/ce3_complete")
# get the labels on the tree properly
clade_position <- data.frame(tree$tip.label,
                             c(1:78)
)
colnames(clade_position)<-c("label", "order")
genomes <- read.table("path/to/raw/pangenome/files/view.txt", sep ="\t", header=T
)
new_order<-match(clade_position$label,genomes$genome_id)
label2_names <- genomes[new_order,2]
d <- data.frame(label = tree$tip.label, label2=label2_names)
tree2 <- full_join(tree, d, by = "label")
test <- ggtree(tree2)
ceorder <- get_taxa_name(test)
p<-ggtree(tree2) +  theme_tree2() + geom_tiplab(aes(label=label2),size=5,align=TRUE) + ggtitle("Complete Genomes Tree: CE3") +theme(legend.position = "none", axis.title.y = element_blank(), plot.title = element_text(size = 12, face = "bold", hjust = 0.5, vjust = 1)) + hexpand(.43)
ggsave(paste0("path/to/trees/Cazymes_ce3_complete.pdf"), p, width=20, height = 20, units = "in")

## import newick file for the GH103 CAZyme tree:
tree <- read.tree("path/to/raw/newick/files/gh103_complete")
# get the labels on the tree properly
clade_position <- data.frame(tree$tip.label,
                             c(1:78)
)
colnames(clade_position)<-c("label", "order")
genomes <- read.table("path/to/raw/pangenome/files/view.txt", sep ="\t", header=T
)
new_order<-match(clade_position$label,genomes$genome_id)
label2_names <- genomes[new_order,2]
d <- data.frame(label = tree$tip.label, label2=label2_names)
tree2 <- full_join(tree, d, by = "label")
test <- ggtree(tree2)
ghorder <- get_taxa_name(test)
p<-ggtree(tree2) +  theme_tree2() + geom_tiplab(aes(label=label2),size=5,align=TRUE) + ggtitle("Complete Genomes Tree: GH103") +theme(legend.position = "none", axis.title.y = element_blank(), plot.title = element_text(size = 12, face = "bold", hjust = 0.5, vjust = 1)) + hexpand(.43)
ggsave(paste0("path/to/trees/Cazymes_gh103_complete.pdf"), p, width=20, height = 20, units = "in")

## import newick file for the GT51 CAZyme tree:
tree <- read.tree("path/to/raw/newick/files/gt51_complete")
# get the labels on the tree properly
clade_position <- data.frame(tree$tip.label,
                             c(1:78)
)
colnames(clade_position)<-c("label", "order")
genomes <- read.table("path/to/raw/pangenome/files/view.txt", sep ="\t", header=T
)
new_order<-match(clade_position$label,genomes$genome_id)
label2_names <- genomes[new_order,2]
d <- data.frame(label = tree$tip.label, label2=label2_names)
tree2 <- full_join(tree, d, by = "label")
test <- ggtree(tree2)
gtorder <- get_taxa_name(test)
p<-ggtree(tree2) +  theme_tree2() + geom_tiplab(aes(label=label2),size=5,align=TRUE) + ggtitle("Complete Genomes Tree: GT51") +theme(legend.position = "none", axis.title.y = element_blank(), plot.title = element_text(size = 12, face = "bold", hjust = 0.5, vjust = 1)) + hexpand(.43)
ggsave(paste0("path/to/trees/Cazymes_gt51_complete.pdf"), p, width=20, height = 20, units = "in")

## import newick file for the accepted phylogeny:
tree <- read.tree("path/to/raw/newick/files/complete_bc")
# get the labels on the tree properly
clade_position <- data.frame(tree$tip.label,
                             c(1:78)
)
colnames(clade_position)<-c("label", "order")
genomes <- read.table("path/to/raw/pangenome/files/view.txt", sep ="\t", header=T
)
new_order<-match(clade_position$label,genomes$genome_id)
label2_names <- genomes[new_order,2]
d <- data.frame(label = tree$tip.label, label2=label2_names)
tree2 <- full_join(tree, d, by = "label")
test <- ggtree(tree2)
comporder <- get_taxa_name(test)
p<-ggtree(tree2) +  theme_tree2() + geom_tiplab(aes(label=label2),size=5,align=TRUE) + ggtitle("Complete Genomes Tree: Better Core") +theme(legend.position = "none", axis.title.y = element_blank(),plot.title = element_text(size = 12, face = "bold",hjust = 0.5,vjust = 1)) + hexpand(.43)
ggsave(paste0("path/to/trees/Cazymes_compelte_bc.pdf"), p, width=20, height = 20, units = "in")
```

From these images, the order of specific clades were noted (see `Cazyme_alluvial_orders.csv`) from the bottom of the tree to the top in order to make the final alluvial plot, which displays movements of the clades throughout the various trees. This was processed in R:

```r
library(easyalluvial)
library(colorBlindness)
## Import orders file
ordersnew <- read.csv("path/to/orders/file/Cazyme_alluvial_orders.csv", row.names = 1)
ordersnew$comporder <- as.factor(ordersnew$comporder)
ordersnew$ceorder <- as.factor(ordersnew$ceorder)
ordersnew$ghorder <- as.factor(ordersnew$ghorder)
ordersnew$gtorder <- as.factor(ordersnew$gtorder)
## Make figure
col_vector <- (c(paletteMartin[c(1:4,4,4,5:15)], "blanchedalmond"))
al <- alluvial_wide( select(ordersnew, comporder, gtorder, ceorder, ghorder),col_vector_flow = col_vector, fill_by = 'first_variable' , stratum_labels=F) +theme_bw() 
ggsave("path/to/final/figures/CAZymeAlluvial4.pdf", al)
```

Final combination of the original phylogeny (complete\_bc\_) and the alluvial plot in addition to color adjustments were made in Adobe Illustrator. Combination of 4 phylogenies was completed in Illustrator.