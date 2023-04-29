# Alteromonas Genus Pangenome, version 2.0.0

**Date:** February 2023
**Primary Author:** Michelle DeMers

In this edition of the pangenome, we focus on adding metagenome-assembled genomes (MAGs) to our dataset. Most of the steps of this pipeline are largely identical to that of v.1, but we combine all forms of annotation and add the phylogeny in one pipeline. In order to reproduce this pipeline, anvio-7.1, anvio-dev, RAxML version 8.2.12, **and** FigTree version 1.4.4 are required. Here, it is important to note that anvio-dev is the only version of anvi'o that is really necessary, but anvio-dev was only used for KOfam and CAZyme annotation in this workflow. The anvi'o version information is below:

	Anvi'o .......................................: hope (v7.1)
	
	Profile database .............................: 38
	Contigs database .............................: 20
	Pan database .................................: 15
	Genome data storage ..........................: 7
	Auxiliary data storage .......................: 2
	Structure database ...........................: 2
	Metabolic modules database ...................: 2
	tRNA-seq database ............................: 2

And:

	Anvi'o .......................................: hope (v7.1-dev)
	
	Profile database .............................: 38
	Contigs database .............................: 20
	Pan database .................................: 16
	Genome data storage ..........................: 7
	Auxiliary data storage .......................: 2
	Structure database ...........................: 2
	Metabolic modules database ...................: 4
	tRNA-seq database ............................: 2

## Downloading Isolate Genomes

The isolate genomes for this pangenome were downloaded from two locations: [JGI's IMG](https://img.jgi.doe.gov/cgi-bin/m/main.cgi) and [NCBI](https://www.ncbi.nlm.nih.gov/genome/microbes/). 

### IMG Genomes

For the IMG genomes, genomes were located by searching for the term 'Alteromonas' through the 'Find Genomes > Genome Search' tool. All search results that specified a genus of 'Alteromonas' were retained. Any of the search results that were identified as a MAG were excluded from this list. This identified 67 genomes. All genomes were exported via IMG. The IMG genome download produced the nucleotide sequences for the genomes in assembled contigs. All nucleotide sequences from IMG were formatted properly by removing spaces in the headers of the files and from there utilizing `anvi-script-reformat-fasta` to format them in a way that anvi'o would be able to understand:

	#reformat IMG genomes
	sed "s/[[:space:]]/_/g" NEW_GENOME_FILE.fna > INTERMEDIATE_GENOME_FILE.fna
	anvi-script-reformat-fasta NTERMEDIATE_GENOME_FILE.fna -o FINAL_GENOME_FILE.fa -l 0 --simplify-names --seq-type NT
	
These genomes were retrieved and are current as of July 15, 2022.

### NCBI Genomes

To cover all of our bases, we also downloaded any complete genomes from NCBI with the commands:

	# Add NCBI genomes
	ncbi-genome-download bacteria --assembly-level complete --genera Alteromonas --metadata NCBI-alter-metadata.txt -F fasta
	gunzip *.gz

These genomes were retrieved and are current as of November 8, 2022. This produced 40 files, but after investigating the meta-data file, it became clear that only 11 of them were new (unique) compared to the 67 that were already downloaded from IMG. These 11 unzipped fasta files were then reformatted in the same manner that the IMG genomes were reformatted above.

## Downloading MAGs

### Paoli et al., 2022

These genomes are from a variety of original sources, but were nicely compiled in [this](https://www.nature.com/articles/s41586-022-04862-3) paper from Paoli et al. [1]. Specifically, this paper compiles MAGs made from differently size-fractionated metagenomes isolated from Tara Oceans, Malaspina, BioGEOTRACES, HOT, and BATS [1]. From this paper, there were 245 MAGs identified as a part of the _Alteromonas_ genus, and the genomes were turned into contigs.dbs and annotated with NCBI COGs by the Meren group (A. Murat Eren and Matthew Schechter).

### IMG MAGs

For the IMG MAGs, genomes were located by searching for the term 'Alteromonas' through the 'Find Genomes > Genome Search' tool. All search results that specified a genus of 'Alteromonas' were retained. Any of the search results that were identified as isolates were excluded from this list. This identified 18 genomes from two projects: 'Water column microbial communities from Red Sea, Saudi Arabia' and 'Draft Metagenome-assembled Genomes from Global Oceans - TOBG MAGs.' All genomes were exported via IMG. The IMG genome download produced the nucleotide sequences for the genomes in assembled contigs. All nucleotide sequences from IMG were formatted properly by removing spaces in the headers of the files and from there utilizing `anvi-script-reformat-fasta` to format them in a way that anvi'o would be able to understand:

	#reformat IMG genomes
	sed "s/[[:space:]]/_/g" NEW_GENOME_FILE.fna > INTERMEDIATE_GENOME_FILE.fna
	anvi-script-reformat-fasta NTERMEDIATE_GENOME_FILE.fna -o FINAL_GENOME_FILE.fa -l 0 --simplify-names --seq-type NT
	
These genomes were retrieved and are current as of July 15, 2022.

## Making Contigs Databases

Contigs databases were then made from all of the files downloaded from IMG, and the databases were hmm-searched using the following commands: 
	
	#in the directory with all of the IMG and NCBI reformatted, "FINAL" fa files, we run the command: 
	for i in `ls *fa | awk 'BEGIN{FS=".fa"}{print $1}'`
	do
		anvi-gen-contigs-database -f $i.fa -o $i.db -T 4
    	anvi-run-hmms -c $i.db
	done

## Annotating Genomes

Our next goal was to annotate the genomes to the best of our ability. The annotations that were added were from COGs, KOfam, dbCAN (CAZymes), Eggnog, and GhostKOALA. We expect from these annotations to produce similar results between the Kofam and GhostKOALA annotations, but we may see slightly different results because they use different search and determination methods to determine if a gene is homologous to a KEGG ID or not. In future versions, the goal would be to incorporate annotations form interproscan, but this source of annotation has not been a top priority thus far. NOTE: The contigs.db from the Meren group were already annotated with COGs, so that portion of the pipeline was not run on those contigs.db.

### COGs

First, the COGs database must be set up within anvi'o [2]: 
	
	anvi-setup-ncbi-cogs --just-do-it
 
To annotate the IMG genomes with COGs, we can conveniently use an anvi'o command, run on each genome separately:

	anvi-run-ncbi-cogs -c ${EACH GENOME}.db --num-threads 20

### KOfam

Within anvio-dev (the developer model of anvi'o), we can annotate the genomes with KOfams, which are HMMs built from the KEGG database [3]. These hmms are searched for within each contigs database, and then a KEGG ID (KOfam), KEGG Module, KEGG Class, and BRITE Hierarchy (developer version only) are assigned to each gene call in the database if the score from the search meets the kofam criteria. To annotate with KOfam, the following is run on each genome separately:

	anvi-run-kegg-kofams -c ${EACH GENOME}.db --just-do-it

### dbCAN (CAZymes)

The anvi'o implementation of this annotation source was written by Matthew Schechter [4] and utilizes dbCAN [5], a database that contains HMMs of each family and subfamily of the CAZymes database. Within anvio-dev, we run the CAZyme hmm search on the genomes in our database, on each contigs.db:

	#set up the cazymes database within the anvi'o environment
	anvi-setup-cazymes
	anvi-run-cazymes -c alteromonas_genomes_databases_v1/${EACH GENOME}.db

### Eggnog

Eggnog is a nice annotation source to use here because it incorporates information from a number of databases into its annotation sources. To annotate with eggnog, there is a little more work that must be done for each genome. The eggnog-mapper source code can be downloaded from [here](https://github.com/eggnogdb/eggnog-mapper). Once installed, the following can be run:
	
	#first the amino acid sequences for each gene call in each genome must be retrieved
	anvi-get-sequences-for-gene-calls -c ${EACH GENOME}.db -o path/to/aminoAcids/protein_sequences_${EACH GENOME}.fa --get-aa-sequences
	
	#the amino acid sequences for each genome can then be run through eggnogmapper:
	path/to/eggnog-mapper/emapper.py -i path/to/aminoAcids/protein_sequences_${EACH GENOME}.fa -o path/to/eggnogResults/${EACH GENOME} --override
	
	#Then, once each genome has been annotated, the annotations must be added to the contigs.db, whichcan conveniently be done with
	anvi-script-run-eggnog-mapper -c ${EACH GENOME}.db --annotation path/to/eggnogResults/${EACH GENOME}.emapper.annotations --use-version 1.0.3

Functional annotation was performed using eggNOG-mapper (version emapper-2.1.9) [6]
 based on eggNOG orthology data [7]. Sequence searches were performed using [8].

### GhostKOALA

Here, annotations with GhostKOALA were also added in the hopes that it would add to the KEGG ID coverage that KOfams did not provide, though we expect this output to be relatively similar. Here, all the amino acid sequences were concatenated into one file, with each gene call identifier containing an indication of which genome it came from:
	
	#add genome labels to the files (run for each genome)
	sed '/>/s/$/${EACH GENOME}/' path/to/aminoAcids/protein_sequences_${EACH GENOME}.fa > test.tmp && mv test.tmp path/to/aminoAcids/protein_sequences_${EACH GENOME}.fa
	#concatenate
	cat path/to/aminoAcids/protein_sequences* > all_genes.fa

This `all_genes.fa` file was then fed into the GhostKOALA web server. Once the GhostKOALA results (named `GhostKoalaAll.txt`) came back, the file was split into files based on the genome from which it came, the results were converted to a form that anvi'o could parse, and then the functions (annotations) were added to the contigs databases, one by one. The python script for parsing the GhostKOALA output was written by Elaina Graham, and documentation on this full process can be found [here](https://merenlab.org/2018/01/17/importing-ghostkoala-annotations/).

	#split GhostKOALA results
	grep "${EACH GENOME}" GhostKoalaAll.txt > Ghost_${EACH GENOME}.txt
	#parse
	python GhostKoalaParser/KEGG-to-anvio --KeggDB KO_Orthology_ko00001.txt -i Ghost_${EACH GENOME}.txt -o AnviImportable_${EACH GENOME}.txt
	#import into anvi'o
	anvi-import-functions -c ${EACH GENOME}.db -i AnviImportable_${EACH GENOME}.txt

## Building the Pangenome

After all the contigs databases are in forms that we are happy with (n=336), we can create the pangenome. This involves creating genomes storage for the set of genomes that we are working with, and from there creating the pangenome.

	#generate genomes storage
	anvi-gen-genomes-storage -e external-genomes-v2.txt -o Alteromonas-v2-GENOMES.db
	#make pangenome
	anvi-pan-genome -g Alteromonas-v2-GENOMES.db -n Alteromonas_Pangenome_2 --num-threads 10 --use-ncbi-blast

Here, it is important to note that we removed 5 MAGs from Paoli et al. when building the genomes storage and pangenome due to issues with the way that anvi'o read them: BGEO\_SAMN07136708\_METAG\_BAHNCNMJ-contigs.db, MALA\_SAMN05421660\_METAG\_DPEMGBDK-contigs.db, TARA\_SAMEA2619974\_METAG\_DFCIFEPC-contigs.db, TARA\_SAMEA2619974\_METAG\_MIENHEPK-contigs.db, and TARA\_SAMEA2622197\_METAG\_COANOGOL-contigs.db. **This brings our number of genomes in the pangenome down from 341 to 336.**

## Bayesian Statistics (mOTUpan)

Once the pangenome was created, we became interested in identifying the core gene clusters within the pangenome, and so we utilized the bayesian statistical method developed by Buck et al. that utilizes mOTUpan.py to determine what gene clusters are likely to be core gene clusters (i.e. found within all of the genomes) based on individual genome completeness [10]. To compute this, the following is run: 
	
	anvi-script-compute-bayesian-pan-core -p Alteromonas_Pangenome_2/Alteromonas_Pangenome_2-PAN.db -g Alteromonas-v2-GENOMES.db --store-in-db

If it is desired to produce a list of gene cluser ids that are determined to be core as opposed to accessory, an out file could be specified with `-o` instead of using `--store-in-db`. ***I*** think that the out file is a bit difficult to process compared to the anvi'o `misc_data_items` directory that is gathered from `anvi-summarize`, so I prefer storing mOTUpan results in the database and using `anvi-summarize` to get a gene-clusters-type summary.

## Computing Genome Similarity

Computing genome similarity (ANI) can be accomplished within anvi'o, and we run these commands in the interest of observing how ANI correlates with phylogeny.

	anvi-compute-genome-similarity --external-genomes external-genomes-v2.txt --program pyANI --output-dir ANI --num-threads 6 --pan-db Alteromonas_Pangenome_2/Alteromonas_Pangenome_2-PAN.db
	
## Adding Phylogeny to the Pangenome

### Identification of Core Genes

There are many ways to visualize the phylogeny of genomes in a pangenome, and so we introduce a number of methods below. In order to build phylogeny in a pangenome, genes must be identified that will be used to assume phylogenetic relationships. There are multiple methods that we can use to determine core genes from our pangenome, and so we'll briefly describe three different versions here before we move onto choosing one that best fits our data, and then incorporating it into our pangenome.

### Bayesian Core

Above we added an assignment of `gene-cluster-type` to all of the gene clusters in the pangenome using bayesian statistics, and then we were able to get a list of the gene cluster ids for all of the 2,354 'CORE' gene clusters. We call this the bayesian core. We format this list of core gene clusters in the file `bayesian_core_gene_list.txt`:

	GC_00000001
	GC_00000002
	GC_00000003
	GC_00000004
	GC_00000006
	GC_00000007
	GC_00000008
	GC_00000009
	GC_00000011
	GC_00000012
	(...)

Then, we obtain the concatenated gene sequences from these gene clusters via an anvi'o command:

	anvi-get-sequences-for-gene-clusters -p Alteromonas_Pangenome_2/Alteromonas_Pangenome_2-PAN.db -g Alteromonas-v2-GENOMES.db -gene-cluster-ids-file bayesian_core_gene_list.txt -o bayesian_MAGs_concatenated_proteins.fa --max-num-genes-from-each-genome 1 --concatenate-gene-clusters

These gene clusters must be concatenated and aligned for downstream analysis, so we use the `--concatenate-gene-clusters` flag. We also must make sure that these use single copy genes, so we include the `--max-num-genes-from-each-genome 1` filter, which brings the number of genes down from 2,354 to 110 genes.

### Bacterial Core

We were also interested in how utilizing universal bacterial core genes in our phylogenomic analyses would affect our results. Here, we utilize a set of hmms that are built into anvi'o, `Bacteria_71`, that were searched when we ran the command, `anvi-run-hmms`. To align and concatenate these genes for use in building phylogeny, we run:

	anvi-get-sequences-for-hmm-hits --external-genomes external-genomes-v2.txt -o bacterial-concatenated-proteins.fa --hmm-source Bacteria_71 --return-best-hit --get-aa-sequences --concatenate
	
### Better Core

Lastly, we can follow anvi'o's method for finding a set of 'better' genes for phylogenomic analysis. This is explained more in depth [here](https://merenlab.org/2017/06/07/phylogenomics/), but briefly, we want to find a set of genes that has sufficiently evolved within all of the genomes. Moreover, we want these genes to be free of any large insertions or deletions such that they can be easily aligned. For that, we use the `--max-functional-homogeneity-index` and `--min-geometric-homogeneity-index` filters, respectively. We identify these genes and concatenate and align in one step:
	
	anvi-get-sequences-for-gene-clusters -p Alteromonas_Pangenome_2/Alteromonas_Pangenome_2-PAN.db -g Alteromonas-v2-GENOMES.db --max-num-genes-from-each-genome 1 --min-num-genomes-gene-cluster-occurs 170 --max-functional-homogeneity-index .95 --min-geometric-homogeneity-index .85 --concatenate-gene-clusters -o better-core-mags-concatenated-proteins.fa

This should result in 80 core gene clusters, or core genes in our case. Here, it is important to note that the filter `--min-num-genomes-gene-cluster-occurs 170` had to be placed on our data because the MAGs severely limited gene clusters that were universal across all genomes in the pangenome. With this in place, this set of genes is not the strongest representation of the core of the pangenome even though it likely has some evolutionary significance for the pangneome.

## Tree Making

We are going to make trees out of all three of the gene sets in the next step. HOWEVER, we are going to use the core genes determined with bayesian statistics (bayesian core) to make the final tree that will be incorporated into the pangenome as we suspect that this is the best set of genes that we can use to make phylogeny. To make the trees, RAxML version 8.2.12 was used [11]. The source code can be downloaded from [here](https://github.com/stamatak/standard-RAxML). Once downloaded, the following commands were run in the same directory as the source code:

	#create bayesian tree
	./standard-RAxML/raxmlHPC-SSE3 -s pathway/to/bayesian_MAGs_concatenated_proteins.fa -m PROTGAMMAAUTO -p 12345 -n bayesian
	#create raw core tree
	./standard-RAxML/raxmlHPC-SSE3 -s pathway/to/bacterial-concatenated-proteins.fa -m PROTGAMMAAUTO -p 12345 -n bacterial
	#create better core tree
	./standard-RAxML/raxmlHPC-SSE3 -s pathway/to/better-core-mags-concatenated-proteins.fa -m PROTGAMMAAUTO -p 12345 -n bc

The filter `-p` specifies a random number for the parsimony inferences. The `-m PROTGAMMAAUTO` is used to choose the model for the dataset based on the model that scores the highest likelihood score on the initial parsimony tree. The `-n` filter specifies the name of the output files (e.g. *.bayesian, *.bacterial, or *.bc).

Using the output file with the `bestTree` assignment for each dataset (calculated as the best tree through maximum likelihood analysis), we display the tree in FigTree version 1.4.4, which can be downloaded [here](https://github.com/rambaut/figtree/releases). In FigTree, we can midpoint root the tree because this dataset does not contain an outgroup (common outgroups for *Alteromonas* are *Pseudoalteromonas atlantica* strains) but does require rooting. Once the trees are midpoint rooted, they are exported in newick file format as currently displayed. These exported trees are named: `midpoint_mags_bayesian`, `midpoint_mags_bacterial` and `midpoint_mags_better_core` for the bayesian core, bacterial core, and better core trees respectively.

## Incorporation of Phylogeny into Pangenome

Once again, the tree that was chosen as the final phylogeny was the tree made from the bayesian core set of genes. All three of the trees are very similar in clade structure, so this decision was made purely based on the reproducibility and non-arbitrary nature of the gene cluster sourcing procedure. This midpoint rooted tree, `midpoint_mags_bayesian`, is the file that will be used in subsequent work.

To incorporate a tree into a pangenome, a data layers file must be imported into anvi'o in a format that anvi'o can read. This is in the form of a tab-delimited file with three columns: `item\_name`, `data\_type`, and `data\_value`. The `item\_name` column contains the name of the tree, which we call `Phylogenomics_from_Bayesian_Core`. The `data\_type` column contains the kind of data that we are importing, which is `newick`. The `data\_value` column contains the newick tree that we choose. Once this three-column tab-delimited file, named `layer_orders.txt`, is created, we can import it into anvi'o:

	anvi-import-misc-data -p Alteromonas_Pangenome_2/Alteromonas_Pangenome_2-PAN.db -t layer_orders layer_orders.txt
	
To view the pangenome in the order of the phylogeny that we added to the database, we can use the anvi'o interactive interface. First we run:
	
	anvi-display-pan -p Alteromonas_Pangenome_2/Alteromonas_Pangenome_2-PAN.db -g Alteromonas-v2-GENOMES.db 
	
The pangenome will be displayed in the last state that was saved, which will not be in the order specified by phylogeny. Thus, we must redisplay the order of genomes by going to Settings > Layers > Order by > Phylogenomics\_from\_Bayesian\_Core, and then redrawing the Pangenome. In our case, it is also of interest to display how ANI correlates with the tree that we built as we would assume that they would be highly related. To display ANI in addition the phylogeny already drawn on the tree, we navigate to Settings > Layers > Order and select the 'Layer Groups' as 'ANI\_percentage\_identity,' and then we redraw the tree. This state of the pangenome can be saved by navigating to the 'Save State' feature of the anvi'o interactive interface.

## Uploading to Zenodo

The contigs databases, genomes storage, genomes information, pangenome files, and phylogeny files are uploaded as version 2.0.0 to Zenodo for community use. Please see the [C-CoMP Zenodo Community page](https://zenodo.org/communities/c-comp/). For questions about this version or future work, please reach out.

## References

[1] Paoli, L., Ruscheweyh, HJ., Forneris, C.C. _et al_. Biosynthetic potential of the global ocean microbiome. _Nature_ **607**, 111–118 (2022). https://doi.org/10.1038/s41586-022-04862-3

[2] Galperin MY, Makarova KS, Wolf YI, Koonin EV. Expanded microbial genome coverage and improved protein family annotation in the COG database. Nucleic Acids Res. 2015 Jan;43(Database issue):D261-9. doi: 10.1093/nar/gku1223. Epub 2014 Nov 26. PMID: 25428365; PMCID: PMC4383993.

[3] Aramaki T., Blanc-Mathieu R., Endo H., Ohkubo K., Kanehisa M., Goto S., Ogata H. KofamKOALA: KEGG ortholog assignment based on profile HMM and adaptive score threshold. Bioinformatics. 2019 Nov 19. pii: btz859. doi: 10.1093/bioinformatics/btz859.

[4] Schechter, M, anvi-run-cazymes. https://anvio.org/help/main/programs/anvi-run-cazymes/

[5] Yin Y, Mao X, Yang JC, Chen X, Mao F and Xu Y, dbCAN: a web resource for automated carbohydrate-active enzyme annotation, Nucleic Acids Res. 2012 Jul;40(Web Server issue):W445-51

[6] eggNOG-mapper v2: functional annotation, orthology assignments, and domain prediction at the metagenomic scale. Carlos P. Cantalapiedra, Ana Hernandez-Plaza, Ivica Letunic, Peer Bork, Jaime Huerta-Cepas. 2021. Molecular Biology and Evolution, msab293, https://doi.org/10.1093/molbev/msab293

[7] eggNOG 5.0: a hierarchical, functionally and phylogenetically annotated orthology resource based on 5090 organisms and 2502 viruses. Jaime Huerta-Cepas, Damian Szklarczyk, Davide Heller, Ana Hernandez-Plaza, Sofia K Forslund, Helen Cook, Daniel R Mende, Ivica Letunic, Thomas Rattei, Lars J Jensen, Christian von Mering and Peer Bork. Nucleic Acids Research, Volume 47, Issue D1, 8 January 2019, Pages D309-D314, https://doi.org/10.1093/nar/gky1085 

[8] Sensitive protein alignments at tree-of-life scale using DIAMOND. Buchfink B, Reuter K, Drost HG. 2021. Nature Methods 18, 366–368 (2021). https://doi.org/10.1038/s41592-021-01101-x

[9] Kanehisa, M., Sato, Y., and Morishima, K. (2016) BlastKOALA and GhostKOALA: KEGG tools for functional characterization of genome and metagenome sequences. J. Mol. Biol. 428, 726-731.

[10] mOTUpan: a robust Bayesian approach to leverage metagenome assembled genomes for core-genome estimation Moritz Buck, Maliheh Mehrshad, and Stefan Bertilsson bioRxiv 2021.06.25.449606; doi: https://doi.org/10.1101/2021.06.25.449606

[11] Stamatakis, A. "RAxML Version 8: A tool for Phylogenetic Analysis and Post-Analysis of Large Phylogenies". In Bioinformatics, 2014, open access.