# Making the _Ruegeria pomeroyi_ Digital Microbe

**Author: Iva Veseli**

**Date: August 2022**
**(Updated: May 2023)**


This workflow describes the generation of a Digital Microbe for _Ruegeria pomeroyi_ DSS-3, a model marine organism. The Moran lab has been curating data on this organism for years, and there are several datasets to be integrated via the framework we have developed using the anvi'o platform. 


We will start with generating the contigs database and importing a curated set of gene annotations. Zac Cooper and Christa Smith have sent me an excel spreadsheet of their annotations (`Ruegeria pomeroyi DSS-3 gene annotations - March 2021.xlsx`), as well as an external gene calls file (`DSS3_external_gene_calls.txt`).

Note that I’m using anvi’o-dev to create the databases. Relevant version info:

```
Anvi'o .......................................: hope (v7.1-dev)

Profile database .............................: 38
Contigs database .............................: 20
Pan database .................................: 16
Genome data storage ..........................: 7
Auxiliary data storage .......................: 2
Structure database ...........................: 2
Metabolic modules database ...................: 4
tRNA-seq database ............................: 2
```

## Making the contigs database
Their data refers to genes with SPOxxxx locus tag numbers, so I made a file mapping SPOXXXX ids to integer gene caller ids, which is called `gcid_to_SPO_locus.txt` and looks like this:

|**gene\_callers_id**|**SPO\_Locus_Tag**|
|:--|:--|
|0|SPO0001|
|1|SPO0002|
|2|SPO0003|
|3|SPO0004|

There are 4448 genes, so the gene caller IDs start from 0 and go all the way to 4447.

I used Python to exchange the SPO ids with the gene caller ids to create a new external gene calls file, `DSS3_external_gene_calls_NEW.txt`, which looks like the following:

|**gene\_callers_id**|**contig**|**start**|**stop**|**direction**|**partial**|**call_type**|**source**|**version**|**aa_sequence**|
|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|
|0|chromosome|0|1869|f|0|1|Glimmer|2.0|VKHSDFDIVVIGAGHAGAEAAHAAARMGMRTALVSLSERDIGVMSCNPAIGGLGKGHLVREIDALDGVMGRVADKAGIQFRLLNRRKGPAVQGPRAQADRKLYRLAMQEEMRNRPGLTIVEGEVTDFRMQGDRVAGVVLADGSEIASQAVILTSGTFLRGIIHIGDVSRPGGRMGDRPSVPLAERLDGFALPMGRLKTGTPPRLDGRTIDWSILERQDGDDDPVLFSFLSKGAYARQIACGITHTNAQTHEIIRKNLSRSAMYGGHIEGVGPRYCPSIEDKIVRFADKDSHQIFLEPEGLEDHTVYPNGISTSLPVDVQEDYVRSIRGLEQVEILQPGYAIEYDYVDPRALTSQLSLPNVPGLYLAGQINGTTGYEEAAAQGMVAGLNAATAILGHEPVPFSRANSYIGVMIDDLTTRGVAEPYRMFTSRAEFRLSLRADNADQRLTPLGLELGCVGDERRDVFARKAEKLATASALLDQSSFSPKEIATAGITISQDGNRRNGFAVLAFPDVRFDDLVPLIPELADTDAETRAQVERDALYANYIARQERDVEAMKRDEALVIPIDFNFSALDGLSNELKQKLTSARPENIAQAGRVEGMTPAALALILARLRRGDRARSA|
|1|chromosome|1865|2480|f|0|1|Glimmer|2.0|MMVPDANTLNVSRETFERLKIFADLVHKWNPRINLVSKRSLEDLWTRHIIDSIQVFRNAPKGDLWVDLGSGGGFPGLICAILAAEEKPETQFICVESDQRKSAFLRSAARECGIACQVISERIEHLDPLDADILSARALTDLTGLLGFAERHLKIGGTALFPKGAAWKKELQDAAKQWNFSYDAVTSLTEPQAVLLKITGVTRV|
|2|chromosome|2472|3276|f|0|1|Glimmer|2.0|VSDLSRPAGPRIIAVANQKGGVGKTTTAINLAAALVESGQRVLVVDLDPQGNASTGLGVDERELTTYELLVDDAPLNSVIQKTSIDGLSIVPATVDLSSADIELISNEKRSFLLHDALRQTAMDAYSWDYILIDCPPSLNLLTVNAMVAAHSVLVPLQSEFFALEGLSQLMLTIREVRQAANPNLRIEGIVLTMYDRRNNLSQQVEKDARDNLGDLVFETKIPRNVRVSEAPSFAMPVLNYDPNSLGAMAYRDLAAELMKKHNKIAA|
|3|chromosome|3298|4189|f|0|1|Glimmer|2.0|MVSNKPRGLGRGLSALMADVTQPAEAAASEAARRPDRTVPIEKLRANPNQPRRTFTEDALQELAASIKEKGVLQPLIVRPVDGDMYEIVAGERRWRAAQLAQLHQVPVLVRELDDTEVLEIAIIENIQRADLNAVEEAAGYRQLMDKFGHTQEKLAEALGKSRSHIANLLRLLSLPDDVQTLVVEGKLSAGHARALITSDNPSELAKIVVRDGLSVRATEALVKKQAEGDKPTATARPKSITADKDADTRALEKDLSAILAMKVAINHKPGTETGQVVLTYENLDQLDDLCAKLSR|

Note for the record: I had to add 1 to each stop position in the file that Zac sent (which had both start and stop positions based on a 0-indexed genome sequence) to make it adhere to the pythonic string slicing convention that we use to indicate gene positions, where the stop position is actually the end of the gene + 1, as described [here](https://anvio.org/help/main/artifacts/external-gene-calls/#gene-startstop-positions).

Then I looked at what the contig names were in the fasta file. There are only two - the main chromosome and the megaplasmid. I changed the headers (the NCBI ones have dots and spaces and dashes which anvi’o won’t accept) to match the simple ‘chromosome’ and ‘megaplasmid’ used in the external gene calls file.

Then I attempted to make the contigs database with the external gene calls file, skipping frame prediction since AA sequences were already provided for all relevant gene calls:

```
anvi-gen-contigs-database -f R-pom_genome.fasta --external-gene-calls DSS3_external_gene_calls_NEW.txt --skip-predict-frame -n R_POMEROYI_DSS3 -o R_POM_DSS3-contigs.db
```

When I was making this database prior to incrementing the gene stop positions, I had gotten the following error about gene length. The error came up because a tRNA gene did not have an AA sequence (as they shouldn’t) but was erroneously tagged as a ‘coding’ gene in the external gene calls file:

```
Config Error: The sequence corresponding to the gene callers id '4067' has 74 nucleotides,
              which is indivisible by 3. This is bad because it is now ambiguous which codon
              frame should be used for translation into an amino acid sequence. Here is the
              culprit sequence:
              GCCCGGGTAGCTCAGGGGTAGAGCAGTGGATTGAAAATCCTCGTGTCGGTGGTTCGATTCCGCCCCCGGGCACC.
              Since you are creating a contigs database, anvi'o is willing to strike you a
              deal, but it will require you to trust her a bit more and give her the power to
              modify the external gene calls you provided. In your external gene calls file
              you have at least one gene call for which you did not provide an amino acid
              sequence, and marked it as a CODING type gene call.
[…..]
```
Here, gcid 4067 corresponds to the sequence with id SPOA_tRNA-Phe-1. There is a 1 in the coding column, and there shouldn’t be. When I checked all tRNA sequences in the external gene calls file, I saw that only the two tRNAs on the megaplasmid were labeled incorrectly (the rest are labeled with a 2 for a noncoding gene). I changed the file for these two (gcids 4067 and 4130) so they have a 2 in the coding column. So, that little mistake in the gene stop positions actually helped me to find and fix these wrongly-tagged genes :)  

After updating the gene calls file, I ran `anvi-gen-contigs-database` again and it worked without errors this time. Here is the parser report:

```
EXTERNAL GENE CALLS PARSER REPORT
===============================================
Num gene calls in file .......................: 4,448
Non-coding gene calls ........................: 63
Partial gene calls ...........................: 0
Num amino acid sequences provided ............: 4,279
  - For complete gene calls ..................: 4,279
  - For partial gene calls ...................: 0
Frames predicted .............................: 0
  - For complete gene calls ..................: 0
  - For partial gene calls ...................: 0
Gene calls marked as NONCODING ...............: 0
  - For complete gene calls ..................: 0
  - For partial gene calls ...................: 0
Gene calls with internal stops ...............: 0
  - For complete gene calls ..................: 0
  - For partial gene calls ...................: 0
```
It correctly identifies the 63 noncoding gene calls at first, but later there is a 0 in the NONCODING category. I suspect this is a reporting bug. I checked the database for the noncoding gene calls and they are in there for sure, and labeled as non-coding.

106 of the coding genes were missing amino acid sequences (there are 4,385 coding gene calls, but only 4,279 AA sequences) in the external gene calls file, but I verified that anvi'o filled in these amino acid sequences after the database was created. 

Here is the creation report (it looks as expected).

```
CONTIGS DB CREATE REPORT
===============================================
Split Length .................................: 20,000
K-mer size ...................................: 4
Skip gene calling? ...........................: False
External gene calls provided? ................: True
External gene calls file have AA sequences? ..: True
Proper frames will be predicted? .............: False
Ignoring internal stop codons? ...............: False
Splitting pays attention to gene calls? ......: True
Contigs with at least one gene call ..........: 2 of 2 (100.0%)
Contigs database .............................: A new database, R_POM_DSS3-contigs.db, has been created.
Number of contigs ............................: 2
Number of splits .............................: 228
Total number of nucleotides ..................: 4,601,048
Gene calling step skipped ....................: False
Splits broke genes (non-mindful mode) ........: False
Desired split length (what the user wanted) ..: 20,000
Average split length (what anvi'o gave back) .: 20,180
```

## Importing annotations

Here is the plan for how to handle each column in the master annotation file (`Ruegeria pomeroyi DSS-3 gene annotations - March 2021.xlsx`). The columns are: locus tag, gene ID, description, RefSeq locus tag, gene ID, and product. In this file the initial locus tag appears to simply be the SPOxxxx ID, and the RefSeq locus tag appears to be equivalent to the NCBI/Genbank locus tag. Similarly, the RefSeq gene ID/product columns appear to be equivalent to the Genbank symbol/name. I will ignore the RefSeq gene ID/product columns and instead focus on the first four columns:

- locus tag (aka SPO ID, one per gene): I will import these as source ’SPO_ID’, with the id itself as the accession and nothing in the ‘function’ column. (all e-values 0)
- GeneID/Description (one per gene): Zac told me these are the important columns as they are the curated annotations. I will import these as source ‘Gene_ID’, with the id itself (often equivalent to SPO id, but with many exceptions) as the accession and the ‘Description’ as the function. (all e-values 0). 
- RefSeq Locus Tag (most genes). I will import these as source ‘Locus_Tag', with the locus tag as the accession and nothing in the ‘function’ column (all e-values 0)

Afterwards, I will run `anvi-run-kegg-kofams` (with BRITE, as Sam requested) and `anvi-run-ncbi-cogs` on the database to get functional annotations from these sources.

Without further ado, here I create the functions-txt files necessary for import (using tab-delimited files saved from the Excel sheets). I will make one separately for each annotation source so that it is easier to diagnose any problems when importing.

#### SPO ID

```
echo -e "gene_callers_id\tsource\taccession\tfunction\te_value" > SPO_external_functions.txt
while read gcid spo; do echo -e "$gcid\tSPO_ID\t$spo\t\t0" >> SPO_external_functions.txt ; done < <(tail -n+2 gcid_to_SPO_locus.txt)
```

#### (RefSeq) Locus Tag + Gene ID/Description

For most genes, I need to match the RefSeq locus tags to gene caller ids via the SPO id. There are also 106 gene sequences that we added to our set later, and these are special cases for which the RefSeq Locus Tag is described in another file, `added_gene_calls.txt`. Since multiple files of information need to be integrated here, I will use Python to do it.

```python
import pandas as pd
spo_to_gcid = {}
with open("gcid_to_SPO_locus.txt", 'r') as f:
     for line in f.readlines()[1:]:
             gcid, spo = line.strip().split("\t")
             spo_to_gcid[spo] = gcid

annotations = pd.read_csv("Ruegeria_pomeroyi_DSS-3_gene_annotations-March_2021.txt", sep="\t", index_col=0)
added = pd.read_csv("added_gene_calls.txt", sep="\t", index_col=0)

locus_df = pd.DataFrame(index=spo_to_gcid.values(), columns=['source', 'accession', 'function', 'e_value'])
locus_df.index.name = 'gene_callers_id'
not_found = []
for spo, gcid in spo_to_gcid.items():
     if spo in added['gene_callers_id'].values:
          locus_df.loc[gcid, 'source'] = 'Locus_Tag'
          locus_df.loc[gcid, 'accession'] = added[added['gene_callers_id'] == spo].index.values[0]
          locus_df.loc[gcid, 'function'] = ''
          locus_df.loc[gcid, 'e_value'] = 0
     elif spo in annotations.index:
          locus_df.loc[gcid, 'source'] = 'Locus_Tag'
          locus_df.loc[gcid, 'accession'] = annotations.loc[spo, 'RefSeq_locus_tag']
          locus_df.loc[gcid, 'function'] = ''
          locus_df.loc[gcid, 'e_value'] = 0
     else:
          not_found.append(spo)
# not found: ['SPOA_tRNA-Phe-1', 'SPOA_tRNA-Thr-1']
locus_df[locus_df.isna().any(axis=1)].shape[0]
# 53 genes have no RefSeq locus tag. Example gcids: 108, 267, 273, 279, etc
locus_df = locus_df.dropna(axis=0, how='any')
locus_df.to_csv("LOCUS_external_functions.txt", sep="\t")

gene_df = pd.DataFrame(index=spo_to_gcid.values(), columns=['source', 'accession', 'function', 'e_value'])
gene_df.index.name = 'gene_callers_id'
not_found = []
for spo, gcid in spo_to_gcid.items():
     if spo in annotations.index:
          gene_df.loc[gcid, 'source'] = 'Gene_ID'
          gene_df.loc[gcid, 'accession'] = annotations.loc[spo, 'gene_ID']
          gene_df.loc[gcid, 'function'] = annotations.loc[spo, 'description']
          gene_df.loc[gcid, 'e_value'] = 0
     else:
          not_found.append(spo)
# not found again contains: 'SPOA_tRNA-Phe-1', and 'SPOA_tRNA-Thr-1' as well as all 106 of the genes we added later (as expected)
gene_df[gene_df.isna().all(axis=1)]
# correspondingly, there are 108 gcids with all NAN in this dataframe, for the 106 added-later genes plus the two tRNAs that were not found in the annotations file

all_nans = gene_df[gene_df.isna().all(axis=1)].index
not_nans = gene_df.index.difference(all_nans)
gene_df = gene_df.loc[not_nans] # drop rows with all NA (for some reason dropna(how='all', axis=1, inplace=True) was not working)
# replace other NA with empty string
gene_df.fillna('', inplace=True)

# add back the two missing tRNAs
gene_df.loc[spo_to_gcid['SPOA_tRNA-Phe-1']] = {'source': 'Gene_ID', 'accession': '', 'function': 'tRNA-Phe', 'e_value': 0}
gene_df.loc[spo_to_gcid['SPOA_tRNA-Thr-1']] = {'source': 'Gene_ID', 'accession': '', 'function': 'tRNA-Thr', 'e_value': 0}

gene_df.to_csv("GENE_ID_external_functions.txt", sep="\t")

# there are 4448 genes w/ SPO ids
# 4395 have a locus tag according to the annotations file (missing 53 genes)
# 4342 have a gene ID/description (so only the 106 genes that were added later are missing descriptions)
```

There are two tRNA genes that were present in the external gene calls file but not in either annotation file: SPOA_tRNA-Phe-1 and SPOA_tRNA-Thr-1. I suspect this is because these are the two tRNAs from the megaplasmid (not the main chromosome). I had to add them in manually to the external functions files. I saw that other tRNAs with similar SPO ids (like SPO_tRNA-Ser-5) had no Gene ID and only the tRNA and the amino acid code as their description (like tRNA-Ser), so I added these two as tRNA-Phe and tRNA-Thr, respectively.

And there are 53 genes without a Genbank/RefSeq locus tag. These include functions like transposases, transcriptional regulators, riboflavin synthase, gluconolactonase, and hypothetical proteins. 

#### Importing all external functions

We could concatenate and import all at once, but I did each annotation source separately so that I could easily find and fix any issues with number of fields (stray tabs, etc).

```
anvi-import-functions -c R_POM_DSS3-contigs.db -i SPO_external_functions.txt
anvi-import-functions -c R_POM_DSS3-contigs.db -i GENE_ID_external_functions.txt
anvi-import-functions -c R_POM_DSS3-contigs.db -i LOCUS_external_functions.txt
```

## Running de novo annotations

To get Pfams, HMM hits for single-copy core genes, COGs and KEGG (including BRITE) into the db.

```
anvi-run-hmms -c R_POM_DSS3-contigs.db
anvi-run-scg-taxonomy -c R_POM_DSS3-contigs.db
anvi-run-ncbi-cogs -c R_POM_DSS3-contigs.db -T 6
anvi-run-kegg-kofams -c R_POM_DSS3-contigs.db -T 6
anvi-run-pfams -c R_POM_DSS3-contigs.db -T 6
```

## Adding TnSeq annotations 

The Moran Lab has a library of TnSeq mutants. They want each gene in the database to be annotated with available mutants to facilitate sharing these mutants with collaborators. Essentially, each gene will be annotated with the Plate and Well number of each mutant in their library that matches to that gene.

I have a file from Lidimarie Trujillo Rodriguez which matches each mutant to the gene's SPO locus. Below, I use Python to match each mutant to the corresponding gene and formulate a functions txt file with the TnSeq annotations.

```python
import pandas as pd
tnseq_library = pd.read_csv("TnSeq_Data/Rpom_reannotation/Results/Re-annotated-MASTER-384-RA-UPDATED.csv")

spo_to_gcid = {}
with open("gcid_to_SPO_locus.txt", 'r') as f:
     for line in f.readlines()[1:]:
             gcid, spo = line.strip().split("\t")
             spo_to_gcid[spo] = gcid

annot_dict = {}
not_found = []
more_than_one = []

for idx in tnseq_library.index:
    t = tnseq_library.loc[idx, 'Locus']
    if t not in spo_to_gcid:
        not_found.append(t)
        continue
    gcid = spo_to_gcid[t]
    plate = tnseq_library.loc[idx, 'Master.Plate.Number']
    well = tnseq_library.loc[idx, 'Destination.Well']
    func_string = f"Plate {plate} Well {well}"
    if gcid in annot_dict.keys():
        more_than_one.append(gcid)
        annot_dict[gcid]['function'] += "; " + func_string
    else:
        annot_dict[gcid] = {'source': 'TN_MUTANT_AVAILABLE', 'accession': t, 'function': func_string, 'e_value': 0}

len(not_found)
# 35
set(not_found)
#{'Multiple genes for Tn mutant- Overlapping CDS', 'Intergenic'}

len(set(more_than_one))
# 1597

annot_df = pd.DataFrame(annot_dict)
annot_df = annot_df.T
annot_df.index.rename('gene_callers_id', inplace=True)
annot_df.to_csv("TnSeq_annotations.txt", sep="\t")
```

There are 2,871 TnSeq mutant annotations to add to the R pom database. Within these, there are 1,597 gene calls with more than one mutant annotated. I added the annotations to the database with the following commmand.

```
anvi-import-functions -c R_POM_DSS3-contigs.db -i TnSeq_annotations.txt
```


## Making profile databases for the Moran Lab's transcriptome samples
Mary Ann, Zac and Christa sent me an Excel spreadsheet of transcriptome samples and their SRA accession numbers. Here I will map these samples to the R. pomeroyi genome to produce profile databases, in a way that is compatible with the R. pomeroyi contigs database.

One note: all the following steps were run on the Meren Lab cluster for high-throughput processing.

### Downloading the samples

I used Matt Schechter's `fasterq-dump` script which is essentially a wrapper script that runs `fastq-dump` on each SRA accession number provided to the sample.

```
cut -f 1 combined_data_cbs.txt | grep SRR > SRR_Acc_List.txt
mkdir 00_LOGS
clusterize -o 00_LOGS/fasterq_dump.log "bash download_fasterq_dump.sh 'SRR_Acc_List.txt'"
# make sure there are no errors during download
grep -i "error" 00_LOGS/* | tr ":" "\t" | cut -f 1 | sed 's|00_LOGS/||' | tr "_" "\t" | cut -f 1 | sort -u
```

This downloads all the samples (in gzipped format) in the `02_FASTA` directory.

### Quality control

Normally I would use anvi'o's snakemake workflow for mapping all the samples to the genome in a high-throughput, automated, and easily reproducible fashion - it would handle everything from QC to profiling. But unfortunately, these are single-end samples, and our snakemake workflow only works with paired-end samples at the moment.

This means I need to run all the steps of the mapping workflow manually.  I'm basically going to write a bunch of loops to run each of the workflow steps on every sample.

The first step is quality control. Previously I have used [FASTX-toolkit](http://hannonlab.cshl.edu/fastx_toolkit/commandline.html#fastq_quality_filter_usage), with the parameters from the [Landa et al paper](https://www.nature.com/articles/ismej2017117) methods section:

>Quality control was performed on 249 million 50-bp reads (10±2 million reads per sample; Supplementary Table S1) using the FASTX toolkit, imposing a minimum quality score of 20 over 80% of read length.

So, `-p 20 -q 80`. I will do the same here.

#### Installing FASTX-toolkit

I will come right out and say that this toolkit is out-dated and troublesome to install. Last time, I didn't record what I did to install it, and I wish I had because I ran into the same problems again. So this time, I am recording everything here.

First step: install the `libgtextutils` package, which is a dependency of FASTX-toolkit.

```
cd /project2/meren/shared
wget https://github.com/agordon/libgtextutils/releases/download/0.7/libgtextutils-0.7.tar.gz
tar xvzf libgtextutils-0.7.tar.gz
cd libgtextutils-0.7/

# need to change some of the code, see explanation below
vi src/gtextutils/text_line_reader.cpp

# install
./configure --prefix /project2/meren/shared/libgtextutils-install
make
make install
```

As alluded to above, this installation initially produced the following error:

```
text_line_reader.cpp:47:9: error: cannot convert 'std::istream' {aka 'std::basic_istream<char>'} to 'bool' in return
   47 |  return input_stream ;
```

I found this error and the solution on [Stack Overflow](https://stackoverflow.com/questions/38659115/make-fails-with-error-cannot-convert-stdistream-aka-stdbasic-istreamchar). To fix it, I changed the `text_line_reader.cpp` file to have

```
return (bool)input_stream ;
```

instead of the original line `input_stream ;`.

Once I changed that, I re-ran the installation commands to install the package in the `/project2/meren/shared/libgtextutils-install` directory.

Step Two: Install FASTX-toolkit

```
cd /project2/meren/shared
wget https://github.com/agordon/fastx_toolkit/releases/download/0.0.14/fastx_toolkit-0.0.14.tar.bz2
tar -xjvf fastx_toolkit-0.0.14.tar.bz2
cd fastx_toolkit-0.0.14

# now you need to make sure the configure script can find the libgtextutils package that was just installed
# otherwise you get the "No package 'gtextutils' found" error
# solution explained here: http://hannonlab.cshl.edu/fastx_toolkit/pkg_config_email.txt
export PKG_CONFIG_PATH=/project2/meren/shared/libgtextutils-install/lib/pkgconfig/:$PKG_CONFIG_PATH

# another code fix, see below
vi src/fasta_formatter/fasta_formatter.cpp

# install
./configure --prefix /project2/meren/shared/fastx-toolkit-install
make
make install
```

Again, the installation produces an error:

```
fasta_formatter.cpp:105:9: error: this statement may fall through [-Werror=implicit-fallthrough=]
```
For this one, the solution is explained here: [agordon/fastx_toolkit#14](https://github.com/agordon/fastx_toolkit/issues/14)
and again requires changing the code in `fasta_formatter.cpp` to get rid of the `usage()` function and replace its only function call with the code that used to be in its function definition.

After that, run the install commands again to install FASTX-toolkit in the directory `/project2/meren/shared/fastx-toolkit-install`

Finally: test that it is working by getting the help output:

```
/project2/meren/shared/fastx-toolkit-install/bin/fastq_quality_filter -h
```

#### Quality control loop

Here is the code to run quality filtering on each sample.

```
for file in 02_FASTA/*.fastq.gz; do   \
  name=$(basename $file | sed 's/.fastq.gz//');   \
  echo $name;   \
  gunzip $file;  \
  clusterize -j $name -x "/project2/meren/shared/fastx-toolkit-install/bin/fastq_quality_filter -q 20 -p 80 -z -i 02_FASTA/${name}.fastq -o 01_QC/${name}-QC.fastq.gz" -o 00_LOGS/${name}_fastq_quality_filter.log; \
done
```

Note that this tool does not accept compressed input, so I had to `gunzip` each sample before passing it to the program. The `-z` flag means that it will produce gzipped output. Note also that here I am using the Meren Lab's `clusterize` wrapper script for submitting slurm jobs to our cluster - anyone without the `clusterize` script can simply run the command in quotation marks, ie `fastq_quality_filter -q 20 -p 80 -z -i 02_FASTA/${name}.fastq -o 01_QC/${name}-QC.fastq.gz`.

UPDATE: After the profiling steps described below, I gzipped these original fasta files again so they take up less space.

### Get a fasta file from the database

To do mapping, you need a fasta file to use as reference, and the contig headers in this file need to match those in the contigs database in order to not get errors during profiling. So I exported the contig sequences as a fasta file directly from the contigs database:

```
anvi-export-contigs -c 03_CONTIGS/R_POM_DSS3-contigs.db -o R_POM_DSS3.fa
```

### Mapping transcriptome samples to the R. pom genome

To prepare for mapping, I indexed the fasta file that I exported from the contigs database:

```
mkdir 04_MAPPING

# index
bowtie2-build R_POM_DSS3.fa 04_MAPPING/R_POM_DSS3
```

Here is the script that runs all the mapping steps for a single sample:

```
#!/bin/bash
# map a single sample of single-end reads and convert the output to BAM

sample_file=$1

name=$(basename $1 | sed 's/.fastq.gz//')
echo "working on sample $name"

bowtie2 --threads 4 -x 04_MAPPING/R_POM_DSS3 -U $sample_file --no-unal -S 04_MAPPING/${name}.sam 2>&1 > 00_LOGS/${name}-bowtie.log

samtools view -F 4 -bS 04_MAPPING/${name}.sam > 04_MAPPING/${name}-RAW.bam
rm 04_MAPPING/${name}.sam

anvi-init-bam 04_MAPPING/${name}-RAW.bam -o 04_MAPPING/${name}.bam
rm 04_MAPPING/${name}-RAW.bam
```

That script is called `map_one_sample.sh`, and here is a loop to run it for each sample:

```
# map
for samp in 01_QC/*; do \
  name=$(basename $samp | sed 's/.fastq.gz//')
  clusterize "./map_one_sample.sh $samp" -j $name -x -o 00_LOGS/map_${name}.log; \
done
```

### Profiling (converting BAM files to profile databases)

This is the part where we take the BAM files with read mapping data and convert them into database format.

```
mkdir 05_PROFILE

for samp in 01_QC/*; do \
  name=$(basename $samp | sed 's/.fastq.gz//'); \
  echo $name; \
  sampname=$(echo $name | sed 's/-QC//g'); \
  clusterize "anvi-profile -c 03_CONTIGS/R_POM_DSS3-contigs.db -i 04_MAPPING/${name}.bam -o 05_PROFILE/${name} -S $sampname -T 4" -j anvi_profile -x -o 00_LOGS/${sampname}_anvi_profile.log
done
```

### Merging profiles

Combining all samples into one database :)

```
clusterize -j anvi-merge -x -o 00_LOGS/anvi_merge.log "anvi-merge -c 03_CONTIGS/R_POM_DSS3-contigs.db -o 06_MERGED 05_PROFILE/*/PROFILE.db"
```

## Miscellaneous data layers 

Sam Miller generated a GENES database from the merged profile by visualizing it in `anvi-interactive` using the `--gene-mode` flag. He then imported miscellaneous data into that GENES database using the anvi'o program `anvi-import-misc-data`. Here is what he had to say about that process:

> I added metadata to the genes database derived from the profile database on the types of experiments from which transcriptomes were sampled. I added this data as a layer order. The experiments seemed to have four major themes that I used to determine sample order: monocultures; co-cultures with phytoplankton; monocultures and co-cultures with up to three bacteria including R. pom; and "seawater invasion" experiments in which R. pom was introduced to seawater. Monocultures, for instance, consist of different experiments with R. pom grown under normal culture conditions and on a phytoplankton filtrate.

> A table relating transcriptome sample (SRA) names to experimental conditions can be found at `/project2/meren/PROJECTS/CCoMP_R_POM/combined_data_cbs_condition_order.txt`. The external data file added to the genes database is `layer_condition_order.txt`.

## Version 00 of the profile database 
At this point in the process we realized that we needed to implement some sort of versioning system for the databases, as we had more samples to map to the _R. pom_ genome. The aforementioned profile database, along with its associated genes database, was labeled as version 00. 

On the Meren Lab cluster, it can be found at the following path: `/project2/meren/PROJECTS/CCoMP_R_POM/06_MERGED/VER_00`. Contact the Meren Lab for access to this version of the data.

## Mapping more transcriptome samples (Version 01)

Zac sent me 11 of his samples and they need to be added to the profile database. This is what he said:

> I’m using BBtools for my own work with the data, and I have already qc’ed it with BBduk to trim adapters. I’ll put all 11 of these files into a zip file and share a drive link with you when it’s ready. Let me know if for some reason you would prefer the non-qc’ed data. Everything I have is also single-ended 50 bp reads. BBmap did a great job with this, but I know you use a different tool in your workflow

So, these are already QC'ed, albeit with a different tool than what I used for the others. I uploaded the files to the Meren Lab cluster and slightly renamed them so that the file extensions are consistent with the other samples we already had (ie, they all now end in `-QC.fastq.gz`).

Now they need to be mapped and merged with the other samples to make a new profile database.

The other samples have been mapped with bowtie, so for consistency's sake I will use bowtie here. I did run a short test in which I mapped a few samples to the _R. pom._ genome with both BBmap and bowtie, and there were no appreciable differences in their mapping patterns.

Here is the loop to run mapping:

```
# map
for samp in 01_QC/DSS3*; do \
  name=$(basename $samp | sed 's/.fastq.gz//')
  clusterize "./map_one_sample.sh $samp" -j $name -x -o 00_LOGS/map_${name}.log; \
done
```

Once this finished (no errors), I ran this loop to make single-sample profile databases:

```
for samp in 01_QC/DSS3*; do \
  name=$(basename $samp | sed 's/.fastq.gz//'); \
  echo $name; \
  sampname=$(echo $name | sed 's/-QC//g'); \
  clusterize "anvi-profile -c 03_CONTIGS/R_POM_DSS3-contigs.db -i 04_MAPPING/${name}.bam -o 05_PROFILE/${name} -S $sampname -T 4" -j anvi_profile -x -o 00_LOGS/${sampname}_anvi_profile.log
done
```

Finally, I had to make a new merged profile. It was versioned as `VER_01`. 

```
clusterize -j anvi-merge -x -o 00_LOGS/anvi_merge_ver_01.log "anvi-merge -c 03_CONTIGS/R_POM_DSS3-contigs.db -o 06_MERGED/VER_01 05_PROFILE/*/PROFILE.db"
```

# Uploading the databases to Zenodo

The contigs database and the merged profile `VER_01` (which includes all transcriptome samples provided by the Moran Lab) are the first version of this dataset that is being uploaded to Zenodo for community-wide sharing. See the [C-CoMP Zenodo Community page](https://zenodo.org/communities/c-comp/).

# Adding proteomics data

Zac provided a table of normalized spectral abundance counts from proteomic samples that are matched to several of the transcriptomes already in the contigs database. We want to visualize these abundances as barplots in the interactive interface. Since this is gene-specific data, we will use gene mode to do it.

Gene mode requires us to specify a collection and a bin (to work with genes within a specific set of splits), so I first added a default collection containing everything in the database:

```
anvi-script-add-default-collection -c R_POM_DSS3-contigs.db -p 06_MERGED/VER_01/PROFILE-VER_01.db
```

Then, I ran the interactive interface in gene mode to obtain a GENES database:

```
anvi-interactive -c R_POM_DSS3-contigs.db -p 06_MERGED/VER_01/PROFILE-VER_01.db -C DEFAULT -b EVERYTHING --gene-mode
```

While this was running, it gave this error indicating that something is up with two gene calls `1598` and `1599`:

```
🚑 SOMETHING WEIRD HAPPENED 🚑
===============================================
Please read this carefully as something sad just happened. While anvi'o was
trying to recover gene level coverage stats, it became clear that a few gene
calls were not found in splits they were meant to be found. These genes will not
be a part of any downstream reporting :/ It is extremely difficult to even
entertain the idea why this might have happened, we suspect it is some sort of
artifact left behind from the use of external gene calls. Regardless, here are
the gene calls that cause you this headache, in case you would like to go after
this and find out why is this happening: 1598, 1599
```

And correspondingly, two gene calls (probably the same ones) don't have any coverage information:

```
WARNING
===============================================
Anvi'o observed something weird while it was processing gene level coverage
statistics. Some of the gene calls stored in your contigs database (2 of them,
precisely) did not have any information in gene level coverage stats dictionary.
It is likely they were added to the contigs database *after* these gene level
coverage stats were computed and stored. One way to address this is to remove
the database file for gene coverage stats and re-run this step. [.....]
```
So, these two gene calls were added later to the database (I did not record this process, shame on me) and as a result have no coverage information from the samples we mapped to the genome. This is not a critical error, it just means that the two affected gene calls will not be shown in the 'genes mode' interface.


The table Zac sent me has the normalized protein abundances already matched to the SPO IDs of each gene. I just have to convert the index to match the gene caller IDs from the contigs database. I quickly extracted a table of the SPO ID annotations for each gene:

```
anvi-export-functions  -c R_POM_DSS3-contigs.db --annotation-sources SPO_ID -o gene_call_to_SPO.txt
```

and then I converted the abundance table to one indexed by gene call. Since the sample (column) names in these files match to their corresponding transcriptome sample names, I altered those to indicate that they are proteome samples instead.

```python
import pandas as pd
prot = pd.read_csv("20221121_DSS3_norm_prot_abund.txt", sep="\t", index_col=0)

# change index to gene caller id
spo = pd.read_csv("gene_call_to_SPO.txt", sep="\t")
spo_to_gcid = dict(zip(spo.accession, spo.gene_callers_id))
prot.rename(index=spo_to_gcid, inplace=True)

# rename samples in columns
new_cols = [x + "_proteome" for x in prot.columns]
prot.columns = new_cols

prot.to_csv("gcid_norm_prot_abund.txt", sep="\t")
```

Finally, I imported the protein abundance data into the genes database, and visualized it again to see that data in the interface:

```
anvi-import-misc-data -p 06_MERGED/VER_01/GENES/DEFAULT-EVERYTHING.db -t items gcid_norm_prot_abund.txt

# visualize with gene mode
anvi-interactive -c R_POM_DSS3-contigs.db -p 06_MERGED/VER_01/PROFILE-VER_01.db -C DEFAULT -b EVERYTHING --gene-mode
```

I then adjusted the settings to show only the proteomic data layers with their paired transcriptomes, and saved the state. To load the interactive interface with these settings, run the following:

```
anvi-interactive -c R_POM_DSS3-contigs.db -p 06_MERGED/VER_01/PROFILE-VER_01.db -C DEFAULT -b EVERYTHING --gene-mode --state-autoload proteomes
```

# Uploading the genes database to Zenodo (and important access information)

I made a new version of the shared _R. pomeroyi_ databases on Zenodo at [https://zenodo.org/deposit/7884613](https://zenodo.org/deposit/7884613) which includes the newly-created genes database (`DEFAULT-EVERYTHING.db`) with the proteomic data layers. The genes database is associated with profile db VER_01 (the one with all the transcriptome samples).

Unfortunately, it seems that Zenodo does not allow uploading of folders, so I could not keep the database in its expected directory structure of `GENES/DEFAULT-EVERYTHING.db`. Thus, to access this database in the anvi'o interactive interface, you need to run the following commands after the Zenodo datapack is downloaded:

```
mkdir GENES
mv DEFAULT-EVERYTHING.db GENES/
```

Only then will you be able to run the visualization command for the proteomics data:

```
anvi-interactive -c R_POM_DSS3-contigs.db -p 06_MERGED/VER_01/PROFILE-VER_01.db -C DEFAULT -b EVERYTHING --gene-mode --state-autoload proteomes
```

If you don't move the genes database to its expected location, then anvi'o will automatically try to generate a new database for you (which will not include the proteomic data in it, since we imported that into the database after we created it).
