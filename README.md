# nf-synteny 

Nextflow pipeline to run Synteny analysis using JCVI

Please cite : "Tang et al. (2008) Synteny and Collinearity in Plant Genomes. Science",
if you use this Nextflow wrapper of JCVI MCSxanX software.

# General information

This is a developmental Nextflow workflow running JCVI, to look at gene synteny with chromosome level assemblies (less than 40 chromosomes or scaffolds ideally- for better visualisation). We will paint chromosomes from one species to another, translate scaffolds between species, and produce dot plots and basic statistics.

All you need is either a genome in fasta format with an annotation file in gff3 (or gff augustus). 
**OR** you can supply a Refseq genome reference ID (which will be automatically downloaded; e.g. GCF_000001215.4).
Do not submit a Genbank reference ID, as these will not work.

There is one main branch.
'main': which can run 2 or more samples against eachother pairwise, producing dotplots and chromosome plots, along with species wise statistics and gene statistics.

To run on different platforms, you may need to create a profile. We recommend using the prebuilt Docker and local profiles (to run locally or through Gitpod), though if you are running on a HPC, you will need to change this. Please open an issue and I can help create a profile for your environment. Use the flag `-profile` to choose the environment in the script command. These are found in the folder `conf`.

# Run with Gitpod (recommended)

Prerequistites : 
- A browser (Ideally, Chrome or Firefox \[tested\]).
- Github account.

Optional: Add a PDF viewer extension in Gitpod. Go to Extensions on left hand side, and install `vscode.pdf`.

The simplest way to run the pipeline is to use Gitpod. This is a free (up to 50 hours a month) cloud environment, which has been loaded with all the tools you need.

Simply click this link: https://gitpod.io/#https://github.com/chriswyatt1/jcvi-nextflow

Then login in to Github, which will open up an environment to run the code, using the same command listed above (nextflow...).

The example run is below (using three public Drosophila genomes):

`nextflow run main.nf -profile docker,local -resume --input data/Example.csv`

# Run locally

Prerequistites : 
- Docker. Make sure it is active log in on your machine.
- Java at least 1.8.
- Nextflow installed (https://www.nextflow.io/; v22 and above [DSL2].)
- Git.

To clone the repo: `git clone https://github.com/Eco-Flow/synteny.git`

Then `cd` into the repository on your machine.

To run Nextflow (locally with docker installed), use the following command:

`nextflow run main.nf -profile docker,local -resume --input data/Example.csv`


or with (if you download these three datasets manually- e.g. http://ftp.ensembl.org/pub/rapid-release/species/Vespula_germanica/GCA_905340365.1/genome/)

`--input example.csv`

#Notice, we use one `-` for Nextflow options, and two `--` for pipeline options.

# Changing the input 

Our example input template looks like this (data/Example.csv):

```
Drosophila_yakuba,GCF_016746365.2
Drosophila_simulans,GCF_016746395.2
Drosophila_santomea,GCF_016746245.2
```

You can also run your own genomes through this program (or mixed with NCBI ones), using the following format:

```
B_impatiens,/path/to/Desktop/B_impatiens_Genome**.f**a**,/path/to/Desktop/B_impatiens**.gff**
A_mellifera,GCF_003254395.2
```

Where NCBI input has two comma separated columns and your own data has three coloumns (Name, Genome.fasta and GFF file). To upload data simply drop an drag your files into the explorer on the left hand side. Or use public data as previously specified (or mix and match them). 

#To run with Gene Ontology information:

You need to provide the transcript Gene Ontology annotations from [chriswyatt1/Goatee](https://github.com/chriswyatt1/Goatee). These should be in the results/Go folder and are the ones labelled *transcript*.

Copy the *transcript* into a new folder, and then point to them with the flag `--go`.e.g. :
 
`nextflow run main.nf -profile docker,local -resume --input data/Example.csv --go /path/to/go/transcripts/folder`

# Run with Gitpod (for development of the pipeline). *For admins*

Prerequistites : 
- A browser (Ideally, Chrome or Firefox \[tested\]).
- Github account.

Optional: Add a PDF viewer extension in Gitpod. Go to Extensions on left hand side, and install `vscode.pdf`. 

The simplest way to run the pipeline is to use Gitpod. This is a free (up to 50 hours a month) cloud environment, which has been loaded with all the tools you need.

Simply click this link: https://gitpod.io/#https://github.com/chriswyatt1/jcvi-nextflow

Then login in to Github, which will open up an environment to run the code, using the same command listed above (nextflow...).

To upload data simply drop an drag your files into the explorer on the left hand side. Or use public data as previously specified. The example run is below:

`nextflow run main.nf -profile docker,local -resume --input data/Example.csv`

# Results

Once completed, you should have a folder called Results/Jcvi_results, in which there should be a:

1. Dot plot (<Species1><Species2>.pdf). Showing the chromosome synteny as a dot plot.
2. Chromosome synteny plot (<Species1><Species2>.macro.pdf). Showing a 1 to 1 chromosome mapping with lines drawn between syntenic chromosomes.
3. Depth plot (<Species1><Species2>.depth.pdf). Number of 1to1 or 1toMany orthologs detected.
4. Chromosome painted mappings (). Showing on graphic chromosomes, which sections are syntenic between two species.
5. Anchors. (<Species1><Species2>.anchors). This shows the syntenic blocks, gene by gene.

# Testing scripts in Docker 

To try out the individual scripts used in this workflow, check out the various containers used in modules/**.nf.

You can enter containers using the docker run command i.e.:

`docker run -it --rm ecoflowucl/jcvi:python-3.10_last-1522 bash`

e.g. in the above container you should have jcvi, so you can execute the following line:

`python -m jcvi.formats.gff ACTION`

use exit to leave the container:

`exit`

You can also enter docker with the same filesystem using:

`docker run -it --rm -v "$PWD":"$PWD" ecoflowucl/jcvi:python-3.10_last-1522 bash`
