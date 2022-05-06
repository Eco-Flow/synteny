# JCVI-nextflow 

Nextflow pipeline to run JCVI

Please cite "Tang et al. (2008) Synteny and Collinearity in Plant Genomes. Science" if you use this Nextflow wrapper of JCVI MCSxanX software.

# General information

This is a developmental workflow running JCVI, to look at gene synteny. 

All you need is either a genome in fasta format with an annotation file in gff3 (or gff augustus). OR you can supply a NCBI genome reference ID (which will be automatically downloaded).

To run on different platforms, you may need to create a profile. We recommend using the prebuilt Docker profile (tp run locally or through Gitpod), though if you are running on a HPC, you will need to change this. Please open an issue and I can help create a profile for your environment. Use the flag `-profile` to choose the environment in the script command. These are found in the folder `conf`

# How to run locally

Prerequistites : Docker account and Docker installed on your machine. Nextflow installed (https://www.nextflow.io/; v22 and above [DSL2]), plus java at at least 1.8.

To run Nextflow (locally with docker installed), use the following command:

`nextflow run main.nf -profile docker -bg -resume --input data/Example.csv`

#Notice one - for Nextflow options, and two -- for pipeline options.


This is what the input template looks like (Example.csv):

`D_melanogaster,GCF_000001215.4
A_mellifera,GCF_003254395.2`

You can also run your own genomes through this program (or mixed with NCBI ones), using the following format:\

`B_impatiens,/Users/cwyatt/Desktop/B_impatiens_Genome**.fasta**,/Users/cwyatt/Desktop/B_impatiens**.gff**
A_mellifera,GCF_003254395.2`

Where NCBI input has two comma separated columns and your own data has three coloumns (Name, Genome.fasta and GFF file).

# Run with Gitpod

Prerequistites : A browser (Ideally, Chrome or Firefox \[tested\]) and a Github account.

The simplest way to run the pipeline is to use Gitpod. This is a free (up to 50 hours a month) cloud environment, which has been loaded with all the tools you need.

Simply click this link: https://gitpod.io/#https://github.com/chriswyatt1/jcvi-nextflow

Then login in to Github, which will open up an environment to run the code, using the same command listed above (nextflow...).

`nextflow run main.nf -profile docker -bg -resume --input data/Example.csv`

