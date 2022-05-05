# JCVI-nextflow 

Nextflow pipeline to run JCVI

# General information

This is a developmental workflow running JCVI, to look at gene synteny. 

All you need is a genome in fasta and a annotation file in gff3 (or gff augustus).

To run on different platforms, you need to create a profile. We recommend using the Docker (local) one, though if you are running on a HPC, you will need to change this. Please open an issue and I can help create a profile for your environment. Use the flag `-profile` to choose the environment in the script command. These are found in the folder `conf`

# How to run locally

Prerequistites : Docker account and Docker installed on your machine. Nextflow installed (https://www.nextflow.io/), plus java at at least 1.8.

To run Nextflow (locally with docker installed), use the following command:

`nextflow run main.nf -profile docker -bg -resume --input data/Example.csv`

#Notice one - for Nextflow options, and two -- for pipeline options.


# Run with Gitpod

Prerequistites : A browser (Ideally, Chrome or Firefox \[tested\]) and a Github account.

The simplest way to run the pipeline is to use Gitpod. This is a free (up to 50 hours a month) cloud environment, which has been loaded with all the tools you need.

Simply click this link: https://gitpod.io/#https://github.com/chriswyatt1/jcvi-nextflow

Then login in to Github, which will open up an environment to run the code, using the same command listed above (nextflow...) .


# How to run with local wasp genomes.

INPUT: Two whole genomes (e.g. Polistes dominula and Vespa crabro)

nextflow run chriswyatt1/jcvi-nextflow -profile docker -bg -resume --input "Polistes_canadensis,Vespa_crabro"
