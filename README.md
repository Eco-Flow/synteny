# MCScanX-nf 
Nextflow pipeline to run MCScanX

# General information

Using the docker container found here: 
Using this version: 

To run on different platforms, you need to create a profile. Docker (local), and myriad (ucl cliuster) are already avaialble, please ask and I can create a profile for your env. Use the flag `-profile` to request this (essential) 

# How to run test

Will attempt to run the basic example script:

`cactus jobStore examples/evolverMammals.txt examples/evolverMammals.hal --root mr`

To do this we use Nextflow (locally wiht docker installed):

`nextflow run main.nf -profile docker --trial 1`

#`trial 1` activates the trail version of the nextflow pipeline
#Notice one - for Nextflow options, and two -- for pipeline options.

# How to run with local wasp genomes.

INPUT: Two whole genomes (e.g. Polistes dominula and Vespa crabro)

nextflow run chriswyatt1/mcsanx-nf --genomes "Polistes_canadensis,Vespa_crabro" -with-singularity ?
