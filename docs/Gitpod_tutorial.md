# Running the synteny pipeline (Quick start)

To run the pipeline quickly (on small datasets, <3 species), you can take advantage of a free online code development platform (Gitpod).

## Prerequistites

A chrome based web browser
A github account (and your credentials).

## Instructions

1. Click on this link: https://gitpod.io/#https://github.com/Eco-Flow/synteny/tree/main
2. Enter your github credentials to access Gitpod (you get 500 free credits a month of the free tier).
3. In the terminal (bottom right) enter the basic test run command:

`nextflow run main.nf -profile docker,test`

4. If this runs sucessfully, you can now try with your own input data.
5. Go to NCBI, and find the Refseq IDs of species you want to compare. 
WARNING: You will require your species to have a chromosome level genome with N50 > 1000000.
6. Make an input comma separated value file to look like the following:
```
Drosophila_yakuba,GCF_016746365.2
Drosophila_simulans,GCF_016746395.2
Drosophila_santomea,GCF_016746245.2
```
7. Save this file to a name you will remember, e.g. Drosophila_Refseq_List.csv
8. Now try to run Nextflow again, but point to this new file you have created:
`nextflow run main.nf -profile docker --input /path/to/Drosophila_Refseq_List.csv`;
(this should take around 15-25 minutes depending on genome sizes).

On Gitpod my "path/to" is `/workspace/synteny/Drosophila_Refseq_List.csv`

NOTICE: We remove the test profile, and just use docker (our container engine, to get the software).

9. Once completed, explore the `Results` folder to see some of the output. On the main README on this repository you will find full explanation of the output files.

