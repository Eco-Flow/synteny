process DOWNLOAD {
    label 'download'
    container = 'chriswyatt/perl_r_e1071:latest'           
          
    output:
        path("evolverMammals.txt") , emit: test_data

    script:
    """
    wget https://raw.githubusercontent.com/ComparativeGenomicsToolkit/cactus/master/examples/evolverMammals.txt

    """
}