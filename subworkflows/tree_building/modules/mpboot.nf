// MPBOOT module

nextflow.enable.dsl = 2

process MPBOOT {

    tag "${sample}: ${job_id}"

    input:
    tuple val(job_id),
        val(sample), 
        val(rep),
        path(afm)

    output:
    tuple val(job_id),
        val(sample), 
        path("final_tree.newick"), emit: tree
    
    script:
    """
    python ${baseDir}/bin/build_tree/create_fasta.py  ${afm}
    mpboot -s genotypes.fa -bb 1000 -con -bi 100 -t 0.5
    mv genotypes.fa.contree final_tree.newick
    """

    stub:
    """
    touch final_tree.newick
    """

}
