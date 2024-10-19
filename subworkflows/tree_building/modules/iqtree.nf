// IQTREE module

nextflow.enable.dsl = 2

process IQTREE {

    tag "${sample}: ${job_id}"

    input:
    tuple val(job_id),
        val(sample), 
        val(bin_key),
        val(tree_key),
        val(rep),
        path(afm)

    output:
    tuple val(job_id),
        val(sample), 
        path("final_tree.newick"), emit: tree
    
    script:
    """
    python ${baseDir}/bin/build_tree/create_fasta.py  ${afm}
    iqtree -s genotypes.fa -B 1000
    python ${baseDir}/build_tree/iqtree_to_newick.py
    """

    stub:
    """
    touch final_tree.newick
    """

}