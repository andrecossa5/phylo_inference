// SUPPORT module
nextflow.enable.dsl = 2

process SUPPORT {

    tag "${sample}: ${filtering_key}"

    // Publish
    publishDir "${params.outdir}/${sample}/${filtering_key}", mode: 'copy'

    input:
    tuple val(sample), 
        val(filtering_key), 
        val(trees)
    
    output:
    tuple val(sample), 
        val(filtering_key), 
        path("final_tree.newick"), emit: tree

    script:
    """
    python ${baseDir}/bin/cassiopeia/support.py \
    --filtering_key ${filtering_key} \
    --solver ${params.cassiopeia_solver} \
    --trees "${trees}" \
    --n_cores ${task.cpus}
    """

    stub:
    """
    touch final_tree.newick
    """

}
