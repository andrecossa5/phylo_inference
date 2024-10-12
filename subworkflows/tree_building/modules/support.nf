// SUPPORT module
nextflow.enable.dsl = 2

process SUPPORT {

    tag "${sample}: ${job_id}"

    input:
    tuple val(job_id),
        val(sample), 
        val(trees)

    output:
    tuple val(job_id),
        val(sample), 
        path("final_tree.newick"), emit: tree

    script:
    """
    python ${baseDir}/bin/build_tree/support.py \
    --trees "${trees}" \
    --support_method ${params.support_method} \
    --n_cores ${task.cpus}
    """

    stub:
    """
    touch final_tree.newick
    """

}
