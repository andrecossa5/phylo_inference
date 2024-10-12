// CASSIOPEIA module

nextflow.enable.dsl = 2


process CASSIOPEIA {

    tag "${sample}: ${job_id}, rep=${rep}"

    input:
    tuple val(job_id),
        val(sample), 
        val(bin_key),
        val(tree_key),
        val(rep),
        val(afm)

    output:
    tuple val(job_id),
        val(sample), 
        path("${rep}.newick"), emit: tree
     
    script:
    """
    python ${baseDir}/bin/build_tree/build_cassiopeia.py \
    --afm ${afm} \
    --path_bin ${params.path_bin} \
    --path_tree ${params.path_distance_tree} \
    --bin_key ${bin_key} \
    --tree_key ${tree_key} \
    --boot_replicate ${rep} \
    --n_cores ${task.cpus}
    """

    stub:
    """
    touch ${rep}.newick
    """

}
