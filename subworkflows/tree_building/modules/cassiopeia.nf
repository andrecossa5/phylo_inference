// CASSIOPEIA module

nextflow.enable.dsl = 2


process CASSIOPEIA {

    tag "${sample}: ${job_id}, rep=${rep}"

    input:
    tuple val(job_id),
        val(sample), 
        val(rep),
        path(afm)

    output:
    tuple val(job_id),
        val(sample), 
        path("*.newick"), emit: tree
     
    script:
    """
    python ${baseDir}/bin/build_tree/build_cassiopeia.py \
    --afm ${afm} \
    --path_pickles ${params.path_pickles} \
    --sample ${sample} \
    --job_id ${job_id} \
    --solver ${params.cassiopeia_solver} \
    --metric ${params.distance_metric} \
    --boot_replicate ${rep} \
    --n_cores ${task.cpus}
    """

    stub:
    """
    touch rep_${rep}.newick
    """

}


//