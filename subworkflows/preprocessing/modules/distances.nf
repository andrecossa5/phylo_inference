// DISTANCES module

nextflow.enable.dsl = 2 

//

process DISTANCES {

    tag "${sample}: ${job_id}, rep ${rep}"

    input:
    tuple val(job_id),
        val(sample), 
        val(bin_key), 
        val(tree_key), 
        path(input_folder), 
        val(rep)

    output:
    tuple val(job_id),
        val(sample), 
        val(bin_key), 
        val(tree_key), 
        val(rep), 
        path('dist.npz'), emit: distances
    
    script:
    """
    python ${baseDir}/bin/pp/distances.py \
    --AD ${input_folder}/AD.npz \
    --DP ${input_folder}/DP.npz \
    --path_bin ${params.path_bin} \
    --path_tree ${params.path_distance_tree} \
    --bin_key ${bin_key} \
    --boot_replicate ${rep} \
    --tree_key ${tree_key} \
    --n_cores ${task.cpus} \
    --boot_strategy ${params.boot_strategy} \
    --frac_char_resampling ${params.frac_char_resampling}
    """

    stub:
    """
    touch dist.npz
    """

}


