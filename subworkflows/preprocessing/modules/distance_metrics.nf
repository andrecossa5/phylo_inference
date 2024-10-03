// DISTANCES module

nextflow.enable.dsl = 2 

//

process DISTANCE_METRICS {

    tag "${sample}: ${job_id}"

    input:
    tuple val(job_id),
        val(sample), 
        val(bin_key), 
        val(tree_key), 
        path(input_folder), 
        val(rep),
        val(dists)

    output:
    tuple val(job_id),
        val(sample),  
        path(input_folder), 
        path("distance_metrics.csv"), emit: distance_metrics
    
    script:
    """
    python ${baseDir}/bin/pp/distance_metrics.py \
    --path_dists "${dists}" \
    --replicates "${rep}" \
    --path_meta ${input_folder}/cell_meta.csv \
    --job_id ${job_id} \
    --K ${params.K} \
    --lineage_column ${params.lineage_column}
    """

    stub:
    """
    touch distance_metrics.csv
    """

}

