 // TREE_METRICS module

nextflow.enable.dsl = 2


process TREE_METRICS {

    tag "${sample}: ${job_id}"
    publishDir "${params.outdir}/${sample}/${job_id}", mode: 'copy'

    input:
    tuple val(job_id),
        val(sample), 
        path(tree)

    output:
    tuple val(job_id),
        val(sample), 
        path("tree_metrics.csv"),
        path("evo_coupling.csv"), emit: metrics
    
    script:
    """
    python ${baseDir}/bin/process_tree/tree_metrics.py \
    --tree ${tree} \
    --lineage_column ${params.lineage_column} \
    --job_id ${job_id}
    """

    stub:
    """
    touch tree_metrics.csv
    touch evo_coupling.csv
    """

}
