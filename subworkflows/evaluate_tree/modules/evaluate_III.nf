// EVALUATE_III module
nextflow.enable.dsl = 2

process EVALUATE_III {

    tag "${sample}_${filtering}"

    publishDir "${params.outdir}/${sample}/${filtering}/${solver}/${metric}/", mode: 'copy'

    input:
    tuple val(sample), 
        val(filtering), 
        val(solver),
        val(metric),
        path(obs_tree),
        path(original_input)
    
    output:
    tuple val(sample), 
        val(filtering), 
        val(solver),
        val(metric),
        path('muts.csv'),
        path('clones.csv'), emit: results
    
    script:
    """
    Rscript ${baseDir}/bin/evaluation_III.r \
    ${obs_tree} \
    ${original_input} \
    ${task.cpus}
    """

    stub:
    """
    touch muts.csv
    touch clones.csv
    """

}