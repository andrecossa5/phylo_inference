// publish_images module

nextflow.enable.dsl = 2

process publish_images {

    tag "${sample}, filtering_key: ${filtering_key}"
    publishDir "${params.outdir}/${sample}", mode: 'copy'

    input:
    tuple val(sample), val(filtering_key), path(dists), path(tree_raw), path(tree_muts)
    
    output:
    tuple val(sample),
        val(filtering_key), 
        path("${sample}_${filtering_key}_distances.png"),
        path("${sample}_${filtering_key}_tree_raw.png"),
        path("${sample}_${filtering_key}_tree_muts.png"), emit: plots
    
    script:
    """
    echo "Publishing plots..."
    """

    stub:
    """
    echo "Publishing plots..."
    """

}
