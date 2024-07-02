// TEST_LINEAGE_ASSOCIATION module
nextflow.enable.dsl = 2 

//

process TEST_LINEAGE_ASSOCIATION {

    tag "${sample}: ${filtering_key}"

    // Publish
    publishDir "${params.outdir}/${sample}/${filtering_key}", mode: 'copy'

    input:
    tuple val(sample), 
        val(filtering_key), 
        path(input_folder),
        path(nodes),
        path(edges)
 
    output:
    tuple val(sample), 
        val(filtering_key), 
        path('lineage_association.csv'), emit: lineage_association_stats
    
    script:
    """
    Rscript ${baseDir}/bin/process_tree/test_lineage_association.r \
    --nodes ${nodes} \
    --edges ${edges} \
    --meta ${input_folder}/meta.csv \
    --lineage_column ${params.lineage_column}
    """

    stub:
    """
    touch lineage_association.csv
    """

}
