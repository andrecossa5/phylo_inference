// PREP_INPUT module

nextflow.enable.dsl = 2 

//

process PREP_INPUT {

    tag "${sample}: ${filtering_key}"

    // Publish
    publishDir "${params.outdir}/${sample}/${filtering_key}/", mode: 'copy'

    input:
    tuple val(sample), val(filtering_key)

    output:
    tuple val(sample), val(filtering_key), path(input_folder), emit: input_folder
    
    script:
    """
    python ${baseDir}/bin/prep_input/prep_input.py \
    --path_data ${params.path_data} \
    --nmads ${params.nmads} \
    --path_meta ${params.path_meta} \
    --path_priors ${params.path_priors} \
    --path_filtering ${params.path_filtering} \
    --filtering_key ${filtering_key} \
    --sample ${sample} \
    --n_cores ${task.cpus} \
    --lineage_column ${params.lineage_column} \
    --cell_file ${params.cell_file}
    """

    stub:
    """
    mkdir input_folder
    """

}
