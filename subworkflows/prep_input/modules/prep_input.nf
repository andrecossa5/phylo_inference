// PREP_INPUT module

nextflow.enable.dsl = 2 

//

process PREP_INPUT {

    tag "${sample}_${filtering}"

    // Publish
    publishDir "${params.outdir}/${sample}/${filtering}/", mode: 'copy'

    input:
    tuple val(sample), val(filtering)

    output:
    tuple val(sample), val(filtering), path(input_folder), emit: input_folder
    
    script:
    """
    python ${baseDir}/bin/prep_input/prep_input.py \
    -p ${params.path_data} \
    --sample ${sample} \
    --ncores ${task.cpus} \
    --filtering ${filtering} \
    --GT_reference ${params.GT_reference}
    """

    stub:
    """
    mkdir input_folder
    """

}
