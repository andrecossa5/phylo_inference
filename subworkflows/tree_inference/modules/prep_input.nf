// PREP_INPUT module

nextflow.enable.dsl = 2 

//

process PREP_INPUT {

    tag "${sample}_${filtering}"

    input:
    tuple val(sample), val(filtering)

    output:
    tuple val(sample), val(filtering), path(input_folder), emit: input_folder
    
    script:
    """
    python ${baseDir}/bin/prep_input.py \
    -p ${params.path_data} \
    --sample ${sample} \
    --filtering ${filtering}
    """

    stub:
    """
    mkdir input_folder
    """

}
