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
        path('morans_muts.csv'),
        path('couplings.csv'), emit: results
    
    script:
    """
    python ${baseDir}/bin/evaluation_III.py \
        --sample_name ${sample} \
        --filtering ${filtering} \
        --solver ${solver} \
        --metric ${metric} \
        --obs_tree ${obs_tree} \
        --input_folder ${original_input}
    """

    stub:
    """
    touch morans_muts.csv
    touch couplings.csv
    """

}