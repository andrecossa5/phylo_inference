// CASSIOPEIA module

nextflow.enable.dsl = 2

//

process CASSIOPEIA {

    tag "${sample}_${filtering}_${metric}_${solver}_${boot_strategy}_${params.n_boot}"
    publishDir "${params.outdir}/", mode: 'copy'

    input:
    tuple val(sample), 
        val(filtering),
        val(metric),
        val(solver),
        val(boot_strategy),
        val(n_boot)

    output:
    path "${sample}_${filtering}_${metric}_${solver}_${boot_strategy}_${params.n_boot}", emit: results
    
    script:
    """
    python ${baseDir}/bin/phylo_cassiopeia.py \
        --p ${params.path_data} \
        --sample ${sample} \
        --filtering ${filtering} \
        --metric ${metric} \
        --solver ${solver} \
        --boot_strategy ${boot_strategy} \
        --n_boot ${params.n_boot} \
        --ncores ${task.cpus} 
    """

    stub:
    """
    touch "${sample}_${filtering}_${metric}_${solver}_${boot_strategy}_${params.n_boot}"
    """

}
