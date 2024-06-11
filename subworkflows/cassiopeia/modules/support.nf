// SUPPORT module
nextflow.enable.dsl = 2

process SUPPORT {

    tag "${sample}: ${filtering_key}, ${boot_method}, ${solver}"

    // Publish
    publishDir "${params.outdir}/${sample}/${filtering_key}/${solver}/", mode: 'copy'

    input:
    tuple val(sample), 
        val(filtering_key), 
        val(boot_input),
        val(boot_method),
        val(boot_replicates),
        val(solver),
        val(trees)
    
    output:
    tuple val(sample), 
        val(filtering_key), 
        val(solver),
        path("${filtering_key}_${solver}_tree.newick"), emit: tree

    script:
    """
    python ${baseDir}/bin/cassiopeia/support.py \
    --filtering_key ${filtering_key} \
    --solver ${solver} \
    --trees "${trees}" \
    --n_cores ${task.cpus}
    """

    stub:
    """
    touch ${filtering_key}_${solver}_tree.newick
    """

}
