// PROCESS_TREE module

nextflow.enable.dsl = 2


process PROCESS_TREE {

    tag "${sample}: ${job_id}"
    publishDir "${params.outdir}/${sample}/${job_id}", mode: 'copy'

    input:
    tuple val(job_id),
        val(sample), 
        path(input_folder), 
        val(observed),
        path(dists),
        path(tree)

    output:
    tuple val(job_id),
        val(sample), 
        path("annotated_tree.pickle"), emit: annotated_tree
    
    script:
    """
    python ${baseDir}/bin/process_tree/annotate_tree.py \
    --tree ${tree} \
    --X_bin ${input_folder}/X_bin.npz \
    --AD ${input_folder}/AD.npz \
    --DP ${input_folder}/DP.npz \
    --cell_meta ${input_folder}/cell_meta.csv \
    --char_meta ${input_folder}/char_meta.csv \
    --dists ${dists}
    """

    stub:
    """
    touch annotated_tree.pickle
    """

}
