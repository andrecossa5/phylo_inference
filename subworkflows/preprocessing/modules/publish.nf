// publish_pp module

nextflow.enable.dsl = 2 

//

process publish_pp {

    tag "${sample}: ${job_id}"

    publishDir "${params.outdir}/${sample}/${job_id}", mode: 'copy'

    input:
    tuple val(job_id),
        val(sample), 
        path(input_folder),
        path(distance_metrics)

    output:
    tuple val(job_id),
        val(sample), 
        path("bin_ops.csv"),
        path("cell_filtering_ops.csv"),
        path("char_filtering_ops.csv"),
        path("char_meta.csv"),
        path("dataset_df.csv"),
        path("tree_ops.csv"),
        path("distance_metrics.csv"), emit: pp

        // path("genotypes.fa")
        // path("cell_meta.csv"),
        // path("cell_meta.csv"),
        // path("X_bin.npz"),
        // path("AD.npz"),
        // path("DP.npz"),

    script:
    """
    cp ${input_folder}/* .
    """

    stub:
    """
    # touch AD.npz
    # touch DP.npz
    # touch X_bin.npz
    touch bin_ops.csv
    touch cell_filtering_ops.csv
    # touch cell_meta.csv
    touch char_filtering_ops.csv
    touch char_meta.csv
    touch dataset_df.csv
    touch tree_ops.csv
    touch distance_metrics.csv
    # touch genotypes.fa
    """

}

