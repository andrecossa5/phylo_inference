// CASSIOPEIA module

nextflow.enable.dsl = 2


process CASSIOPEIA {

    tag "${sample}: ${filtering_key}, rep=${boot_replicate}"

    input: 
    tuple val(sample), 
        val(filtering_key), 
        path(input_folder),
        val(boot_replicate),

    output:
    tuple val(sample), 
        val(filtering_key), 
        path(input_folder),
        val(boot_replicate),
        path("rep_${boot_replicate}.newick"), emit: tree
    
    script:
    """
    python ${baseDir}/bin/cassiopeia/build_cassiopeia.py \
    -p ${input_folder} \
    --sample ${sample} \
    --solver ${params.cassiopeia_solver} \
    --metric ${params.metric} \
    --name rep_${boot_replicate} \
    --ncores ${task.cpus} \
    --path_filtering ${params.path_filtering} \
    --filtering_key ${filtering_key} \
    --path_meta ${params.path_meta} \
    --path_priors ${params.path_priors} \
    --lineage_column ${params.lineage_column}
    """

    stub:
    """
    touch rep_${boot_replicate}.newick
    """

}
