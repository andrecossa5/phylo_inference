// BUILD TREE module

nextflow.enable.dsl = 2


process BUILD_TREE {

    input:
    tuple val(sample), val(filtering), val(boot_replicate)
    
    output:
    path ("${sample}_${filtering}_rep${boot_replicate}.txt"), emit: tree
    
    script:
    """
    echo aaa > ${sample}_${filtering}_rep${boot_replicate}.txt
    """

    stub:
    """
    touch ${sample}_${filtering}_rep${boot_replicate}.txt
    """

}
