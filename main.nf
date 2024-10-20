// phylo_inference FINAL VERSION

nextflow.enable.dsl = 2
include { preprocess } from "./subworkflows/preprocessing/main"
include { build_tree } from "./subworkflows/tree_building/main"
include { process_tree } from "./subworkflows/tree_processing/main"

//

// Input channel
ch_jobs = Channel.fromPath(params.path_input)
            .splitCsv(header: true)
            .map { row -> [ row.job_id, row.sample, row.ch_matrix, row.char_filtering_key,
                            row.cell_filtering_key, row.bin_key, row.tree_key, row.cell_file ]}

//

//----------------------------------------------------------------------------//
// phylo_inference workflow main entry point
//----------------------------------------------------------------------------//

workflow phylo {

    preprocess(ch_jobs)
    build_tree(preprocess.out.input)
    process_tree(preprocess.out.input, build_tree.out.final_tree)
    process_tree.out.metrics.view()

}


// Mock

// Default message
workflow {
    
    println "\n"
    println "Hi there! This is the new version of the phylo_inference toolkit. The phylo entry point is currently supported."
    println "Usage: nextflow run main.nf -c <config> -params-file <params> -profile <profile> -entry phylo"
    println "See https://github.com/.../main.nf ./config and ./params for configurations and options available."
    println "N.B. BETA version under active development."
    println "\n"

}
