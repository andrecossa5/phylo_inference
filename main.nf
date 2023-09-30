// phylo_inference pipeline
nextflow.enable.dsl = 2
include { tree_inference_workflow } from "./subworkflows/tree_inference/main"

// Samples channel
ch_samples = Channel
    .fromPath("${params.path_data}/*", type:'dir') 
    .map{ it.getName() }

//

//----------------------------------------------------------------------------//
// phylo_inference entry points
//----------------------------------------------------------------------------//

//

workflow phylo {

    tree_inference_workflow(ch_samples)
    tree_inference_workflow.out.results.view()

}

//

// Mock
workflow  {
    
    Channel.of(1,2,3,4) | view

}