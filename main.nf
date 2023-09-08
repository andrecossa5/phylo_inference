// phylo_inference pipeline
nextflow.enable.dsl = 2
include { cassiopeia_workflow } from "./subworkflows/cassiopeia/main"

// Samples channel
ch_samples = Channel
    .fromPath("${params.path_data}/*", type:'dir') 
    .map{ it.getName() }

//

//----------------------------------------------------------------------------//
// phylo_inference entry points
//----------------------------------------------------------------------------//

//

workflow cassiopeia {

    cassiopeia_workflow(ch_samples)
    cassiopeia_workflow.out.results.view()

}

//

// Mock
workflow  {
    
    Channel.of(1,2,3,4) | view

}