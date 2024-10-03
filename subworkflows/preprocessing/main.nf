// preprocess

// Include here
nextflow.enable.dsl = 2
include { MAESTER } from "./modules/maester.nf"
include { DISTANCES } from "./modules/distances.nf"
include { DISTANCE_METRICS } from "./modules/distance_metrics.nf"
include { publish_pp } from "./modules/publish.nf"

//
 
//----------------------------------------------------------------------------//
// preprocess subworkflow
//----------------------------------------------------------------------------//

// preprocess subworkflow
workflow preprocess {
    
    take:
        ch_jobs 

    main:
     
        // Process input channel based on tracer system
        if (params.tracing_system == "MAESTER") {
            input_folder = MAESTER(ch_jobs)   
        } 
        
        // ...
        
        else {
            println('Only valid option is MAESTER so far...')
            // println('Provide valid tracing system option! (e.g., MAESTER, RedeeM, scmtATAC, Cas9, scWGS)')
        }

        // Calculate distances
        replicates = Channel.of( 1..(params.n_boot_distances-1) ).concat(Channel.of( "observed" )) 
        DISTANCES(input_folder.combine(replicates))
        ch_input = input_folder.combine(DISTANCES.out.distances, by:[0,1,2,3]).groupTuple(by: [0,1,2,3,4])
        DISTANCE_METRICS(ch_input)
        publish_pp(DISTANCE_METRICS.out.distance_metrics)

    emit:

        input = ch_input
        
} 
