// hyper_tuning

// Include here
nextflow.enable.dsl = 2
include { ONESAMPLE } from "./modules/one_sample.nf"

//
 
//----------------------------------------------------------------------------//
// hyper_tuning subworkflow
//----------------------------------------------------------------------------//

workflow hyper_tuning {
    
    take:
        ch_input

    main: 
      
        // ch_input = ch_afm
        // ONESAMPLE(ch_input)

    emit:
        stats = ch_input
        // stats = ONESAMPLE.out.stats
        
} 
