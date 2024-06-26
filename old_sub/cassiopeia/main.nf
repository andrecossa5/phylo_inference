// CASSIOPEIA

// Include here
nextflow.enable.dsl = 2
include { CASSIOPEIA } from "./modules/cassiopeia.nf"

// 
 
//----------------------------------------------------------------------------//
// CASSIOPEIA subworkflow
//----------------------------------------------------------------------------//

//

// CASSIOPEIA subworkflow
workflow CASSIOPEIA {
    
    take:
        ch_input   

    main: 

        

 
    emit:

        tree = SUPPORT.out.tree

}


 