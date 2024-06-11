// MPBOOT

// Include here
nextflow.enable.dsl = 2
include { BUILD_MPBOOT } from "./modules/build_mpboot.nf"

// 
 
//----------------------------------------------------------------------------//
// MPBOOT subworkflow
//----------------------------------------------------------------------------//

//

// MPBOOT subworkflow
workflow MPBOOT {
    
    take:
        ch_input   

    main: 
 
        BUILD_MPBOOT(ch_input)

    emit:
        
        tree = BUILD_MPBOOT.out.tree

}



  