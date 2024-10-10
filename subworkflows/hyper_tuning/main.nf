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

        def index = 0
        ch_input = ch_input.map{ it -> tuple(it[1], it[2]) }
                .combine(Channel.fromList(params.min_n_positive))
                .combine(Channel.fromList(params.af_confident_detection))
                .combine(Channel.fromList(params.min_n_confidently_detected))
                .combine(Channel.fromList(params.min_mean_AD_in_positives))
                .combine(Channel.fromList(params.min_AD))
                .combine(Channel.fromList(params.bin_method))
                .combine(Channel.fromList(params.t_prob))
                .combine(Channel.fromList(params.min_cell_prevalence))
                .map{
                    def (a, b, c, d, e, f, g, h, i, l) = it
                    tuple(a, b, c, d, e, f, g, h, i, l, index++)
                }

        ONESAMPLE(ch_input)

    emit:
        stats = ONESAMPLE.out.stats
        
} 
