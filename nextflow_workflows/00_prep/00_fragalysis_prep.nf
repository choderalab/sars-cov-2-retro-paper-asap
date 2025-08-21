#!/usr/bin/env nextflow
include {
    PREP_FRAGALYSIS
} from "./modules.nf"

workflow {
    frag_ch = Channel.fromPath("${params.curatedFragalysis}", type: 'dir')
    PREP_FRAGALYSIS(frag_ch)
}
