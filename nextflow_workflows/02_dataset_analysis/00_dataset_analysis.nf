#!/usr/bin/env nextflow
include {
    GENERATE_DATE_DICTIONARY
    CALCULATE_ECFP_TANIMOTO
    CALCULATE_MCS_TANIMOTO
    CALCULATE_TANIMOTO_COMBO
    COMBINE_CHEMICAL_SIMILARITY_DATA
    RUN_BEMIS_MURCKO_CLUSTERING
} from "./modules.nf"

workflow {
    // Your input channel for ligand file
    ligand_file_3d = Channel.fromPath("${params.ligandFiles}/${params.ligandFile3d}")
    ligand_file_2d = Channel.fromPath("${params.ligandFiles}/${params.ligandFile2d}")

    // Run calculation processes and collect their CSV outputs
    ecfp_results = CALCULATE_ECFP_TANIMOTO(ligand_file_3d)
        .ecfp_tanimoto
        .map { dir -> dir.listFiles().find { it.name.endsWith('.csv') } }

    mcs_results = CALCULATE_MCS_TANIMOTO(ligand_file_3d)
        .mcs_tanimoto
        .map { dir -> dir.listFiles().find { it.name.endsWith('.csv') } }

    tanimoto_combo_results = CALCULATE_TANIMOTO_COMBO(ligand_file_3d)
        .tanimoto_combo
        .map { dir -> dir.listFiles().find { it.name.endsWith('.csv') } }

    // Combine all CSV files into a single channel
    all_csv_files = ecfp_results
        .mix(mcs_results)
        .mix(tanimoto_combo_results)
        .collect()

    // Pass collected CSV files to combine process
    COMBINE_CHEMICAL_SIMILARITY_DATA(all_csv_files)

    // Generate date dictionary
    GENERATE_DATE_DICTIONARY()

    // Run scaffolding
    RUN_BEMIS_MURCKO_CLUSTERING(ligand_file_2d)
}

// Entry point: Calculate ECFP Tanimoto similarity only
workflow ECFP_ANALYSIS {
    ligand_file_3d = Channel.fromPath("${params.ligandFiles}/${params.ligandFile3d}")
    CALCULATE_ECFP_TANIMOTO(ligand_file_3d)
}

// Entry point: Combine existing CSV files (assumes CSV files already exist)
workflow COMBINE_SIMILARITY_DATA {
    // This assumes CSV files already exist in the expected locations
    // You might need to adjust paths based on your directory structure
    csv_files = Channel.fromPath("${params.chemicalSimilarityData}/*/*.csv").collect()
    COMBINE_CHEMICAL_SIMILARITY_DATA(csv_files)
}