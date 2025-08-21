process GENERATE_DATE_DICTIONARY {
    publishDir "${params.dataPath}", mode: 'copy', overwrite: true
    conda "${params.asap}"
    tag "generate-date-dictionary"

    output:
    path "cmpd_date_dict"
    path "cmpd_date_dict/date_dict.json", emit: structure_to_date_dict
    path "cmpd_date_dict/structure_to_cmpd_dict.json", emit: structure_to_cmpd_dict

    script:
    """
    python3 "${params.scripts}"/generate_date_dict.py --fragalysis-dir "${params.curatedFragalysis}" --output-dir cmpd_date_dict
    """
}
process CALCULATE_ECFP_TANIMOTO {
    publishDir "${params.chemicalSimilarityData}", mode: 'copy', overwrite: true
    conda "${params.asap}"
    tag "calculate-ecfp-tanimoto"

    input:
    path(ligand_file_3d)

    output:
    path("ecfp_tanimoto"), emit: ecfp_tanimoto

    script:
    """
    python3 "${params.scripts}"/calculate_ecfp_tanimoto.py --ref-ligand-sdf "${ligand_file_3d}" --output-dir ecfp_tanimoto
    """
}
process CALCULATE_MCS_TANIMOTO {
    publishDir "${params.chemicalSimilarityData}", mode: 'copy', overwrite: true
    conda "${params.asap}"
    tag "calculate-mcs-tanimoto"
    clusterOptions '--partition "cpu" --time=24:00:00 --mem=64GB --cpus-per-task=32'

    input:
    path(ligand_file_3d)

    output:
    path("mcs_tanimoto"), emit: mcs_tanimoto

    script:
    """
    python3 "${params.scripts}"/calculate_mcs_tanimoto.py --ref-ligand-sdf "${ligand_file_3d}" --output-dir mcs_tanimoto --ncpus 32
    """
}
process CALCULATE_TANIMOTO_COMBO {
    publishDir "${params.chemicalSimilarityData}", mode: 'copy', overwrite: true
    conda "${params.asap}"
    tag "calculate-tanimoto-combo"
    clusterOptions '--partition "cpu" --time=06:00:00 --mem=64GB --cpus-per-task=32'

    input:
    path(ligand_file_3d)

    output:
    path("tanimoto_combo"), emit: tanimoto_combo

    script:
    """
    python3 "${params.scripts}"/calculate_tanimoto_combo.py --ref-ligand-sdf "${ligand_file_3d}" --output-dir tanimoto_combo
    """
}
process COMBINE_CHEMICAL_SIMILARITY_DATA {
    publishDir "${params.chemicalSimilarityData}", mode: 'copy', overwrite: true
    conda "${params.asap}"
    tag "combine-chemical-similarity-data"

    input:
    path csv_files

    output:
    path "combined_chemical_similarity_data.csv", emit: combined_chemical_similarity_data

    script:
    """
    python3 "${params.scripts}/combine_chemical_similarity_data.py" ${csv_files.join(' ')}
    """
}
process RUN_BEMIS_MURCKO_CLUSTERING {
    publishDir "${params.chemicalSimilarityData}", mode: 'copy', overwrite: true
    conda "${params.asap}"
    tag "run-bemis-murcko-clustering"

    input:
    path(ligand_file_2d)

    output:
    path "${params.scaffoldDataName}", emit: bemis_murcko_clustering

    script:
    """
    python "${params.scripts}"/run_bemis_murcko_clustering.py --sdf-2d ${ligand_file_2d} --output-dir "${params.scaffoldDataName}"
    """
}