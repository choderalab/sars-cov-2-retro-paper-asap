process PREP_FRAGALYSIS {
    publishDir "${params.dataPath}", mode: 'copy', overwrite: true, saveAs: {fn -> "${params.fragalysisCache}"}
    conda "${params.drugforge}"
    tag "prep-fragalysis"
    clusterOptions '-c 32 --mem=128G --time=04:00:00'

    input:
    path curatedFragalysis

    output:
    path("output"), emit: fragalysisCache

    script:
    """
    asap-cli protein-prep \
      --target SARS-CoV-2-Mpro \
      --fragalysis-dir "${curatedFragalysis}" \
      --loop-db ${params.loopDB} \
      --ref-chain A \
      --active-site-chain A \
      --use-dask \
      --dask-n-workers 32 \
      --dask-type local \
    """
}
process PREP_CACHE_FOR_DOCKING {
    publishDir "${params.dataPath}", mode: 'copy', overwrite: true, saveAs: {fn -> "${params.fixedFragalysisCache}"}
    conda "${params.drugforge}"
    tag "prep-cache-for-docking"
    label 'local'

    input:
    path cache_dir

    output:
    path "*_fixed", emit: fixed_cache

    script:
    """
    python3 "${params.scripts}"/prep_cache_for_docking.py --input_cache "${cache_dir}"
    """
}
process GENERATE_COMBINED_LIGAND_FILES {
    publishDir "${params.ligandFiles}", mode: 'copy', overwrite: true
    conda "${params.drugforge}"
    tag "generate-ligand-files"
    label 'local'

    input:
    path cache_dir

    output:
    path "${params.ligandFile3d}", emit: ligandFile3d
    path "${params.ligandFile2d}", emit: ligandFile2d

    script:
    """
    python3 ${params.scripts}/combined_sdf_from_cache.py --input_cache "${cache_dir}"
    python3 ${params.scripts}/combined_sdf_from_cache.py --input_cache "${cache_dir}" --flatten
    """
}
process GENERATE_SPLIT_LIGAND_FILES {
    publishDir "${params.ligandFiles}", mode: 'copy', overwrite: true
    conda "${params.drugforge}"
    tag "generate-ligand-files"
    label 'local'

    input:
    path(ligandFile3d)
    path(ligandFile2d)

    output:
    path "${params.split3dligandFiles}", emit: split3dligandFiles
    path "${params.split2dligandFiles}", emit: split2dligandFiles

    script:
    """
    python3 ${params.scripts}/split_sdf.py --sdf_fn ${ligandFile3d} --out_dir ${params.split3dligandFiles} --chunk_size 1 --name_convention "integer"
    python3 ${params.scripts}/split_sdf.py --sdf_fn ${ligandFile2d} --out_dir ${params.split2dligandFiles} --chunk_size 1 --name_convention "integer"
    """
}

process DEDUPLICATE_LIGANDS {
    publishDir "${params.dataPath}", mode: 'copy', overwrite: true, saveAs: {fn -> "${params.fixedFragalysisCache}"}
    conda "${params.drugforge}"
    tag "deduplicate-ligands"
    label 'local'

    input:
    path fragalysis_dir
    path prepped_path

    output:
    path "deduped_cache", emit: fixed_cache

    script:
    """
    python3 ${params.scripts}/deduplicate_ligands.py \
    --fragalysis-dir ${fragalysis_dir} \
    --prepped-path ${prepped_path} \
    --output-dir deduped_cache \
    --remove-covalent
    """
}