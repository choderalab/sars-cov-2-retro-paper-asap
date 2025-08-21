process CROSS_DOCK_BY_LIGAND {
    publishDir "${params.dockedFiles}/${posit_method}_${num_poses}_poses", mode: 'link', overwrite: true
    conda "${params.drugforge}"
    tag "cross-dock ${compound_name}"
    clusterOptions '--partition "cpu" --cpus-per-task=1'
    memory 128.GB
    time 100.h

    input:
    tuple path(input_dir), path(prepped_dir), val(compound_name), path(ligandFile2d)
    val posit_method
    val selector
    val num_poses

    output:
    path("*_docked/*"), emit: docked

    script:
    """
    asap-docking cross-docking \
    --target SARS-CoV-2-Mpro \
    --use-omega \
    --omega-dense \
    --allow-retries \
    --allow-final-clash \
    --relax-mode clash \
    --posit-method "${posit_method}" \
    --structure-selector "${selector}" \
    --fragalysis-dir ${input_dir} \
    --ligands "${ligandFile2d}" \
    --cache-dir "${prepped_dir}" \
    --output-dir "${compound_name}_docked" \
    --overwrite \
    --no-save-to-cache \
    --use-only-cache \
    --num-poses "${num_poses}" \
    """
}

process CROSS_DOCK_BY_LIGAND_MULTIPOSE {
    publishDir "${params.dockedFiles}/${posit_method}_${num_poses}_poses", mode: 'link', overwrite: true
    conda "${params.drugforge}"
    tag "cross-dock ${compound_name}"
    clusterOptions '--partition "cpu" --cpus-per-task=1'
    memory 128.GB
    time 100.h

    input:
    tuple path(input_dir), path(prepped_dir), val(compound_name), path(ligandFile2d)
    val posit_method
    val selector
    val num_poses

    output:
    path("*_docked/*"), emit: docked

    script:
    """
    asap-docking cross-docking \
    --target SARS-CoV-2-Mpro \
    --use-omega \
    --omega-dense \
    --allow-final-clash \
    --posit-method "${posit_method}" \
    --structure-selector "${selector}" \
    --fragalysis-dir ${input_dir} \
    --ligands "${ligandFile2d}" \
    --cache-dir "${prepped_dir}" \
    --output-dir "${compound_name}_docked" \
    --overwrite \
    --no-save-to-cache \
    --use-only-cache \
    --num-poses "${num_poses}" \
    """
}