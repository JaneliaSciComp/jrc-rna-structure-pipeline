#!/bin/bash
# Get the relative path to the script directory from the current script even if link or run from another directory
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
SCRIPT_DIR=$(realpath "${SCRIPT_DIR}/..")
# All of those create outputs in the current directory
# Get all RNA structures from PDB in the given date range
#
START="1978-01-01"
STOP="2025-12-17"
SPLIT_DATE="2025-05-29"
VERSION="v330"


RNA3DDB_CLUSTER="${SCRIPT_DIR}/../rna3db/rna3db-jsons/cluster.json"
RNA3DHUB_CSV="${SCRIPT_DIR}/../rna_3d_hub/nrlist_4.16_all.csv"

# Fail if any command fails
set -e

# Step control flags - set to 1 to enable, 0 to disable
STEP_FETCH_DATA=${STEP_FETCH_DATA:-1}
STEP_EXTRACT_METADATA=${STEP_EXTRACT_METADATA:-1}
STEP_GROUP_SEQUENCES=${STEP_GROUP_SEQUENCES:-1}
STEP_EXTRACT_NONREDUNDANT=${STEP_EXTRACT_NONREDUNDANT:-1}
STEP_CLUSTER_MMSEQS=${STEP_CLUSTER_MMSEQS:-1}
STEP_CALCULATE_STRUCTUREDNESS=${STEP_CALCULATE_STRUCTUREDNESS:-1}
STEP_ADD_STRUCTUREDNESS=${STEP_ADD_STRUCTUREDNESS:-1}
STEP_FILTER_METADATA=${STEP_FILTER_METADATA:-1}
STEP_SELECT_REPRESENTATIVES=${STEP_SELECT_REPRESENTATIVES:-1}
STEP_EXTRACT_COORDINATES=${STEP_EXTRACT_COORDINATES:-1}
STEP_SPLIT_TRAIN_TEST=${STEP_SPLIT_TRAIN_TEST:-1}
STEP_EXTRACT_BIOASSEMBLY=${STEP_EXTRACT_BIOASSEMBLY:-1}
STEP_KAGGLE=${STEP_KAGGLE:-1}

# DISABLE ALL STEPS FOR TESTING
# STEP_FETCH_DATA=0
# STEP_EXTRACT_METADATA=0
# STEP_GROUP_SEQUENCES=0
# STEP_EXTRACT_NONREDUNDANT=0
# STEP_CLUSTER_MMSEQS=0
# STEP_CALCULATE_STRUCTUREDNESS=0
# STEP_ADD_STRUCTUREDNESS=0
# STEP_FILTER_METADATA=0
# STEP_SELECT_REPRESENTATIVES=0
# STEP_EXTRACT_COORDINATES=0
# STEP_SPLIT_TRAIN_TEST=0
# STEP_EXTRACT_BIOASSEMBLY=0
# STEP_KAGGLE=0

# Parse command line arguments
show_help() {
    cat << EOF
Usage: $(basename "$0") [OPTIONS]

RNA structure processing workflow v330

Options:
  --help                Show this help message
  --dry-run             Show which steps will be run without executing them

  Individual step control:
  --skip-fetch          Skip fetching data from PDB
  --skip-extract-metadata   Skip metadata extraction
  --skip-group          Skip sequence grouping
  --skip-nonredundant   Skip non-redundant sequence extraction
  --skip-cluster        Skip MMseqs2 clustering
  --skip-calc-struct    Skip structuredness calculation
  --skip-add-struct     Skip adding structuredness to metadata
  --skip-filter         Skip metadata filtering
  --skip-select         Skip representative selection
  --skip-coords         Skip coordinate extraction
  --skip-split          Skip train/test split
  --skip-bioassembly-meta    Skip bioassembly metadata extraction
  --skip-kaggle        Skip Kaggle format conversion

  Group control (skip multiple steps):
  --skip-metadata       Skip metadata group (extract-metadata, group, cluster)
  --skip-structuredness Skip structuredness group (calc-struct, add-struct)
  --skip-filter-group   Skip filter group (filter, select)
  --skip-coordinates    Skip coordinates group (coords)
  --skip-split-group    Skip split group (split)
  --skip-bioassembly    Skip all bioassembly extraction

  Only run specific groups:
  --only-fetch          Only run data fetching
  --only-metadata       Only run metadata group
  --only-structuredness Only run structuredness group
  --only-filter         Only run filter group
  --only-coordinates    Only run coordinates group
  --only-split          Only run split group
  --only-bioassembly    Only run bioassembly group
  --only-kaggle         Only run Kaggle format conversion
EOF
}

# Initialize dry run flag
DRY_RUN=0

# Process command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --help)
            show_help
            exit 0
            ;;
        --dry-run)
            DRY_RUN=1
            ;;
        --skip-fetch)
            STEP_FETCH_DATA=0
            ;;
        --skip-extract-metadata)
            STEP_EXTRACT_METADATA=0
            ;;
        --skip-group)
            STEP_GROUP_SEQUENCES=0
            ;;
        --skip-nonredundant)
            STEP_EXTRACT_NONREDUNDANT=0
            ;;
        --skip-cluster)
            STEP_CLUSTER_MMSEQS=0
            ;;
        --skip-calc-struct)
            STEP_CALCULATE_STRUCTUREDNESS=0
            ;;
        --skip-add-struct)
            STEP_ADD_STRUCTUREDNESS=0
            ;;
        --skip-filter)
            STEP_FILTER_METADATA=0
            ;;
        --skip-select)
            STEP_SELECT_REPRESENTATIVES=0
            ;;
        --skip-coords)
            STEP_EXTRACT_COORDINATES=0
            ;;
        --skip-split)
            STEP_SPLIT_TRAIN_TEST=0
            ;;
        --skip-bioassembly-meta)
            STEP_EXTRACT_BIOASSEMBLY=0
            ;;
        --skip-kaggle)
            STEP_KAGGLE=0
            ;;
        --skip-metadata)
            STEP_EXTRACT_METADATA=0
            STEP_GROUP_SEQUENCES=0
            STEP_CLUSTER_MMSEQS=0
            ;;
        --skip-structuredness)
            STEP_CALCULATE_STRUCTUREDNESS=0
            STEP_ADD_STRUCTUREDNESS=0
            ;;
        --skip-filter-group)
            STEP_FILTER_METADATA=0
            STEP_SELECT_REPRESENTATIVES=0
            ;;
        --skip-coordinates)
            STEP_EXTRACT_COORDINATES=0
            ;;
        --skip-split-group)
            STEP_SPLIT_TRAIN_TEST=0
            ;;
        --skip-bioassembly)
            STEP_EXTRACT_BIOASSEMBLY=0
            ;;
        --only-fetch)
            STEP_FETCH_DATA=1
            STEP_EXTRACT_METADATA=0
            STEP_GROUP_SEQUENCES=0
            STEP_EXTRACT_NONREDUNDANT=0
            STEP_CLUSTER_MMSEQS=0
            STEP_CALCULATE_STRUCTUREDNESS=0
            STEP_ADD_STRUCTUREDNESS=0
            STEP_FILTER_METADATA=0
            STEP_SELECT_REPRESENTATIVES=0
            STEP_EXTRACT_COORDINATES=0
            STEP_SPLIT_TRAIN_TEST=0
            STEP_EXTRACT_BIOASSEMBLY=0
            STEP_KAGGLE=0
            ;;
        --only-metadata)
            STEP_FETCH_DATA=0
            STEP_EXTRACT_METADATA=1
            STEP_GROUP_SEQUENCES=1
            STEP_EXTRACT_NONREDUNDANT=0
            STEP_CLUSTER_MMSEQS=1
            STEP_CALCULATE_STRUCTUREDNESS=0
            STEP_ADD_STRUCTUREDNESS=0
            STEP_FILTER_METADATA=0
            STEP_SELECT_REPRESENTATIVES=0
            STEP_EXTRACT_COORDINATES=0
            STEP_SPLIT_TRAIN_TEST=0
            STEP_EXTRACT_BIOASSEMBLY=0
            STEP_KAGGLE=0
            ;;
        --only-structuredness)
            STEP_FETCH_DATA=0
            STEP_EXTRACT_METADATA=0
            STEP_GROUP_SEQUENCES=0
            STEP_EXTRACT_NONREDUNDANT=0
            STEP_CLUSTER_MMSEQS=0
            STEP_CALCULATE_STRUCTUREDNESS=1
            STEP_ADD_STRUCTUREDNESS=1
            STEP_FILTER_METADATA=0
            STEP_SELECT_REPRESENTATIVES=0
            STEP_EXTRACT_COORDINATES=0
            STEP_SPLIT_TRAIN_TEST=0
            STEP_EXTRACT_BIOASSEMBLY=0
            STEP_KAGGLE=0
            ;;
        --only-filter)
            STEP_FETCH_DATA=0
            STEP_EXTRACT_METADATA=0
            STEP_GROUP_SEQUENCES=0
            STEP_EXTRACT_NONREDUNDANT=0
            STEP_CLUSTER_MMSEQS=0
            STEP_CALCULATE_STRUCTUREDNESS=0
            STEP_ADD_STRUCTUREDNESS=0
            STEP_FILTER_METADATA=1
            STEP_SELECT_REPRESENTATIVES=1
            STEP_EXTRACT_COORDINATES=0
            STEP_SPLIT_TRAIN_TEST=0
            STEP_EXTRACT_BIOASSEMBLY=0
            STEP_KAGGLE=0
            ;;
        --only-coordinates)
            STEP_FETCH_DATA=0
            STEP_EXTRACT_METADATA=0
            STEP_GROUP_SEQUENCES=0
            STEP_EXTRACT_NONREDUNDANT=0
            STEP_CLUSTER_MMSEQS=0
            STEP_CALCULATE_STRUCTUREDNESS=0
            STEP_ADD_STRUCTUREDNESS=0
            STEP_FILTER_METADATA=0
            STEP_SELECT_REPRESENTATIVES=0
            STEP_EXTRACT_COORDINATES=1
            STEP_SPLIT_TRAIN_TEST=0
            STEP_EXTRACT_BIOASSEMBLY=0
            STEP_KAGGLE=0
            ;;
        --only-split)
            STEP_FETCH_DATA=0
            STEP_EXTRACT_METADATA=0
            STEP_GROUP_SEQUENCES=0
            STEP_EXTRACT_NONREDUNDANT=0
            STEP_CLUSTER_MMSEQS=0
            STEP_CALCULATE_STRUCTUREDNESS=0
            STEP_ADD_STRUCTUREDNESS=0
            STEP_FILTER_METADATA=0
            STEP_SELECT_REPRESENTATIVES=0
            STEP_EXTRACT_COORDINATES=0
            STEP_SPLIT_TRAIN_TEST=1
            STEP_EXTRACT_BIOASSEMBLY=0
            STEP_KAGGLE=0
            ;;
        --only-bioassembly)
            STEP_FETCH_DATA=0
            STEP_EXTRACT_METADATA=0
            STEP_GROUP_SEQUENCES=0
            STEP_EXTRACT_NONREDUNDANT=0
            STEP_CLUSTER_MMSEQS=0
            STEP_CALCULATE_STRUCTUREDNESS=0
            STEP_ADD_STRUCTUREDNESS=0
            STEP_FILTER_METADATA=0
            STEP_SELECT_REPRESENTATIVES=0
            STEP_EXTRACT_COORDINATES=0
            STEP_SPLIT_TRAIN_TEST=0
            STEP_EXTRACT_BIOASSEMBLY=1
            STEP_KAGGLE=0
            ;;
        --only-kaggle)
            STEP_FETCH_DATA=0
            STEP_EXTRACT_METADATA=0
            STEP_GROUP_SEQUENCES=0
            STEP_EXTRACT_NONREDUNDANT=0
            STEP_CLUSTER_MMSEQS=0
            STEP_CALCULATE_STRUCTUREDNESS=0
            STEP_ADD_STRUCTUREDNESS=0
            STEP_FILTER_METADATA=0
            STEP_SELECT_REPRESENTATIVES=0
            STEP_EXTRACT_COORDINATES=0
            STEP_SPLIT_TRAIN_TEST=0
            STEP_EXTRACT_BIOASSEMBLY=0
            STEP_KAGGLE=1
            ;;  

        *)
            echo "Unknown option: $1"
            show_help
            exit 1
            ;;
    esac
    shift
done

NAME_PREFIX="${VERSION}_${START}_${STOP}"
TEST_NAME_PREFIX="${VERSION}_${SPLIT_DATE}_${STOP}"
TRAIN_NAME_PREFIX="${VERSION}_${START}_${SPLIT_DATE}"

echo "Using name prefix: ${NAME_PREFIX}"
echo "Enabled steps:"
echo "  STEP_FETCH_DATA: ${STEP_FETCH_DATA}"
echo "  STEP_EXTRACT_METADATA: ${STEP_EXTRACT_METADATA}"
echo "  STEP_GROUP_SEQUENCES: ${STEP_GROUP_SEQUENCES}"
echo "  STEP_EXTRACT_NONREDUNDANT: ${STEP_EXTRACT_NONREDUNDANT}"
echo "  STEP_CLUSTER_MMSEQS: ${STEP_CLUSTER_MMSEQS}"
echo "  STEP_CALCULATE_STRUCTUREDNESS: ${STEP_CALCULATE_STRUCTUREDNESS}"
echo "  STEP_ADD_STRUCTUREDNESS: ${STEP_ADD_STRUCTUREDNESS}"
echo "  STEP_FILTER_METADATA: ${STEP_FILTER_METADATA}"
echo "  STEP_SELECT_REPRESENTATIVES: ${STEP_SELECT_REPRESENTATIVES}"
echo "  STEP_EXTRACT_COORDINATES: ${STEP_EXTRACT_COORDINATES}"
echo "  STEP_SPLIT_TRAIN_TEST: ${STEP_SPLIT_TRAIN_TEST}"
echo "  STEP_EXTRACT_BIOASSEMBLY: ${STEP_EXTRACT_BIOASSEMBLY}"
echo "  STEP_KAGGLE: ${STEP_KAGGLE}"
echo ""

if [ ${DRY_RUN} -eq 1 ]; then
    echo "DRY RUN MODE - No commands will be executed"
    echo ""
    echo "The following steps will be run:"
    [ ${STEP_FETCH_DATA} -eq 1 ] && echo "  ✓ Fetch data from PDB"
    [ ${STEP_EXTRACT_METADATA} -eq 1 ] && echo "  ✓ Extract metadata"
    [ ${STEP_GROUP_SEQUENCES} -eq 1 ] && echo "  ✓ Group sequences"
    [ ${STEP_EXTRACT_NONREDUNDANT} -eq 1 ] && echo "  ✓ Extract non-redundant sequences"
    [ ${STEP_CLUSTER_MMSEQS} -eq 1 ] && echo "  ✓ Cluster with MMseqs2"
    [ ${STEP_CALCULATE_STRUCTUREDNESS} -eq 1 ] && echo "  ✓ Calculate structuredness"
    [ ${STEP_ADD_STRUCTUREDNESS} -eq 1 ] && echo "  ✓ Add structuredness to metadata"
    [ ${STEP_FILTER_METADATA} -eq 1 ] && echo "  ✓ Filter metadata"
    [ ${STEP_SELECT_REPRESENTATIVES} -eq 1 ] && echo "  ✓ Select representatives"
    [ ${STEP_EXTRACT_COORDINATES} -eq 1 ] && echo "  ✓ Extract coordinates"
    [ ${STEP_SPLIT_TRAIN_TEST} -eq 1 ] && echo "  ✓ Split train/test"
    [ ${STEP_EXTRACT_BIOASSEMBLY} -eq 1 ] && echo "  ✓ Extract bioassembly"
    [ ${STEP_KAGGLE} -eq 1 ] && echo "  ✓ Convert to Kaggle format"
    echo ""
    echo "Skipped steps:"
    [ ${STEP_FETCH_DATA} -eq 0 ] && echo "  ✗ Fetch data from PDB"
    [ ${STEP_EXTRACT_METADATA} -eq 0 ] && echo "  ✗ Extract metadata"
    [ ${STEP_GROUP_SEQUENCES} -eq 0 ] && echo "  ✗ Group sequences"
    [ ${STEP_EXTRACT_NONREDUNDANT} -eq 0 ] && echo "  ✗ Extract non-redundant sequences"
    [ ${STEP_CLUSTER_MMSEQS} -eq 0 ] && echo "  ✗ Cluster with MMseqs2"
    [ ${STEP_CALCULATE_STRUCTUREDNESS} -eq 0 ] && echo "  ✗ Calculate structuredness"
    [ ${STEP_ADD_STRUCTUREDNESS} -eq 0 ] && echo "  ✗ Add structuredness to metadata"
    [ ${STEP_FILTER_METADATA} -eq 0 ] && echo "  ✗ Filter metadata"
    [ ${STEP_SELECT_REPRESENTATIVES} -eq 0 ] && echo "  ✗ Select representatives"
    [ ${STEP_EXTRACT_COORDINATES} -eq 0 ] && echo "  ✗ Extract coordinates"
    [ ${STEP_SPLIT_TRAIN_TEST} -eq 0 ] && echo "  ✗ Split train/test"
    [ ${STEP_EXTRACT_BIOASSEMBLY} -eq 0 ] && echo "  ✗ Extract bioassembly"
    [ ${STEP_KAGGLE} -eq 0 ] && echo "  ✗ Convert to Kaggle format"
    exit 0
fi

if [ ${STEP_FETCH_DATA} -eq 1 ]; then
echo '================================'
echo "Fetching data from ${START} to ${STOP}"
python "${SCRIPT_DIR}/get_pdb_data.py" --start_date ${START}  --stop_date ${STOP} --search_type polymer_type --search_term RNA NA-hybrid --output_dir rna_structures_all --cache_dir ${SCRIPT_DIR}/output/all_cache
echo '================================'
fi

if [ ${STEP_EXTRACT_METADATA} -eq 1 ]; then
echo '================================'
echo "Extract metadata for all RNA structures"
echo '================================'
python ${SCRIPT_DIR}/metadata_extraction/extract_metadata.py ./rna_structures_all/ metadata.csv --keyword ribosom ribosome --multiprocessing
echo '================================'

echo '================================'
echo "Extract ligands SMILES for all RNA structures"
echo '================================'
python ${SCRIPT_DIR}/metadata_extraction/extract_smiles.py rna_structures_all/ ligands.csv --multiprocessing --smiles-file ${SCRIPT_DIR}/data/Components-smiles-stereo-oe_20251221.smi
echo '================================'
fi

if [ ${STEP_GROUP_SEQUENCES} -eq 1 ]; then
echo '================================'
echo "Group RNA structures by sequence identity to remove redundancy at 100% identity (and 90% coverage)"
echo '================================'
python "${SCRIPT_DIR}/groupby_sequence_csv.py" --group_output_name=group_id --min_match_ratio=0.9 ./metadata.csv ./metadata_grouped1.csv
python "${SCRIPT_DIR}/groupby_sequence_csv.py" --group_output_name=seq_group_id --min_match_ratio=1.0 ./metadata_grouped1.csv ./metadata_grouped_full.csv
echo '================================'
fi

if [ ${STEP_EXTRACT_NONREDUNDANT} -eq 1 ]; then
echo '================================'
echo "Extract non-redundant sequences from grouped metadata"
echo '================================'
python ${SCRIPT_DIR}/tools/csv_to_fasta.py --id-column target_id --group-by group_id metadata_grouped_full.csv -o nonredundant_sequences.fasta
echo '================================'
fi

if [ ${STEP_CLUSTER_MMSEQS} -eq 1 ]; then
echo '================================'
echo "Cluster RNA structures using MMseqs2 at different identity thresholds"
echo '================================'
# clean previous runs
if [ -d mmseqs_clustering ]; then
    rm -rf mmseqs_clustering
fi
mkdir -p mmseqs_clustering
(cd mmseqs_clustering 
  "${SCRIPT_DIR}/tools/cluster_mmseqs.sh" -i ../nonredundant_sequences.fasta
  "${SCRIPT_DIR}/tools/process_clustering_rounds.py" --cluster-files $(ls *.tsv | sort -r | tr '\n' ' ') --initial-members ../metadata_grouped_full.csv --output merged_clusters.csv --rename-columns target_id group_id $(ls *.tsv | sort -r | sed 's/_clusters.*//' | tr '\n' ' ')

)
echo "Merge files"
echo '================================'


python "${SCRIPT_DIR}/metadata_extraction/merge_csv_files.py" metadata_grouped_full.csv mmseqs_clustering/merged_clusters.csv metadata_merged_clustering.csv --on 'target_id' 'group_id'
python "${SCRIPT_DIR}/metadata_extraction/add_rna3ddb_mapping.py" metadata_merged_clustering.csv ${RNA3DDB_CLUSTER} 
python "${SCRIPT_DIR}/metadata_extraction/add_rna3dhub_mapping.py" metadata_merged_clustering_with_rna3ddb.csv ${RNA3DHUB_CSV}
fi

if [ ${STEP_CALCULATE_STRUCTUREDNESS} -eq 1 ]; then
echo '================================' 
echo "Calculate structuredness for all files"
echo '================================'
# Adjust for local resources
( 
  mkdir -p processing && cd processing 
  find ../rna_structures_all/ -name '*.cif' |  parallel --memsuspend 30G --resume --bar -j 60 --joblog job.log "ls {/.}_structuredness.csv > /dev/null 2> /dev/null || python ${SCRIPT_DIR}/processing/structuredness.py {} {/.}_structuredness.csv"
)
echo '================================' 
fi

if [ ${STEP_ADD_STRUCTUREDNESS} -eq 1 ]; then
echo "Add structuredness metrics to metadata"
echo '================================'
python ${SCRIPT_DIR}/merge_metadata.py metadata_merged_clustering_with_rna3ddb_with_rna3dhub.csv 'processing/*.csv' jrc_${VERSION}_metadata_with_structuredness.csv --join-columns pdb_id chain_id
python ${SCRIPT_DIR}/postprocess/001_add_length_adjust_structuredness.py jrc_${VERSION}_metadata_with_structuredness.csv jrc_${VERSION}_metadata.csv
echo '================================'
fi

if [ ${STEP_FILTER_METADATA} -eq 1 ]; then
# Filter
echo '================================'
echo "Filter metadata based on quality criteria"
echo '================================'
duckdb  -c "copy (
        select *
        from \"jrc_${VERSION}_metadata.csv\"
        where 
                unmapped_nakb is false
            and undefined_residues is false
            and unexpected_residues is false
            and mapped_fraction <= 0.25 
            and fraction_observed >= 0.5
            and len(sequence) > 10
            and intra_chain_structuredness_adjusted > 0.2
) to \"jrc_${VERSION}_metadata_filtered.csv\" with csv header delimiter as ',';
"

# Bioassembly
# Filter for bioassemblies
echo '================================'
echo "Filter metadata for bioassembly based on quality criteria"
echo '================================'
duckdb  -c "
copy (
    with 
        -- filter out short RNA chains
        rna as 
            (select *
                from \"jrc_${VERSION}_metadata.csv\"
            ),
        -- filter out PDB ids that have poor RNA chains
        poorly_defined_pdb as 
            (select distinct pdb_id
                from rna
                where 
                    unmapped_nakb is true
                or undefined_residues is true
                or unexpected_residues is true
                or mapped_fraction > 0.25 
                or fraction_observed < 0.5
                or total_structuredness_adjusted < 0.2
            )
    select *
        from rna
        where pdb_id not in (select pdb_id from poorly_defined_pdb)
) to \"jrc_${VERSION}_metadata_for_bioassembly.csv\" with csv header delimiter as ',';
"
fi

if [ ${STEP_EXTRACT_BIOASSEMBLY} -eq 1 ]; then
echo '================================' 
echo "Extract bioassembly metadata and select single representative per sequence"
echo '================================'

python "${SCRIPT_DIR}/metadata_extraction/extract_bioassembly_metadata.py" rna_structures_all jrc_${VERSION}_metadata_for_bioassembly.csv jrc_${VERSION}_bioassembly_tmp.csv --pattern {pdb_id}-assembly1.cif.gz
python ../../metadata_extraction/merge_groupings.py jrc_${VERSION}_bioassembly_tmp.csv jrc_${VERSION}_metadata.csv jrc_${VERSION}_bioassembly_tmp2.csv

# Filter bioassemblys to only include multichain assemblies or single chains with high structuredness
# This is necessary due to bug where RNA structuredness is calculated with DNA.
duckdb -c "
copy (
    select * from (
            select *,
                string_split(stoichiometry, ';') as sparts
            from 'jrc_${VERSION}_bioassembly_tmp2.csv'
        )
    where len(sparts) > 1
        or (len(sparts) = 1 and string_split(sparts[1], ':')[2]::int >  1)
        or pdb_id in (
            select pdb_id
            from 'jrc_${VERSION}_metadata.csv'
            where intra_chain_structuredness_adjusted >= 0.2
        )
    ) to 'jrc_${VERSION}_bioassembly_filtered.csv' with csv header delimiter as ',';
"
fi

# Add group by combined sequence
python "${SCRIPT_DIR}/groupby_sequence_csv.py" --group_output_name=seq_group_id --min_match_ratio=1.0 ./jrc_${VERSION}_bioassembly_filtered.csv ./jrc_${VERSION}_bioassembly.csv

if [ ${STEP_EXTRACT_COORDINATES} -eq 1 ]; then
echo '================================'
echo "Extract coordinates for all RNA structures"
echo '================================'
python ${SCRIPT_DIR}/processing/extract_coordinates.py jrc_${VERSION}_bioassembly.csv rna_structures_all/ jrc_${VERSION}_bioassembly_coordinates.csv --file-pattern {pdb_id}-assembly1.cif.gz --multiple-chains --skip-mismatch-validation

fi

if [ ${STEP_SPLIT_TRAIN_TEST} -eq 1 ]; then
echo '================================'
echo "Split filtered metadata and coordinates into train/test by temporal cutoff at ${SPLIT_DATE}, grouped by 30% identity clusters"
echo '================================'

python "${SCRIPT_DIR}/split_by_cutoff.py" jrc_${VERSION}_metadata.csv jrc_${VERSION}_bioassembly.csv --cutoff_date ${SPLIT_DATE} --group_key "mmseqs_0.300" --reference_key='pdb_id'

for part in train test; do
    duckdb -c "copy (
        with coords as (
        select * from 'jrc_${VERSION}_bioassembly_coordinates.csv'
        ) 
        select * from coords where target_id in (
        select target_id from '${part}_jrc_${VERSION}_bioassembly.csv'
        ) 
    ) to '${part}_jrc_${VERSION}_bioassembly_coordinates.csv';"
done

echo '================================'
fi

if [ ${STEP_KAGGLE} -eq 1 ]; then
# Convert to Kaggle format
echo '================================'
echo "Convert to Kaggle format"
echo '================================'

KAGGLE_OUTPUT_DIR="kaggle_jrc_${VERSION}"
rm -rf ${KAGGLE_OUTPUT_DIR}
mkdir -p ${KAGGLE_OUTPUT_DIR}

# Prepare kaggle specfic coordinate files
# Test/validation set is non-redundant at 90% identity, but group solutions by identical sequences (100% identity) for coordinate files
python ${SCRIPT_DIR}/add_group_index.py test_jrc_${VERSION}_bioassembly_coordinates.csv test_jrc_${VERSION}_bioassembly_coordinates_gid.csv --reference test_jrc_${VERSION}_bioassembly.csv --group_index_name solution_id --sort_by temporal_cutoff target_id --group_key seq_group_id --sort_by_after seq_group_id solution_id ID

sed -E "1s/\\<([xyz])_1\\>/\"C1'_\\1\"/g" test_jrc_${VERSION}_bioassembly_coordinates_gid.csv > test_jrc_${VERSION}_bioassembly_coordinates_gid_c1p.csv 
python "${SCRIPT_DIR}/convert_to_multisolution_wide.py" "test_jrc_${VERSION}_bioassembly_coordinates_gid_c1p.csv" "test_jrc_${VERSION}_bioassembly_C1prime.csv" --select "C1'" --add-column Usage Public --group-key seq_group_id --drop-column seq_group_id --fill-coordinates="-1e18" --pad 40

python "${SCRIPT_DIR}/get_single_target_per_group.py" test_jrc_${VERSION}_bioassembly.csv test_jrc_${VERSION}_bioassembly_single_per_cluster.csv --group-key mmseqs_0.900 --multichain --sort_by temporal_cutoff target_id --merge description

# Create sequence files in Kaggle format 
# Filter out short sequences, less than 10 nucleotides
# Filter out sequences that don't have coordinates (missing C1' atoms)
# Add ligand SMILES 
duckdb -c "copy (
    select 
        m.pdb_id as target_id, 
        sequence, temporal_cutoff, 
        title as description,
        stoichiometry, 
        all_sequences,
        l.ligand_comp_id as ligand_ids, 
        l.ligand_SMILES as ligand_SMILES
    from 'train_jrc_${VERSION}_bioassembly.csv' as m
    left join 'ligands.csv' as l on m.pdb_id == l.pdb_id
    where len(sequence) >= 10
    and m.pdb_id in (
        select distinct target_id 
        from 'train_jrc_${VERSION}_bioassembly_coordinates.csv'
        )
    order by temporal_cutoff, target_id
    )
    to '${KAGGLE_OUTPUT_DIR}/train_sequences.csv'
"

# Filter bioassemblies to only include ones where RNA is more than 40% of total mass
# This is necessary reduce number of test/validation targets
# Add ligand SMILES
duckdb -c "copy (
    select 
        m.pdb_id as target_id, 
        sequence, temporal_cutoff, 
        title as description,
        stoichiometry, 
        all_sequences, 
        l.ligand_comp_id as ligand_ids, 
        l.ligand_SMILES as ligand_SMILES
    from 'test_jrc_${VERSION}_bioassembly_single_per_cluster.csv' as m
    left join 'ligands.csv' as l on m.pdb_id == l.pdb_id
    where len(sequence) >= 10 
    and m.pdb_id in (
            select pdb_id
            from 'jrc_${VERSION}_metadata.csv'
            where composition_rna_fraction >= 0.4 
        )
    order by temporal_cutoff, target_id
    )

    to '${KAGGLE_OUTPUT_DIR}/test_sequences.csv'
"
cp ${KAGGLE_OUTPUT_DIR}/test_sequences.csv ${KAGGLE_OUTPUT_DIR}/validation_sequences.csv

## Create label (coordinate) files in Kaggle format
duckdb -c "copy (
    with coords as (
    select * from 'test_jrc_${VERSION}_bioassembly_C1prime.csv'
    ) 
    select ID,resname,resid,COLUMNS('[xyz]_.*') from (
        select * from coords where ID[:4] in (
            select target_id from '${KAGGLE_OUTPUT_DIR}/test_sequences.csv'
        )
    )
) to '${KAGGLE_OUTPUT_DIR}/test_labels.csv';"

duckdb -c "copy (
    with coords as (
    select * from 'jrc_${VERSION}_bioassembly_coordinates.csv'
    ) 
    select ID,resname,resid,COLUMNS('[xyz]_.*') from (
        select * from coords where target_id in (
            select target_id from '${KAGGLE_OUTPUT_DIR}/train_sequences.csv'
        )
    )
) to '${KAGGLE_OUTPUT_DIR}/validation_labels.csv';"


# Validate Kaggle files
for part in train test; do
    fail=0
    echo "Validating ${part} Kaggle files"
    ${SCRIPT_DIR}/validation/validate_bioassembly_metadata.py ${KAGGLE_OUTPUT_DIR}/${part}_sequences.csv || fail=1
    ${SCRIPT_DIR}/validation/validate_solution.py ${KAGGLE_OUTPUT_DIR}/${part}_labels.csv ${KAGGLE_OUTPUT_DIR}/${part}_sequences.csv || fail=1
    if [ $fail -ne 0 ]; then
        echo "✗ Validation failed for ${part} Kaggle files"
    else
        echo "✓ ${part} Kaggle files validated successfully"
    fi
done

fi
