#!/bin/bash
# Get the relative path to the script directory from the current script even if link or run from another directory
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
SCRIPT_DIR=$(realpath "${SCRIPT_DIR}/..")
# All of those create outputs in the current directory
# Get all RNA structures from PDB in the given date range
#
START="1978-01-01"
STOP="2025-12-03"
SPLIT_DATE="2025-05-29"
VERSION="v320rc8"


RNA3DDB_CLUSTER="${SCRIPT_DIR}/../rna3db/rna3db-jsons/cluster.json"
RNA3DHUB_CSV="${SCRIPT_DIR}/../rna_3d_hub/nrlist_4.16_all.csv"

# Fail if any command fails
set -e

NAME_PREFIX="${VERSION}_${START}_${STOP}"
TEST_NAME_PREFIX="${VERSION}_${SPLIT_DATE}_${STOP}"
TRAIN_NAME_PREFIX="${VERSION}_${START}_${SPLIT_DATE}"
echo "Using name prefix: ${NAME_PREFIX}"
echo '================================'
echo "Fetching data from ${START} to ${STOP}"
python "${SCRIPT_DIR}/get_pdb_data.py" --start_date ${START}  --stop_date ${STOP} --search_type polymer_type --search_term RNA NA-hybrid --output_dir rna_structures_all --cache_dir ${SCRIPT_DIR}/output/all_cache
echo '================================'

echo '================================'
echo "Extract metadata for all RNA structures"
echo '================================'
python ${SCRIPT_DIR}/metadata_extraction/extract_metadata.py ./rna_structures_all/ metadata.csv --keyword ribosom ribosome --multiprocessing
echo '================================'


echo '================================'
echo "Group RNA structures by sequence identity to remove redundancy at 100% identity (and 90% coverage)"
echo '================================'
python "${SCRIPT_DIR}/groupby_sequence_csv.py" --group_output_name=group_id --min_match_ratio=0.9 ./metadata.csv ./metadata_grouped1.csv
python "${SCRIPT_DIR}/groupby_sequence_csv.py" --group_output_name=seq_group_id --min_match_ratio=1.0 ./metadata_grouped1.csv ./metadata_grouped_full.csv
echo '================================'

echo '================================'
echo "Extract non-redundant sequences from grouped metadata"
echo '================================'
python ${SCRIPT_DIR}/tools/csv_to_fasta.py --id-column target_id --group-by group_id metadata_grouped_full.csv -o nonredundant_sequences.fasta
echo '================================'

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

echo '================================' 
echo "Calculate structuredness for all files"
echo '================================'
# Adjust for local resources
( 
  mkdir -p processing && cd processing 
  find ../rna_structures_all/ -name '*.cif' |  parallel --memsuspend 30G --resume --bar -j 60 --joblog job.log "ls {/.}_structuredness.csv > /dev/null 2> /dev/null || python ${SCRIPT_DIR}/processing/structuredness.py {} {/.}_structuredness.csv"
)
echo '================================' 
echo "Add structuredness metrics to metadata"
echo '================================'
python ${SCRIPT_DIR}/merge_metadata.py metadata_merged_clustering_with_rna3ddb_with_rna3dhub.csv 'processing/*.csv' jrc_${VERSION}_metadata_with_structuredness.csv --join-columns pdb_id chain_id
python ${SCRIPT_DIR}/postprocess/001_add_length_adjust_structuredness.py jrc_${VERSION}_metadata_with_structuredness.csv jrc_${VERSION}_metadata.csv
echo '================================'
# 
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

echo '================================'
echo "Select single representative per sequence and per group from metadata"
echo '================================'
python "${SCRIPT_DIR}/get_single_target_per_group.py" jrc_${VERSION}_metadata.csv jrc_${VERSION}_single_per_sequence.csv --group-key seq_group_id --sort_by temporal_cutoff target_id --rename target_id all_pdb_ids --rename seq_group_id target_id --merge target_id description
python "${SCRIPT_DIR}/get_single_target_per_group.py" jrc_${VERSION}_metadata_filtered.csv jrc_${VERSION}_filtered_single_per_sequence.csv --group-key seq_group_id --sort_by temporal_cutoff target_id --rename target_id all_pdb_ids --rename seq_group_id target_id --merge target_id description
python "${SCRIPT_DIR}/get_single_target_per_group.py" jrc_${VERSION}_metadata.csv jrc_${VERSION}_single_per_group.csv --sort_by temporal_cutoff target_id --rename target_id all_pdb_ids --rename group_id target_id --merge target_id description 
python "${SCRIPT_DIR}/get_single_target_per_group.py" jrc_${VERSION}_metadata_filtered.csv jrc_${VERSION}_filtered_single_per_group.csv --sort_by temporal_cutoff target_id --rename target_id all_pdb_ids --rename group_id target_id --merge target_id description 

echo '================================'
echo "Extract coordinates for all filtered RNA structures"
echo '================================'
python ${SCRIPT_DIR}/processing/extract_coordinates.py jrc_${VERSION}_metadata_filtered.csv rna_structures_all jrc_${VERSION}_filtered_coords.csv

echo '================================'
echo "Split filtered metadata and coordinates into train/test by temporal cutoff at ${SPLIT_DATE}, grouped by 30% identity clusters"
echo '================================'

python "${SCRIPT_DIR}/split_by_cutoff.py" jrc_${VERSION}_metadata.csv jrc_${VERSION}_metadata_filtered.csv --cutoff_date ${SPLIT_DATE} --group_key "mmseqs_0.300"

for part in train test; do
duckdb -c "copy (
    with coords as (
      select * from 'jrc_${VERSION}_filtered_coords.csv'
    ) 
    select * from coords where target_id in (
      select target_id from '${part}_jrc_${VERSION}_metadata_filtered.csv'
    ) 
) to '${part}_jrc_${VERSION}_coordinates.csv';"
done 
echo '================================'

python ${SCRIPT_DIR}/add_group_index.py test_jrc_${VERSION}_coordinates.csv test_jrc_${VERSION}_coordinates_gid.csv --reference test_jrc_${VERSION}_metadata_filtered.csv --group_index_name solution_id --group_key seq_group_id --sort_by_after seq_group_id solution_id ID
sed -E "1s/\\<([xyz])_1\\>/\"C1'_\\1\"/g" test_jrc_${VERSION}_coordinates_gid.csv > test_jrc_${VERSION}_coordinates_gid_c1p.csv 
python "${SCRIPT_DIR}/convert_to_multisolution_wide.py" "test_jrc_${VERSION}_coordinates_gid.csv" "test_jrc_${VERSION}_C1prime.csv" --select "C1'" --add-column Usage Public --drop-column group_id --fill-coordinates="-1e18" --pad 40