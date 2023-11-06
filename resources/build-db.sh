#!/bin/bash

OS=$(uname -s)
if [[ $OS == "Darwin" ]]; then
  set +o nomatch
fi

cd ${BUILD_PATH}


# Print version info
printf "+ Version info:\n"
${BUILD_EXE} -v


# Get species IDs from Reactome
REACTOME_API="https://reactome.org/ContentService/data/species/all"
SPECIES=( $( curl -k ${REACTOME_API} | jq -r '.[].abbreviation' ) )

printf "Processing database curation for:\n"
for X in ${SPECIES[@]};
  do mkdir -p ${BUILD_PATH}/${X} ;
done

# Run
parallel ${BUILD_EXE} curate --force_new_curation --output ${BUILD_PATH}/{} --organism_id {} ::: "${SPECIES[@]}"

# Print metadata from run 
printf "+ Processing complete...\n"
printf "Metadata for bulk Metaboverse .mvdb curation:\n" >> ${BUILD_PATH}/README.txt
printf "\nDate: " >> ${BUILD_PATH}/README.txt
date '+%Y-%m-%d %H:%M:%S' >> ${BUILD_PATH}/README.txt
printf "\nMetaboverse version: " >> ${BUILD_PATH}/README.txt
${BUILD_EXE} -v >> ${BUILD_PATH}/README.txt
printf "\nReactome version: " >> ${BUILD_PATH}/README.txt
curl -kX GET https://reactome.org/ContentService/data/database/version >> ${BUILD_PATH}/README.txt
printf "\n\nOrganisms curated:" >> ${BUILD_PATH}/README.txt

# Print species list successfully
printf "\nSTART" >> ${BUILD_PATH}/README.txt
INCLUDE_PATTERN=""
EXCLUDE_PATTERN=""
for X in ${SPECIES[@]}; do
  if [ -f "${BUILD_PATH}/${X}/${X}.mvrs" ]; then
    printf "\n ${X}" >> ${BUILD_PATH}/README.txt
    INCLUDE_PATTERN+=("${X}/")
    INCLUDE_PATTERN+=("${X}/***")
  fi
done
printf "\nEND" >> ${BUILD_PATH}/README.txt

# Print conda env 
# Try printing conda env if can, if not, skip 
printf "\n\nConda environment:" >> ${BUILD_PATH}/build_env.txt
# Activate conda environment
source "${CONDA_PATH}"
conda activate pyinstaller 
conda list > ${BUILD_PATH}/build_env.txt
# If unable to view conda environment, then print message to build_env.txt 
if [ $? -ne 0 ]; then
  printf "\n\nUnable to view conda environment." >> ${BUILD_PATH}/build_env.txt
fi


# Afterwards, upload to host
cd ${BUILD_PATH}
#chmod -R 755 ${BUILD_PATH}

# Include the specific directories and their content
INCLUDE_PATTERN=()
for X in "${SPECIES[@]}"; do
    if [ -f "${BUILD_PATH}/${X}/${X}.mvrs" ]; then
        INCLUDE_PATTERN+=("--include=${X}/")
        INCLUDE_PATTERN+=("--include=${X}/***")
    fi
done

# Include other necessary files
INCLUDE_PATTERN+=("--include=build_env.txt")
INCLUDE_PATTERN+=("--include=README.txt")
INCLUDE_PATTERN+=("--include=metaboverse-cli-nix")

# Complete rsync command
RSYNC_COMMAND=("rsync" "-avzv" "-e" "ssh")
RSYNC_COMMAND+=("${INCLUDE_PATTERN[@]}")
RSYNC_COMMAND+=("--exclude=*") # Exclude other directories and files
RSYNC_COMMAND+=("${BUILD_PATH}/" "${BD_DEST}")

# Execute rsync command
eval "${RSYNC_COMMAND[@]}"

conda deactivate