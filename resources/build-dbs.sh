#!/bin/bash

cd ${BUILD_PATH}

# Activate conda environment
source "${CONDA_PATH}"
conda activate pyinstaller 


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
printf "\nSTART" >> ${BUILD_PATH}/README.txt
for X in ${SPECIES[@]}; do
  if [ -f "${BUILD_PATH}/${X}/${X}.mvrs" ]; then
    printf "\n ${X}" >> ${BUILD_PATH}/README.txt
done
printf "\nEND" >> ${BUILD_PATH}/README.txt

# Print conda env 
conda list > ${BUILD_PATH}/build_env.txt


# Afterwards, upload to host
cd ${BUILD_PATH}
rsync -avz -e ssh ${BUILD_PATH}/* ${BD_DEST}/v$VERSION

#scp -r */*.mvdb ${BD_DEST}/v$VERSION/mvdb
#scp -r */*_template.mvrs ${BD_DEST}/v$VERSION/mvrs
#scp -r */*.nbdb ${BD_DEST}/v$VERSION/nbdb
#scp ${BUILD_PATH}/README.txt ${BD_DEST}/v${VERSION}
#scp ${BUILD_PATH}/build_versions.txt ${BD_DEST}/v${VERSION}