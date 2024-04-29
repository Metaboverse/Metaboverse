#!/bin/bash

#VERSION=0.11.2
#BUILD_PATH="../build"
#BUILD_EXE="${BUILD_PATH}/metaboverse-cli-nix"
#BD_DEST=j-berg@frs.sourceforge.net:/home/frs/project/metaboverse/v${VERSION}

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
SBML_URL="https://reactome.org/download/current/all_species.3.1.sbml.tgz"
SBML_FILE="all_species.sbml.tgz"
SBML_DIR="sbml"

# Get species list from Reactome
SPECIES=($(curl -s $REACTOME_API | jq -r '.[].abbreviation'))

# Download SBML archive 
curl -L $SBML_URL -o $SBML_FILE

# Extract SBML files
if [ -d "$SBML_DIR" ]; then
  rm -rf $SBML_DIR
fi
mkdir $SBML_DIR
tar xzf $SBML_FILE -C $SBML_DIR

# Build list of organisms with SBML files
declare -A species_seen=()
# Loop through SBML files
for sbml_file in $SBML_DIR/*.sbml; do

  # Get abbreviation from file name 
  IFS='-' read -ra ADDR <<< "$(basename $sbml_file .sbml)"
  abbv=${ADDR[1]}
  
  # Check if in species list
  if [[ " ${SPECIES[@]} " =~ " $abbv " ]]; then
    species_seen[$abbv]=1
  fi
done
PROCESSED_SPECIES=("${!species_seen[@]}")


printf "Processing database curation for:\n"
for X in ${PROCESSED_SPECIES[@]}; do
  echo -e "${BUILD_PATH}/${X}";
  mkdir -p ${BUILD_PATH}/${X} ;
done
rm -rf $SBML_DIR
rm -rf $SBML_FILE


# Run
printf "+ Processing database curation...\n"
parallel -j 4 ${BUILD_EXE} curate --force_new_curation --output ${BUILD_PATH}/{} --organism_id {} ::: ${PROCESSED_SPECIES[@]}


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
INCLUDE_PATTERN=()
for X in ${SPECIES[@]}; do
  if [ -f "${BUILD_PATH}/${X}/${X}_template.mvrs" ]; then
    printf "\n ${X}" >> ${BUILD_PATH}/README.txt
    INCLUDE_PATTERN+=("--include=${BUILD_PATH}/${X}/")
    INCLUDE_PATTERN+=("--include=${BUILD_PATH}/${X}/***")
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
echo -e "\nUploading to host..."

# Include other necessary files
INCLUDE_PATTERN+=("--include=build_env.txt")
INCLUDE_PATTERN+=("--include=README.txt")
INCLUDE_PATTERN+=("--include=metaboverse-cli-nix")

# Complete rsync command
RSYNC_COMMAND=("rsync" "-avzv" "-e" "ssh")
RSYNC_COMMAND+=("${INCLUDE_PATTERN[@]}")
#RSYNC_COMMAND+=("--exclude=*") # Exclude other directories and files
RSYNC_COMMAND+=("${BUILD_PATH}/" "${BD_DEST}")

# Execute rsync command
eval "${RSYNC_COMMAND[@]}"

conda deactivate