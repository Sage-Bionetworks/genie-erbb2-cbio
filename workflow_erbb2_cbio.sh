# Description: Workflow for regenerating ERBB2 cBioPortal files.
# Author: Haley Hunter-Zinck
# Date: 2022-05-03

# file names
file_raw=genie_erbb2_uncoded_raw.csv
file_fix=genie_erbb2_uncoded_fixed.csv

# files
synid_file_uncoded_raw=syn30257223
synid_file_uncoded_fix=syn30257779
synid_file_dd=syn29880018
synid_file_uncoded=syn29990375
synid_file_ng_sam=syn30041987
synid_file_ng_bed=syn30135792
synid_file_ng_maf=syn30041988

# tables
synid_table_map=syn29989712

# folders
synid_folder_update=syn30041919
synid_folder_staging=syn30041935
synid_folder_mg=syn13247707
synid_folder_cbio_staging=syn30041961
synid_folder_cbio_final=syn30156238

# uncode redcap export
Rscript uncode_redcap_export.R -e $synid_file_export -d $synid_file_dd -o $synid_folder_staging -f $file_raw -v
Rscript fix_uncoded_erbb2.R -i $synid_file_uncoded_raw -o $synid_folder_staging -f $file_fix -v

# create cbio 
Rscript create_clinical.R -i $synid_file_uncoded_fix -m $synid_table_map -o $synid_folder_cbio_staging -v 
Rscript create_genomic.R -i $synid_file_uncoded_fix -g $synid_folder_mg -o $synid_folder_cbio_staging -v 
Rscript create_case.R -i $synid_folder_cbio_staging -o $synid_folder_cbio_staging -v
Rscript create_panel.R -i $synid_folder_cbio_staging -g $synid_folder_mg -o $synid_folder_cbio_staging -v

# validate genie files
rm -rf erbb2/
synapse get --downloadLocation erbb2/ -r $synid_folder_cbio_staging
rm erbb2/SYNAPSE_METADATA_MANIFEST.tsv erbb2/*/SYNAPSE_METADATA_MANIFEST.tsv
python $HOME/cbioportal/core/src/main/scripts/importer/validateData.py -s erbb2/ -v -n | grep ERROR

# copy, modify, and upload
synapse cp --updateExisting --destinationId $synid_folder_update $synid_folder_cbio_staging
Rscript add_nongenie.R -i $synid_folder_cbio_staging -o $synid_folder_cbio_final -v

# regenerate case files now that non-GENIE mutation data has been added
Rscript create_case.R -i $synid_folder_cbio_final -o $synid_folder_cbio_final -v

# validate modified files
rm -rf erbb2/
synapse get --downloadLocation erbb2/ -r syn30156238
rm erbb2/SYNAPSE_METADATA_MANIFEST.tsv erbb2/*/SYNAPSE_METADATA_MANIFEST.tsv
python $HOME/cbioportal/core/src/main/scripts/importer/validateData.py -s erbb2/ -v -n | grep ERROR

