# Description: Workflow for regenerating ERBB2 cBioPortal files.
# Author: Haley Hunter-Zinck
# Date: 2022-05-03

# parameters
synid_file_export=syn30042048
synid_file_dd=syn29880018
synid_folder_out=syn30041935
filename=genie_erbb2_uncoded_export.csv
synid_file_uncoded=syn29990375
synid_table_map=syn29989712
synid_folder_cbio_genie=syn30041961
synid_folder_mg=syn26706564
synid_folder_cbio_final=syn30156238

# uncode redcap export
Rscript uncode_redcap_export.R -e syn30042048 -d syn29880018 -o syn30041935 -f genie_erbb2_uncoded_export.csv -v

# create cbio 
Rscript create_clinical.R -i syn29990375 -m syn29989712 -o syn30041961 -v 
Rscript create_genomic.R -i syn29990375 -g syn26706564 -o syn30041961 -v 
Rscript create_case.R -i syn30041961 -o syn30041961 -v
Rscript create_panel.R -i syn30041961 -g syn26706564 -o syn30041961 -v

# validate genie files
rm -rf erbb2/
synapse get --downloadLocation erbb2/ -r syn30041961
rm erbb2/SYNAPSE_METADATA_MANIFEST.tsv erbb2/*/SYNAPSE_METADATA_MANIFEST.tsv
python $HOME/cbioportal/core/src/main/scripts/importer/validateData.py -s erbb2/ -v -n | grep ERROR

# copy, modify, and upload
synapse cp --updateExisting --destinationId syn30041919 syn30041961
Rscript add_nongenie.R -i syn30041961 -o syn30156238 -v

# regenerate case files now that non-GENIE mutation data has been added
Rscript create_case.R -i syn30156238 -o syn30156238 -v

# validate modified files
rm -rf erbb2/
synapse get --downloadLocation erbb2/ -r syn30156238
rm erbb2/SYNAPSE_METADATA_MANIFEST.tsv erbb2/*/SYNAPSE_METADATA_MANIFEST.tsv
python $HOME/cbioportal/core/src/main/scripts/importer/validateData.py -s erbb2/ -v -n | grep ERROR

# clean up local 
