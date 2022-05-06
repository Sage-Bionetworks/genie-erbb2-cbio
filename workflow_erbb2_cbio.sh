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
synid_folder_cbio=syn30041961
synid_folder_mg=syn26706564

# uncode redcap export
Rscript uncode_redcap_export.R -e syn30042048 -d syn29880018 -o syn30041935 -f genie_erbb2_uncoded_export.csv -v

# create cbio 
Rscript create_clinical.R -i syn29990375 -m syn29989712 -o syn30041961 -v 
Rscript create_genomic.R -i syn29990375 -g syn26706564 -o syn30041961 -v 
Rscript create_case.R -i syn30041961 -o syn30041961 -v
Rscript create_panel.R -i syn30041961 -g syn26706564 -o syn30041961 -v

# validate using cBioPortal validator
synapse get --downloadLocation erbb2/ -r syn30041961
rm erbb2/SYNAPSE_METADATA_MANIFEST.tsv erbb2/*/SYNAPSE_METADATA_MANIFEST.tsv
python $HOME/cbioportal/core/src/main/scripts/importer/validateData.py -s erbb2/ -v -n

# clean up local 
