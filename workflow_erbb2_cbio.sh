# Description: Workflow for regenerating ERBB2 cBioPortal files.
# Author: Haley Hunter-Zinck
# Date: 2022-05-03

Rscript uncode_redcap_export.R -e syn30042048 -d syn29880018 -o syn30041935 -f genie_erbb2_uncoded_export.csv -v
Rscript create_clinical.R -i syn29990375 -m syn29989712 -o syn30041961 -v 

