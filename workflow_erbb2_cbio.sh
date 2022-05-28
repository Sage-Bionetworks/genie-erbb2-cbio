# Description: Workflow for regenerating ERBB2 cBioPortal files.
# Author: Haley Hunter-Zinck
# Date: 2022-05-03

#########################
#       PARAMETERS      #
#########################

# file names
file_raw=genie_erbb2_uncoded_raw.csv
file_fix=genie_erbb2_uncoded_fixed.csv
file_pre=nonGENIE_data_mutations_extended_GRCC_UHN_preann.txt
file_ann=nonGENIE_data_mutations_extended_GRCC_UHN_ann.txt

# files
synid_file_uncoded_raw=syn30257223
synid_file_uncoded_fix=syn30257779
synid_file_dd=syn29880018
synid_file_uncoded=syn29990375
synid_file_ng_sam=syn30041987
synid_file_ng_bed=syn30135792
synid_file_ng_maf=syn30041988
synid_file_ng_maf2=syn30394472
synid_file_ng_maf_pre=syn30394657
synid_file_ng_maf_ann=syn30396139

# tables
synid_table_map=syn29989712

# folders
synid_folder_update=syn30041919
synid_folder_staging=syn30041935
synid_folder_mg=syn13247707
synid_folder_cbio_staging=syn30041961
synid_folder_cbio_final=syn30156238

########################
#        UNCODE        #
########################

# uncode redcap export
Rscript uncode_redcap_export.R -e $synid_file_export -d $synid_file_dd -o $synid_folder_staging -f $file_raw -v
Rscript fix_uncoded_erbb2.R -i $synid_file_uncoded_raw -o $synid_folder_staging -f $file_fix -v

#########################
#       GENIE-ONLY      #
#########################

# create cbio 
Rscript create_clinical.R -i $synid_file_uncoded_fix -m $synid_table_map -o $synid_folder_cbio_staging -v 
Rscript create_genomic.R -i $synid_file_uncoded_fix -g $synid_folder_mg -o $synid_folder_cbio_staging -v 
Rscript create_case.R -i $synid_folder_cbio_staging -o $synid_folder_cbio_staging -v
Rscript create_panel.R -i $synid_folder_cbio_staging -g $synid_folder_mg -o $synid_folder_cbio_staging -v
Rscript create_meta.R -i $synid_folder_cbio_staging -o $synid_folder_cbio_staging -v 

# validate genie files
rm -rf erbb2/
synapse get --downloadLocation erbb2/ -r $synid_folder_cbio_staging
rm erbb2/SYNAPSE_METADATA_MANIFEST.tsv erbb2/*/SYNAPSE_METADATA_MANIFEST.tsv
python ../cbioportal/core/src/main/scripts/importer/validateData.py -s erbb2/ -v -n | grep ERROR

#########################
#       NON-GENIE      #
#########################

# annotate additional genomic data
Rscript fix_supp_maf_erbb2.R -i $synid_file_ng_maf -s $synid_file_ng_maf2 -o $synid_folder_staging -f $file_pre -v
synapse get $synid_file_ng_maf_pre 
java -jar ../genome-nexus-annotation-pipeline/annotationPipeline/target/annotationPipeline-*.jar \
    -r \
    --filename $file_pre  \
    --output-filename $file_ann \
    --isoform-override uniprot \
    --post-interval-size 1000
synapse store --parentid $synid_folder_staging $file_ann
rm $file_pre $file_ann

# copy, modify, and upload
synapse cp --updateExisting --destinationId $synid_folder_update $synid_folder_cbio_staging
Rscript add_nongenie.R -i $synid_folder_cbio_staging \
  -o $synid_folder_cbio_final \
  --synid_file_ng_sam $synid_file_ng_sam \
  --synid_file_ng_bed $synid_file_ng_bed \
  --synid_file_ng_maf $synid_file_ng_maf_ann \
  -v

# regenerate case files now that non-GENIE mutation data has been added
Rscript create_case.R -i $synid_folder_cbio_final -o $synid_folder_cbio_final -v

# validate modified files
rm -rf erbb2/
synapse get --downloadLocation erbb2/ -r syn30156238
rm erbb2/SYNAPSE_METADATA_MANIFEST.tsv erbb2/*/SYNAPSE_METADATA_MANIFEST.tsv
python ../cbioportal/core/src/main/scripts/importer/validateData.py -s erbb2/ -v -n | grep ERROR
