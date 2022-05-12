# GENIE ERBB2 Sponsored Project cBioPortal File Update

## Overview
The GENIE ERBB2 Sponsored Project required incorporation of data updates into the cBioPortal release files.  This repository takes the supplied GENIE and non-GENIE clinical and genomic data and converts these datasets into cBioPortal file format.  

## Installation

Ensure the base requirements are installed on your system:
- [R](https://www.r-project.org/)
- [Java](https://www.java.com/en/download/help/download_options.html)
- [maven](https://maven.apache.org/install.html)

Clone and install the cBioPortal package:
```
git clone git@github.com:cBioPortal/cbioportal.git
```

Clone Genome Nexus annotation pipeline and following instructions in the README for installation:
```
git clone git@github.com:genome-nexus/genome-nexus-annotation-pipeline.git
```

Clone this repository:
```
git clone git@github.com:Sage-Bionetworks/genie-erbb2-cbio.git
```

Install all required R packages:
```
cd genie-erbb2-cbio/
R -e 'renv::restore()'
```

## Synapse credentials

Cache your Synapse personal access token (PAT) as an environmental variable:
```
export SYNAPSE_AUTH_TOKEN={your_personal_access_token_here}
```

or store in `~/.synapseConfig` with the following format:
```
[authentication]

# either authtoken OR username and password
authtoken = {your_personal_access_token_here}
```

## Usage

Run the workflow:
```
sh workflow_erbb2_cbio.sh
```
