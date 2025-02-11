# nifh_amplicons_analysis

____

***nifh_amplicons_analysis*** contains the bioinformatic workflow of the analysis of the [*nifH* database](https://figshare.com/articles/dataset/Global_biogeography_of_N_sub_2_sub_-fixing_microbes_i_i_i_nifH_i_amplicon_database_and_analytics_workflow/23795943/1?file=46033371), most of which is described in [(Morando, Magasin) et al. 2025](https://essd.copernicus.org/articles/17/393/2025/essd-17-393-2025.html). This workflow is managed by a Snakefile that creates a conda environment and then executes multiple R scripts that load and process the data before preforming the ecological and biogeochemical analysis of nifH ASVs and their metadata.

> **Note**: This project is actively maintained and will expand over time.

## Table of Contents

- [nifh\_amplicons\_analysis](#nifh_amplicons_analysis)
  - [Table of Contents](#table-of-contents)
  - [Data](#data)
  - [Workflows that generated database used in analysis](#workflows-that-generated-database-used-in-analysis)
  - [Installation](#installation)
  - [Running the analysis](#running-the-analysis)
  - [Output](#output)
  - [Troubleshooting and error handling](#troubleshooting-and-error-handling)

Outputs are figures, tables, and data (csv) and log (txt) files that represent both work found in [(Morando, Magasin) et al. 2025](https://essd.copernicus.org/articles/17/393/2025/essd-17-393-2025.html) as well as other information providing a broader view of the [*nifH* database](https://figshare.com/articles/dataset/Global_biogeography_of_N_sub_2_sub_-fixing_microbes_i_i_i_nifH_i_amplicon_database_and_analytics_workflow/23795943/1?file=46033371).

## Data

____

The [*nifH* database](https://figshare.com/articles/dataset/Global_biogeography_of_N_sub_2_sub_-fixing_microbes_i_i_i_nifH_i_amplicon_database_and_analytics_workflow/23795943/1?file=46033371) comprises:

- Nearly all published _nifH_ amplicon MiSeq data sets that existed at the time of publication
- Two new data sets produced by the [Zehr Lab](https://www.jzehrlab.com/) at [UC Santa Cruz](https://www.ucsc.edu/). 

Click on the map for an interactive Google Map with detailed information (study names, sample IDs, collection information)

[![Map of studies used in Morando, Magasin et al. 2024](images_for_readme/Morando_Magasin_et_al_2024_studies_used.png)](https://www.google.com/maps/d/u/0/edit?mid=1OlWftvxU_o7Fy3nFsSJDcUlbEWSX_U0&usp=sharing)

## Workflows that generated database used in analysis

____

Two separate bioinformatic workflows were used to generate the data (the *nifH* database) used in this analysis, i.e., the ASV database and associated metadata.

1. [*nifH-amplicons-DADA2*](https://github.com/jdmagasin/nifH_amplicons_DADA2): nifH DADA2 pipeline that aggregated and produced the intial ASV table
2. [*nifH-ASV-workflow*](https://github.com/jdmagasin/nifH-ASV-workflow): Post-processing pipeline for quality filtering, data validation, and metadata acquisition via CMAP. **Ultimately generated the *nifH* database**

**Workflow overview**:
![Overview of DADA2 niH workflow](images_for_readme/workflow_overview.png)

**Workflow Steps**:
1. [DADA2](https://benjjneb.github.io/dada2/) ASVs were created by our [DADA2 _nifH_ pipeline](https://github.com/jdmagasin/nifH_amplicons_DADA2) (green).
2. Post-pipeline stages (lavender) executed by Makefile or Snakefile:
   - Gather ASVs from all studies
   - Filter ASVs for quality
   - Annotate ASVs
   - Download sample-colocated environmental data from [Simons Collaborative Marine Atlas Project (CMAP)](https://simonscmap.com)

The resulting _nifH_ ASV database supports future research into N<sub>2</sub>-fixing marine microbes.

> **Access the Database**:
> - [WorkspaceStartup directory](https://github.com/jdmagasin/nifH-ASV-workflow/tree/master/WorkspaceStartup):
>   - `nifH_ASV_database.tgz`
>   - R image `workspace.RData`
> - [Figshare](https://doi.org/10.6084/m9.figshare.23795943.v1)

## Installation

____

1. Clone the *nifh_amplicons_analysis* workflow repository:

```bash
git clone https://github.com/mo-morando/nifh_amplicons_analysis
cd nifh_amplicons_analysis
```

2. Install Snakemake (version 8.27.1 recommended):

The analysis can then be carried out by executing a Snakefile located in the scripts directory. This requires the installation of Snakemake. We recommend using a package manager, e.g., conda/mamba, and creating a contained environment. This ensures that the analysis runs smoothly and avoids potential conflicts with other programs and packages. 
> **Note**: This analysis workflow has been tested with Snakemake version *8.27.1* and so we suggest using this, however, newer versions may work.

To create a conda environment named *snakemake*:

```bash
conda create --name snakemake snakemake=8.27.1
conda activate snakemake
```

However, if you are already using the [*nifH-ASV-workflow*](https://github.com/jdmagasin/nifH-ASV-workflow), it contains a conda environment with Snakemake already, so the above step is not necessary. You can simply active this existing environment.

```bash
conda activate nifH_ASV_workflow
```

## Running the analysis

____

Once this environment is created and activated, the entire analysis is managed by a Snakemake that:
1. Creates a conda environment (`nifh_amplicons_analysis`)
2. Executes 8 R scripts for data analysis and the production of figures and tables
3. Closes the conda environment upon completion

Enter into the `scripts` directory and execute the Snakefile:

```bash
cd scripts
snakemake -c1 --use-conda
```

## Output

The workflow generates:
- Figures and tables from [(Morando, Magasin) et al. 2025](https://essd.copernicus.org/articles/17/393/2025/essd-17-393-2025.html)
- Additional figures and tables
- Log files for each process executed with detailed information

> **Note:** These log files contain potential error messages as well as detailed information regarding the processing of the data and results of interest. 
>
> In particular, *basic_sample_stats.log* and *samp_n_tax_breakdown.log* have a lot of useful information regarding stats on the samples in general, e.g., number of DNA or replicate samples, as well as breakdowns of their taxonomy supplied within its log file.

Output locations:

```bash

analysis/out_files/ # The majority of the tables (csv)
analysis/out_files/logs # All log files
analysis/out_files/plots # All figures
analysis/out_files/tables # Important tables specific to Morando, Magasin et al. 2025
```

## Troubleshooting and error handling

____

The workflow is designed to be self-contained with comprehensive error handling. For custom analyses or data changes, refer to the documentation and error messages in the log files for debugging assistance.
