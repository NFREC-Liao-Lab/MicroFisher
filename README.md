![MicroFisher workflow](https://github.com/NFREC-Liao-Lab/MicroFisher/actions/workflows/python-app.yml/badge.svg)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)



![microfisher_2](https://github.com/user-attachments/assets/81a83a18-654b-4e56-a533-d6ab2f97e8ef)


# MicroFisher
Profiling the taxonomic and functional composition of microbes using metagenomic and metatranscriptomic sequencing is advancing our understanding of microbial functions. However, the sensitivity and accuracy of microbial classification using genome- or core protein-based approaches, especially the classification of eukaryotic organisms, is limited by the availability of reference genomes and the resolution of sequence databases. To address this, we propose the MicroFisher, a novel tool to filter out taxonomically useful reads from metagenomic or metatranscriptomic data, enabling taxonomic identification of community members based on multiple marker regions. We applied MicroFisher to profile the simulated mock fungal communities to assess the performance of the developed tool, and found high performance in fungal prediction and abundance estimation. In addition, we also used metagenomes from forest soil and metatranscriptomics of root eukaryotic microbes to test our method and found that MicroFisher provided more accurate profiling of environmental microbiomes compared to other classification tools. Overall, MicroFisher serves as a novel pipeline for classification of fungal communities from metagenomes and metatranscriptomics.

## Introduction
We developed the MicroFisher package, a comprehensive bioinformatics tool for analyzing fungal community composition from metagenomic and metatranscriptomic sequencing datasets using multiple hypervariable marker databases. Our approach integrates 104,072 fungi from 6,545 genera, utilizing carefully curated hypervariable regions from ITS1, ITS2, LSU D1, and LSU D2 markers with sequence lengths ranging from 120 to 350 bp. The pipeline employs a weight-based integration algorithm that maps metagenomic and metatranscriptomic reads across these hypervariable marker databases using Centrifuge, strategically addressing key challenges in fungal taxonomic classification from short-read sequencing datasets. This innovative approach provides researchers with a robust, multiple hypervariable marker strategy for comprehensive and precise fungal taxonomic analysis, overcoming limitations of traditional classification methods.

## Requirement
**Python:** https://www.python.org/

**Centrifuge:** Classifier for metagenomic sequences
- GitHub: https://github.com/DaehwanKimLab/centrifuge
- Website: https://ccb.jhu.edu/software/centrifuge/

## Installation
To use the MicroFisher, installation is required. Recommended installing MicroFisher using pip - package installer for Python [pip.pypa.io](https://pip.pypa.io).
```bash
# Create a new conda environment
# Note: The `cgi` package has been removed in Python3.13 https://peps.python.org/pep-0594/
conda create -n MicroFisher Python=3.12

# Activate MicroFisher environment
conda activate MicroFisher
conda install bioconda::centrifuge

# Clone this repository
git clone https://github.com/NFREC-Liao-Lab/MicroFisher.git
cd MicroFisher/microfisher

# Install dependencies and MicroFisher
pip3 install -r requirements.txt
pip3 install .
```


## Quick start
1. Install MicroFisher. See [Installation](#installation) section
2. Go to the `MicroFisher/microfisher` folder
    ```bash
    cd MicroFisher/microfisher
    ```
1. Initialize the database. The prebuild database is available at the `default_db` folder.
    ```bash
    MicroFisher init_db
    ```
1. Run example dataset
    ```bash
    MicroFisher preset --db_path default_db --paired example/example_R1.fastq.gz example/example_R2.fastq.gz
    # OR
    MicroFisher preset --db_path default_db --workspace example --paired example_R1.fastq.gz example_R2.fastq.gz
    ```

## Usage
### Initialize centrifuge database
To get started with MicroFisher, you will have to get the database files on your system, you can download prebuilt database files, database files are available [online](https://figshare.com/articles/dataset/MicroFisher_DBs/19679595). 
`MicroFisher init_db --help`
#### Arguments
`--db_loc`: Custom location/folder for the prebuilt centrifuge databases.

#### Examples
Fetch the prebuilt database and store it in a custom folder `$PATH_TO_NEW_DATABASE_FOLDER`
```bash
MicroFisher init_db --db_loc $PATH_TO_NEW_DATABASE_FOLDER
```

### Preset pipeline
Run both steps in the MicroFisher with default configurations.
```bash
MicroFisher preset --help
```

#### Different preset databases
- `IST+LSU`: Search against four databases: ITS1, IST2, LSU_D1, and LSU_D2.
- `IST`: ITS1 and IST2
- `LSU`: LSU_D1 and LSU_D2

#### Examples
- Basic example
    ```bash
    MicroFisher preset --preset_db ITS+LSU \
    --db_path default_db \
    --workspace example \
    --prefix example \
    --out_dir merged_results
    ```

- Full configuration
    ```bash
    MicroFisher preset --preset_db ITS+LSU --verbose \
    --db_path default_db \
    --min 120 \
    --workspace example \
    --paired example_R1.fastq.gz example_R2.fastq.gz \
    --out_dir merged_result_folder \
    --out_prefix results_prefix \
    --threads 4
    ```

    Explanation for the full configuration
    ```bash
    MicroFisher preset --preset_db ITS+LSU --verbose \  # The selected databases used for the job (Metagenomic data: ITS+LSU; Metatranscriptomic data: LSU)
    --db_path $PATH_TO_DATABASE \  # Path to the DATABASE of MicroFisher
    --min 120 \  # Minimum matching length (Default: 120).
    --workspace example \  # Path to the work folder (output the searching result)
    --paired example_R1.fastq.gz example_R2.fastq.gz \  # Path to fastq file(s) (using --single if the data is single end reads)
    --out_dir merged_result_folder \  # Path to folder output the results files
    --out_prefix results_prefix \  # Prefix of the result output files
    --threads 4 #Number of threads
    ```

- Customized the path for `centrifuge`
    ```bash
    # e.g.
    # FULL_PATH_TO_CENTRIFUGE="/home/user/software/centrifuge/
    MicroFisher preset --preset_db ITS+LSU \
    --db_path default_db \
    --centrifuge_path $FULL_PATH_TO_CENTRIFUGE \
    --workspace example \
    --paired example_R1.fastq.gz example_R2.fastq.gz
    --out_dir merged_results
    ```


### Search taxonomy with centrifuge
```bash
MicroFisher search --help
```

#### Arguments
- `--prefix`: One prefix for two paired-end files.
- `--paired`: Two paired-end files.
- `--single`: One single-end file.
- `--preset_db`: Which centrifuge database to search against.
- `--min`: Minimum matching length.


#### Examples
```bash
MicroFisher search -v
  --db_path default_db --db LSU_D2
  --workspace example \
  --prefix example \
  --min 120 \
```



### Combine reports from multiple databases
```bash
MicroFisher combine --help
```


#### Arguments
- `--combine`: List of results files to combine.
- `--mode`: Different combining mode.. See the section below.
- `--rank`: Output results for these taxonomy ranks. Default ranks: `family,genus,species`
- `--include_all`: Include unfiltered results.
- `--filter`: Filter out taxa with low proportion in the results. Default: 0.00001.


#### Different combing mode
- `weighted`: normalised by the total number of reads, and
    minimum length used in centrifuge (`--cent_length`), and the average
    length of the database (`--db_length`).
- `boolean`: Present or absent of the taxa (optional: `--min_overlap`).
- `raw`: sum of the number of reads.
<!-- -
`weighted_abundance_only`: normalised by the total number of reads (testing-only).
- `weighted_centlength_only`: normalised by the total number of reads and minimum length used in centrifuge (--cent_length) (testing-only).
-->


#### Examples
Note: `--combine` argument list multiple result files. Users might need to change these filenames.
- `boolean` mode.
    ```bash
    MicroFisher combine -v \
    --workspace example \
    --combine result_*_dbITS1_report.tsv result_*_dbITS2_report.tsv result_*_dbLSU_D1_report.tsv result_*_dbLSU_D2_report.tsv \
    --mode boolean --min_overlap 3
    ```

- `raw` mode with custom output folder
    ```bash
    MicroFisher combine \
    --combine result_*_dbITS1_report.tsv result_*_dbITS2_report.tsv \
    --mode raw --out_dir custom_output
    ```

- `weighted` mode with custom filter
    ```bash
    MicroFisher combine \
    --combine result_*_dbLSU_D1_report.tsv result_*_dbLSU_D2_report.tsv \
    --mode weighted --filter 1e-8
    ```

- `weighted` mode with custom length
    ```bash
    MicroFisher combine \
    --combine report_1.tsv report_2.tsv report_3.tsv \
    --mode weighted --cent_length 90 100 110
    ```


### Test
To run the tests for MicroFisher, use the following command:
```bash
pytest
```

## Database generation
The `database_generation` contains scripts for curating customised database for MicroFisher.


## Manuscript
The `manuscript` folder contains source codes used to perform analyses and validation for the manuscript.

Wang H., S. Wu, K. Zhang, K. Chen, R. Vilgalys, H. Liao. MicroFisher: Fungal taxonomic classification for metatranscriptomic and metagenomic data using multiple short hypervariable markers.

