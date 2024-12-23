![MicroFisher workflow](https://github.com/NFREC-Liao-Lab/MicroFisher/actions/workflows/python-app.yml/badge.svg)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)



![microfisher_2](https://github.com/user-attachments/assets/81a83a18-654b-4e56-a533-d6ab2f97e8ef)


# MicroFisher
Profiling the taxonomic and functional composition of microbes using metagenomic and metatranscriptomic sequencing is advancing our understanding of microbial functions. However, the sensitivity and accuracy of microbial classification using genome- or core protein-based approaches, especially the classification of eukaryotic organisms, is limited by the availability of reference genomes and the resolution of sequence databases. To address this, we propose the MicroFisher, a novel tool to filter out taxonomically useful reads from metagenomic or metatranscriptomic data, enabling taxonomic identification of community members based on multiple marker regions. We applied MicroFisher to profile the simulated mock fungal communities to assess the performance of the developed tool, and found high performance in fungal prediction and abundance estimation. In addition, we also used metagenomes from forest soil and metatranscriptomics of root eukaryotic microbes to test our method and found that MicroFisher provided more accurate profiling of environmental microbiomes compared to other classification tools. Overall, MicroFisher serves as a novel pipeline for the classification of fungal communities from metagenomes and metatranscriptomics.

## Introduction
We developed the MicroFisher package, a comprehensive bioinformatics tool for analyzing fungal community composition from metagenomic and metatranscriptomic sequencing datasets using multiple hypervariable marker databases. Our approach integrates 104,072 fungi from 6,545 genera, utilizing carefully curated hypervariable regions from ITS1, ITS2, LSU D1, and LSU D2 markers with sequence lengths ranging from 120 to 350 bp. The pipeline employs a weight-based integration algorithm that maps metagenomic and metatranscriptomic reads across these hypervariable marker databases using Centrifuge, strategically addressing key challenges in fungal taxonomic classification from short-read sequencing datasets. This innovative approach provides researchers with a robust, multiple hypervariable marker strategy for comprehensive and precise fungal taxonomic analysis, overcoming the limitations of traditional classification methods.
### Database constructions
The MicroFisher pipeline compiled databases of hypervariable regions within fungal ITS and LSU sequences, covering 104,072 fungi within 6,545 genera (Fig. 1A). The hypervariable marker databases provide high-resolution (by applying hypervariable sequence regions) and multiple set references, which reduces the number of mismatches (false positives) and false negatives in fungal classification. The trimming steps (Steps 6-7 in Fig. 1A) based upon sequence consensus allow us to keep only the hypervariable region for marker databases, which improves the accuracy of fungal taxa predictions. The flanking regions of the trimmed sequences in marker databases are less than 10 bp to minimize the false positives in the fungal classification. The marker sequences have a wide range of lengths (120 - 350 bp) to enhance variability.

### MicroFisher pipeline and algorithms
MicroFisher employs multiple hypervariable marker databases (including ITS1, ITS2, LSU D1, and LSU D2) to analyze the taxonomic composition of fungal communities from metagenomics and metatranscriptomic sequence data to classify fungal composition and estimate taxa abundance from metagenomic and metatranscriptomic sequencing (Figs. 1B and 1C). The first step of the MicroFisher pipeline is to map reads to several marker databases HMDs and generate corresponding abundance reports using Centrifuge. A weight-based algorithm is further employed to optimize and integrate those Centrifuge reports into a final abundance table of detected taxa. The weight-based integration algorithm takes into account the total number of mapped reads, MiniHit length, and average sequence length of the mapped data. Thus, the weight-based integration algorithm effectively reduces mismatches (false positives), and the combination of multiple marker databases HMDs reduces the chance of missing taxa (false negatives). 

![image](https://github.com/user-attachments/assets/b065c2c4-2a7f-48f7-bc66-47cb6cefdc32)
Fig.1 Overview of hypervariable marker database generation and the principle of MicroFisher classification.



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
1. Initialize the database. The prebuild database is available at the `$MicroFisher/microfisher/default_db` folder.
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
### Initialize MicroFisher database
To get started with MicroFisher, you will have to get the database files on your system, database files are available [online](https://figshare.com/articles/dataset/MicroFisher_DBs/19679595), and you can download prebuilt database files using the code below.

Fetch the prebuilt database and store it in a custom folder `$PATH_TO_NEW_DATABASE_FOLDER`
```bash
MicroFisher init_db --help
MicroFisher init_db --db_loc $PATH_TO_DATABASE_FOLDER
```
#### Arguments
`--db_loc`: Custom location/folder for the prebuilt MicroFisher databases.



### Running MicroFisher for fungal taxonomic classification from Metagenomic and Metatranscriptomic sequencing datasets
The MicroFisher pipeline identifies the fungal taxa from Metagenomic and Metatranscriptomic datasets using the ITS and LSU rRNA sequences. Thus, the datasets that filtered the adaptors and low-quality reads were used for fungal taxonomic classification (Note: do not remove the rRNA sequences). MicroFisher pipeline applies two steps to classify the fungal community composition: (1) Alignment: Align the qualified reads of metagenomic and metagenomic to marker databases; and (2) Integration and Optimization: Integrate abundance reports using weight-based algorithms. You can run the fungal classification in one step or two steps.

### (1) Performing the fungal community classification in one step using "--preset"
Run both steps in the MicroFisher with default configurations.
```bash
MicroFisher preset --help
```

#### Different preset databases in "--preset_db"
- `ITS+LSU`: Search against four databases: ITS1, IST2, LSU_D1, and LSU_D2.
- `ITS`: Search against four databases: ITS1 and IST2
- `LSU`: Search against four databases: LSU_D1 and LSU_D2
For Metagenomic sequencing datasets, use "ITS+LSU"; for Metatranscriptomic datasets, use "LSU".

#### Examples
- Run with the test samples
   ```bash
   cd MicroFisher/microfisher
   MicroFisher preset --preset_db ITS+LSU --verbose \
   --db_path default_db \
   --min 120 \
   --workspace example \
   --paired example_R1.fastq.gz example_R2.fastq.gz \
   --out_dir merged_result_folder \
   --out_prefix results_prefix \
   --threads 4
   ```

- Basic example: the MicroFisher will run in default settings
    ```bash
    MicroFisher preset --preset_db ITS+LSU \
    --db_path $PATH_TO_DATABASE_FOLDER \
    --workspace $PATH_TO_WORKSAPCE_FOLDER \
    --paired $PATH_TO_INPUT_FILE_FOLDER/example_R1.fastq.gz $PATH_TO_INPUT_FILE_FOLDER/example_R2.fastq.gz \
    --out_dir $PATH_TO_OUTPUT_FILE_FOLDER \
    --prefix OUTPUT_FILE_NAME    
    ```

- Full configuration: modify the parameters of MicroFisher 
    ```bash
    MicroFisher preset --preset_db ITS+LSU --verbose \
    --db_path $PATH_TO_DATABASE_FOLDER \
    --min 120 \
    --workspace $PATH_TO_WORKSAPCE_FOLDER \
    --paired $PATH_TO_INPUT_FILE_FOLDER/example_R1.fastq.gz $PATH_TO_INPUT_FILE_FOLDER/example_R2.fastq.gz \
    --out_dir $PATH_TO_OUTPUT_FILE_FOLDER \
    --out_prefix OUTPUT_FILE_NAME \
    --threads 4
    ```

 - Explanation for the full configuration
    ```bash
    MicroFisher preset --preset_db ITS+LSU --verbose \  # The selected databases used for the job (Metagenomic data: ITS+LSU; Metatranscriptomic data: LSU)
    --db_path $PATH_TO_DATABASE_FOLDER \  # Path to the DATABASE of MicroFisher
    --min 120 \  # Minimum matching length (Default: 120).
    --workspace $PATH_TO_WORKSAPCE_FOLDER \  # Path to the work folder (The folder contains input files and output the temporal searching result)
    --paired $PATH_TO_INPUT_FILE_FOLDER/example_R1.fastq.gz $PATH_TO_INPUT_FILE_FOLDER/example_R2.fastq.gz \  # Path to .fastq file(s) (using --single if the data is single end reads)
    --out_dir $PATH_TO_OUTPUT_FILE_FOLDER \  # Path to folder output the results files
    --out_prefix OUTPUT_FILE_NAME \  # Prefix of the result output files
    --threads 4 #Number of threads
    ```

- Customized the path for `centrifuge`: when the package "Centrifuge" was manually installed.
    ```bash
    # e.g.
    # FULL_PATH_TO_CENTRIFUGE="/home/user/software/centrifuge/
    MicroFisher preset --preset_db ITS+LSU \
    --db_path $PATH_TO_DATABASE_FOLDER \
    --centrifuge_path $FULL_PATH_TO_CENTRIFUGE \
    --workspace $PATH_TO_WORKSAPCE_FOLDER \
    --paired $PATH_TO_INPUT_FILE_FOLDER/example_R1.fastq.gz $PATH_TO_INPUT_FILE_FOLDER/example_R2.fastq.gz \
    --out_dir $PATH_TO_OUTPUT_FILE_FOLDER \
    --out_prefix OUTPUT_FILE_NAME
    ```

  
 
 
### (2) Performing the fungal community classification in two steps using "--search" and "--combine"
#### Step 1: Alignment: Align the qualified reads of metagenomic and metagenomic to marker databases (Search fungal taxa with centrifuge).
```bash
MicroFisher search --help
```

- Example using test files
```bash
cd MicroFisher/microfisher
MicroFisher search -v
  --db_path default_db --db LSU_D2
  --workspace example \
  --prefix example \
  --min 120 \
```
- Full configuration: modify the parameters of MicroFisher 
```bash
MicroFisher search -v
  --db_path  $PATH_TO_DATABASE_FOLDER --db LSU_D2
  --workspace $PATH_TO_WORKSAPCE_FOLDER \
  --prefix OUTPUT_FILE_NAME \
  --min 120 \
```
#### Arguments
- `--prefix`: One prefix for two paired-end files.
- `--paired`: Two paired-end files.
- `--single`: One single-end file.
- `--preset_db`: Which centrifuge database to search against.
- `--min`: Minimum matching length.


  
#### Step2: Integration and Optimization: Integrate abundance reports using weight-based algorithms (Combine reports from Results using multiple databases).
```bash
MicroFisher combine --help
```

- Example using test files
```bash
  MicroFisher combine -v \
    --workspace example \
    --combine result_*_dbITS1_report.tsv result_*_dbITS2_report.tsv result_*_dbLSU_D1_report.tsv result_*_dbLSU_D2_report.tsv \
    --mode weighted --min_overlap 3
```

- Full configuration: modify the parameters of MicroFisher 
Note: `--combine` argument list multiple result files. Users might need to change these filenames.
- `weighted` mode with custom filter
    ```bash
    MicroFisher combine \
    --combine result_*_dbLSU_D1_report.tsv result_*_dbLSU_D2_report.tsv \
    --mode weighted --filter 1e-8
    ```
    
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


- `weighted` mode with custom length
    ```bash
    MicroFisher combine \
    --combine report_1.tsv report_2.tsv report_3.tsv \
    --mode weighted --cent_length 90 100 110
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

 

 

 
 



### Test
To run the tests for MicroFisher, use the following command:
```bash
pytest
```

## Database generation
The `database_generation` contains scripts for curating customized databases for MicroFisher.


## Manuscript
The `manuscript` folder contains source codes used to perform analyses and validation for the manuscript.

Citation: Wang H., S. Wu, K. Zhang, K. Chen, R. Vilgalys, H. Liao. MicroFisher: Fungal taxonomic classification for metatranscriptomic and metagenomic data using multiple short hypervariable markers.

