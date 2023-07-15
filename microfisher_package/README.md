# MicroFisher

## Requirement
**Centrifuge:** Classifier for metagenomic sequences
GitHub: https://github.com/DaehwanKimLab/centrifuge
Website: https://ccb.jhu.edu/software/centrifuge/

**Python**

Python >= 3.9



## Installation
Recommended install MicroFisher using pip - the package installer for Python [pip.pypa.io](pip.pypa.io).
```bash
cd microfisher_packages
pip3 install -r requirements.txt
pip3 install .
```
**OR**
```bash
pip3 install -r microfisher_package/requirements.txt
pip3 install microfisher_package
```

## Usage
### Initialize centrifuge database
`MicroFisher init_db --help`
#### Arguments
`--db_loc`: Custom location/folder for the prebuilt centrifuge databases.

#### Examples
Fetch prebuilt databases and store at default location.
```bash
MicroFisher init_db
```
Fetch prebuilt database and store is custom folder (`PATH_TO_NEW_DATABASE_FOLDER`)
```bash
MicroFisher init_db --db_loc PATH_TO_NEW_DATABASE_FOLDER
```

### Preset pipeline
Run both steps in the MicroFisher with default configurations.
`MicroFisher preset --help` \

#### Different preset databases
- IST+LSU: Search against four databases: ITS1, IST2, LSU_D1, and LSU_D2.
- IST: ITS1 and IST2
- LSU: LSU_D1 and LSU_D2

#### Examples
```bash
MicroFisher preset --preset_db ITS+LSU \
-w PATH_TO_WORKSPACE \
--prefix example_reads \
--out_dir OUTPUT_DIR 
```
Full configuration
```bash
MicroFisher preset --preset_db ITS+LSU --verbose\
--workspace PATH_TO_WORKSPACE
--paired example_R1.fastq.gz example_R2.fastq.gz \ 
--out_dir merged_result_folder --out_prefix results_prefix \
--centrifuge_path PATH_TO_CENTRIFUGE --db_path DB_PATH --threads 2
```

### Search taxonomy with centrifuge
`MicroFisher search --help` \
`python3 -m microfisher search --help`

#### Arguments
- `--prefix`: One prefix for two paired-end files.
- `--paired`: Two paired-end files.
- `--single`: One single-end file.
- `--db`: Which centrifuge database to search against.
- `--min`: Minimum matching length.


#### Examples
```
MicroFisher search -v -w PATH_TO_WORKSPACE \
  --prefix example_ --min 120 \
  --centrifuge_path PATH_TO_CENTRIFUGE \
  --db_path PATH_TO_DATABASE --db LSU_D2 \
```



### Combine reports from multiple databases
`MicroFisher combine --help` \
`python3 -m microfisher combine --help`


#### Arguments
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


#### Examples

```
# boolean mode: verbose output, and all reports are in PATH_TO_WORKSPACE folder
MicroFisher combine -v -w PATH_TO_WORKSPACE \
--combine report_1 report_2 report_3 \
--mode boolean --min_overlap 3

# raw mode with custom output folder
MicroFisher combine --combine report_1 report_2 report_3 \
--mode raw --out_dir custom_output

# weighted by abundance only and with custom filter
MicroFisher combine --combine report_1 report_2 report_3 \
--mode weighted_abundance_only --filter 1e-8

# weighted_length
MicroFisher combine --combine report_1 report_2 report_3 \
--mode weighted_length --length 90 100 110
```
