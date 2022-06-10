# MicroFisher

## Installation
```bash
pip3 install microfisher_packages/
```
**OR**
```bash
cd microfisher_packages
pip3 install .
```

## Usage
### Initialise centrifuge database
`MicroFisher init_db --help`
#### Arguments
`--db_loc`: Custome location/folder for the prebuild centrifuge databases.

#### Examples
Fetch prebuild databases and store at default location.
```bash
MicroFisher init_db
```
Fetch prebuild database and store at any folder (`PATH_TO_NEW_DATABASE_FOLDER`)
```bash
MicroFisher init_db --db_loc PATH_TO_NEW_DATABASE_FOLDER
```

### Search taxonmy with centrifuge
`MicroFisher search --help`
`python3 -m microfisher search --help`

#### Arguments
`--prefix`: Prefix for the input file
`--min`: Minimum matching length
`--db`: Which centrifuge database to search against.
#### Examples
```
MicroFisher search -vv  
  -w PATH_TO_WORKSPACE \
  --prefix example_ --min 120 \
  --centrifuge_path PATH_TO_CENTRIFUGE \
  --db_path PATH_TO_DATABASE --db LSU_D2 \
```



### Combine reports from multiple databases
`MicroFisher combine --help`
`python3 -m microfisher combine --help`


```
# boolean mode: verbose output, and all reports are in PATH_TO_WORKSPACE folder
MicroFisher combine -v -w PATH_TO_WORKSPACE \
--combine report_1 report_2 report_3 \
--mode boolean --min_overlap 3

# raw mode with custom output folder
MicroFisher combine --combine report_1 report_2 report_3 \
--mode raw --out_dir custom_output

# weighted with custom filter
MicroFisher combine --combine report_1 report_2 report_3 \
--mode weighted --filter 1e-8

# weighted_length
MicroFisher combine --combine report_1 report_2 report_3 \
--mode weighted_length --length 90 100 110
```

#### Different combing mode
- boolean: Present or absent of the taxa (optional: --min_overlap).
- raw: sum of the number of reads.
- weighted: normalised by the total number of reads.
- weighted_length: normalised by the total number of reads and minimum length (requires --length).
