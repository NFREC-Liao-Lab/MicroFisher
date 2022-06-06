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

```
python3 -m microfisher search ---help
```


### Combine reports from multiple databases
`MicroFisher combine --help`
`python3 -m microfisher combine --help`

```
MicroFisher combine -vv -w PATH_TO_WORKSPACE \
--combine report_1 report_2 report_3 \
--min_overlap 3 --output combined_result.report.tsv
```

#### Combining methods
**TODO**
- boolean
- raw  
- weight
- wegiht_length



Testing parameters.
```bash
python3 MicroFisher.py -vv  -w /home/microbiome/data_storage/SATA2/Fisher_test/ --centrifuge_path '' --db_path short_DBs/LSU_D1D2_DBs_new/ --prefix simulating_100species_r5.short_read  --min 120 --db ITS1 --dry
```

Testing Centrifuge parameters on workstationrkstation
```bash
python /home/microbiome/data_storage/SATA2/Fisher_test/MicroFisher/MicroFisher-Fungal-Profiling/python/MicroFisher.py  search -vv  -w /home/microbiome/data_storage/SATA2/Fisher_test/Test_Fungi_RefSeq/NovaSeq_test/hitlength_test/simulating_100species_5/  --db_path /home/microbiome/data_storage/SATA2/Fisher_test/short_DBs/ITS_DBs/ --prefix simulating_100species_5.short_read --dry --min 120 --db ITS1

python /home/microbiome/data_storage/SATA2/Fisher_test/MicroFisher/MicroFisher-Fungal-Profiling/python/MicroFisher.py  search -vv  -w /home/microbiome/data_storage/SATA2/Fisher_test/Test_Fungi_RefSeq/NovaSeq_test/hitlength_test/simulating_100species_5/  --db_path /home/microbiome/data_storage/SATA2/Fisher_test/short_DBs/ITS_DBs/ --prefix simulating_100species_5.short_read --dry --min 120 --db ITS2_fisher

python /home/microbiome/data_storage/SATA2/Fisher_test/MicroFisher/MicroFisher-Fungal-Profiling/python/MicroFisher.py  search -vv  -w /home/microbiome/data_storage/SATA2/Fisher_test/Test_Fungi_RefSeq/NovaSeq_test/hitlength_test/simulating_100species_5/  --db_path /home/microbiome/data_storage/SATA2/Fisher_test/short_DBs/LSU_D1D2_DBs_new/ --prefix simulating_100species_5.short_read --dry --min 120 --db LSU_D1_fisher_new

python /home/microbiome/data_storage/SATA2/Fisher_test/MicroFisher/MicroFisher-Fungal-Profiling/python/MicroFisher.py  search -vv  -w /home/microbiome/data_storage/SATA2/Fisher_test/Test_Fungi_RefSeq/NovaSeq_test/hitlength_test/simulating_100species_5/  --db_path /home/microbiome/data_storage/SATA2/Fisher_test/short_DBs/LSU_D1D2_DBs_new/ --prefix simulating_100species_5.short_read --dry --min 120 --db LSU_D2_fisher_new
```
Testing combination on workstation
```
python /home/microbiome/data_storage/SATA2/Fisher_test/MicroFisher/MicroFisher-Fungal-Profiling/python/MicroFisher.py  combine -vv -w /home/microbiome/data_storage/SATA2/Fisher_test/Test_Fungi_RefSeq/NovaSeq_test/hitlength_test/simulating_100species_5/ --combine result_simulating_100species_5.short_read_min120_dbITS1_report.tsv result_simulating_100species_5.short_read_min120_dbITS2_fisher_report.tsv result_simulating_100species_5.short_read_min120_dbLSU_D1_fisher_new_report.tsv result_simulating_100species_5.short_read_min120_dbLSU_D2_fisher_new_report.tsv --min_overlap 3 --output combined_result.report.tsv --combine_db ITS1 ITS2
```


testing combination with taxonmy ranks

python3 run_merge_reports.py  --combine eg1.report.tsv eg2.report.tsv --mode raw

python3 run_merge_reports.py  --combine eg1.report.tsv eg2.report.tsv --mode weighted

```    
python3 /home/microbiome/data_storage/SATA2/Fisher_test/MicroFisher/MicroFisher-Fungal-Profiling/python/run_merge_reports.py  \
       --combine result_simulating_100species_5.short_read_min120_dbITS1_report.tsv result_simulating_100species_5.short_read_min120_dbITS2_fisher_report.tsv \
       --mode raw

python3 /home/microbiome/data_storage/SATA2/Fisher_test/MicroFisher/MicroFisher-Fungal-Profiling/python/run_merge_reports.py  \
       --combine result_simulating_100species_5.short_read_min120_dbITS1_report.tsv result_simulating_100species_5.short_read_min120_dbITS2_fisher_report.tsv \
       --mode weighted


python3 /home/microbiome/data_storage/SATA2/Fisher_test/MicroFisher/MicroFisher-Fungal-Profiling/python/run_merge_reports.py  \
       --combine result_simulating_100species_5.short_read_min120_dbITS1_report.tsv result_simulating_100species_5.short_read_min120_dbITS2_fisher_report.tsv  result_simulating_100species_5.short_read_min120_dbLSU_D1_fisher_new_report.tsv result_simulating_100species_5.short_read_min120_dbLSU_D2_fisher_new_report.tsv \
       --mode raw
```
