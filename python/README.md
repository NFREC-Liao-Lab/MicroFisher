## Work in progress

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
```

