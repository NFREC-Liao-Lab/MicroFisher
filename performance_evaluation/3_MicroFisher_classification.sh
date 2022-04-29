
#######################################################
#classification using MicroFisher
#######################################################
#!/bin/bash

conda activate metagenome

run_mock_classification_test () {
# set the work directory
Test_Fungi_RefSeq=/home/microbiome/data_storage/SATA2/Fisher_test/Test_Fungi_RefSeq
DB_taxonomy=/home/microbiome/data_storage/SATA2/Fisher_test/short_DBs/DB_taxonomy
MicroFisher_DBs=/home/microbiome/data_storage/SATA2/Fisher_test/short_DBs/MicroFisher_DBs

mkdir $Test_Fungi_RefSeq/mock_community
mkdir $Test_Fungi_RefSeq/mock_community/classification_result

#set the minimum hit legth
for length in 70 80 90 100 110 120 130 140 150; do
   #set the number of fungal species in the simulating dataset
   for num in 50 100 200; do
       #set the number of replicates
       for replicate in 1 2 3 4 5; do      
           result_dir=$Test_Fungi_RefSeq/mock_community/simulating_${num}species_${replicate}
           seq_dir=$Test_Fungi_RefSeq/mock_community/simulating_${num}species_${replicate}/splite_seq
           stat_result=$Test_Fungi_RefSeq/mock_community/classification_result/stat_result_simulating_${num}species_${length}_${replicate}
           
           #fungal classification applying Fisher database
            mkdir $stat_result
            for DB_cat in  ITS1 ITS2 LsuD1 LsuD2 ; do
                   if [ "$DB_cat" = "ITS1" ]; then
                       DB=ITS1_fisher
                     elif [ "$DB_cat" = "ITS2" ]; then
                       DB=ITS2_fisher
                     elif [ "$DB_cat" = "LsuD1" ]; then
                       DB=LSU_D1_fisher_new
                     elif [ "$DB_cat" = "LsuD2" ]; then
                       DB=LSU_D1_fisher_new
                    fi
                   echo $DB
                   #generate the classification result by using centrifuge
                   if [ ! -f "$stat_result/simulating_${num}species_${length}_${replicate}_${DB_cat}.short_read.reprot.tsv" ]; then
                        python /home/microbiome/data_storage/SATA2/Fisher_test/MicroFisher/MicroFisher-Fungal-Profiling/python/MicroFisher.py  search -vv  \
                               -w $result_dir  \
                               --db_path $MicroFisher_DBs \
                               --prefix simulating_${num}species_${replicate}.short_read \
                               --min $length --db $DB
                    fi
               done
             
               #combine the classification results
               python /home/microbiome/data_storage/SATA2/Fisher_test/MicroFisher/MicroFisher-Fungal-Profiling/python/MicroFisher.py  combine -vv \
                   -w /home/microbiome/data_storage/SATA2/Fisher_test/Test_Fungi_RefSeq/NovaSeq_test/hitlength_test/simulating_100species_5/ 
                   --combine result_simulating_100species_5.short_read_min120_dbITS1_report.tsv result_simulating_100species_5.short_read_min120_dbITS2_fisher_report.tsv result_simulating_100species_5.short_read_min120_dbLSU_D1_fisher_new_report.tsv result_simulating_100species_5.short_read_min120_dbLSU_D2_fisher_new_report.tsv 
                   --min_overlap 3 
                   --output combined_result.report.tsv 
                   --combine_db ITS1 ITS2
                        
                        
                        
                        
                        

