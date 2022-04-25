
#######################################################
#classification using MicroFisher
#######################################################
#!/bin/bash

conda activate metagenome

run_miniLength_test () {
# set the work directory
Test_Fungi_RefSeq=/home/microbiome/data_storage/SATA2/Fisher_test/Test_Fungi_RefSeq
DB_taxonomy=/home/microbiome/data_storage/SATA2/Fisher_test/short_DBs/DB_taxonomy
ITS_DBs=/home/microbiome/data_storage/SATA2/Fisher_test/short_DBs/ITS_DBs
LSU_D1D2_DBs=/home/microbiome/data_storage/SATA2/Fisher_test/short_DBs/LSU_D1D2_DBs
LSU_D1D2_DBs_new=/home/microbiome/data_storage/SATA2/Fisher_test/short_DBs/LSU_D1D2_DBs_new
DB_combination=/home/microbiome/data_storage/SATA2/Fisher_test/short_DBs/DB_combination
ITS1_fisher=$ITS_DBs/ITS1_fisher
ITS2_fisher=$ITS_DBs/ITS2_fisher
LsuD1_fisher=$LSU_D1D2_DBs_new/LSU_D1_fisher_new
LsuD2_fisher=$LSU_D1D2_DBs_new/LSU_D2_fisher_new
DB_combination_fisher=$DB_combination/DB_combination_fisher
mkdir $Test_Fungi_RefSeq/hitlength_test
mkdir $Test_Fungi_RefSeq/hitlength_test/stat_result
#set the minimum hit legth
for length in 70 80 90 100 110 120 130 140 150; do
   #set the number of fungal species in the simulating dataset
   for num in 100; do
       #set the number of replicates
       for replicate in 1 2 3 4 5; do      
           result_dir=$Test_Fungi_RefSeq/hitlength_test/simulating_${num}species_${replicate}
           seq_dir=$Test_Fungi_RefSeq/hitlength_test/simulating_${num}species_${replicate}/splite_seq
           stat_result=$Test_Fungi_RefSeq/hitlength_test/stat_result/stat_result_simulating_${num}species_${length}_${replicate}
           
           #fungal classification applying Fisher database
            mkdir $stat_result
            for DB_cat in  ITS1 ITS2 LsuD1 LsuD2 ; do
                   if [ "$DB_cat" = "ITS1" ]; then
                       DB=$ITS1_fisher
                     elif [ "$DB_cat" = "ITS2" ]; then
                       DB=$ITS2_fisher
                     elif [ "$DB_cat" = "LsuD1" ]; then
                       DB=$LsuD1_fisher
                     elif [ "$DB_cat" = "LsuD2" ]; then
                       DB=$LsuD2_fisher
                    fi
                   echo $DB
                   #generate the classification result by using centrifuge
                   if [ ! -f "$stat_result/simulating_${num}species_${length}_${replicate}_${DB_cat}.short_read.reprot.tsv" ]; then
                        time centrifuge -p 24 -k 1  --min-hitlen $length -x $DB  \
                            -1  $result_dir/simulating_${num}species_${replicate}.short_read_R1.fastq.gz  \
                            -2  $result_dir/simulating_${num}species_${replicate}.short_read_R2.fastq.gz  \
                            -S  $stat_result/simulating_${num}species_${length}_${replicate}_${DB_cat}.short_read.result \
                            --report-file $stat_result/simulating_${num}species_${length}_${replicate}_${DB_cat}.short_read.reprot.tsv 
                   fi
                        time centrifuge-kreport -x $DB  \
                             $stat_result/simulating_${num}species_${length}_${replicate}_${DB_cat}.short_read.reprot.tsv \
                               > $stat_result/simulating_${num}species_${length}_${replicate}_${DB_cat}.short_read.kreprot.tsv
