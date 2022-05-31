#######################################################
#classification using MicroFisher
#######################################################
#!/bin/bash

conda activate metagenome

run_mock_classiciation_evaluation_test () {
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
           cd $stat_result
         for db in ITS LSU ITS_LSU; do
           for rank in class order family genus species; do
               if [ "$rank" = "class" ]; then
                  rank_lel=3
               elif [ "$rank" = "order" ]; then
                  rank_lel=4
               elif [ "$rank" = "family" ]; then
                  rank_lel=5
               elif [ "$rank" = "genus" ]; then
                  rank_lel=6
               elif [ "$rank" = "species" ]; then
                  rank_lel=7
           cat ${db}_merged_output_taxa_${rank}.tsv |awk '{print $2}' |while read line; do
               var1=$cat
