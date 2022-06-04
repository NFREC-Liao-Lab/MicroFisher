
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
mkdir $Test_Fungi_RefSeq/mock_community/classification_result_collection
classification_result_collection=$Test_Fungi_RefSeq/mock_community/classification_result_collection

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
            yes | cp $result_dir/*.txt $stat_result
            yes | cp $result_dir/*.txt $classification_result_collection
            for DB_cat in  ITS1 ITS2 LsuD1 LsuD2 ; do
                   if [ "$DB_cat" = "ITS1" ]; then
                       DB=ITS1_fisher
                     elif [ "$DB_cat" = "ITS2" ]; then
                       DB=ITS2_fisher
                     elif [ "$DB_cat" = "LsuD1" ]; then
                       DB=LSU_D1_fisher_new
                     elif [ "$DB_cat" = "LsuD2" ]; then
                       DB=LSU_D2_fisher_new
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
                   yes | cp $result_dir/*_report.tsv  $classification_result_collection
                   yes | cp $result_dir/*report_kreport.tsv $classification_result_collection
                   mv -f $result_dir/*_report*    $stat_result
                   mv -f $result_dir/*_output*   $stat_result
             done
             
             
               #combine the classification results
              # python /home/microbiome/data_storage/SATA2/Fisher_test/MicroFisher/MicroFisher-Fungal-Profiling/python/MicroFisher.py  combine -vv \
               #    -w $stat_result \
               #    --combine result_simulating_${num}species_${replicate}.short_read_min${length}_dbITS1_fisher_report.tsv result_simulating_${num}species_${replicate}.short_read_min${length}_dbITS2_fisher_report.tsv  \
               #    --min_overlap 2 \
               #    --output $stat_result/combine_result_ITS.tsv \
               #    --combine_db ITS
                   
               # python /home/microbiome/data_storage/SATA2/Fisher_test/MicroFisher/MicroFisher-Fungal-Profiling/python/MicroFisher.py  combine -vv \
               #    -w $stat_result \
               #    --combine result_simulating_${num}species_${replicate}.short_read_min${length}_dbLSU_D1_fisher_new_report.tsv result_simulating_${num}species_${replicate}.short_read_min${length}_dbLSU_D2_fisher_new_report.tsv  \
               #    --min_overlap 2 \
                #   --output $stat_result/combine_result_LSU.tsv \
               #    --combine_db LSU
                   
             #   python /home/microbiome/data_storage/SATA2/Fisher_test/MicroFisher/MicroFisher-Fungal-Profiling/python/MicroFisher.py  combine -vv \
             #      -w $stat_result \
               #    --combine result_simulating_${num}species_${replicate}.short_read_min${length}_dbITS1_fisher_report.tsv result_simulating_${num}species_${replicate}.short_read_min${length}_dbITS2_fisher_report.tsv result_simulating_${num}species_${replicate}.short_read_min${length}_dbLSU_D1_fisher_new_report.tsv result_simulating_${num}species_${replicate}.short_read_min${length}_dbLSU_D2_fisher_new_report.tsv  \
               #    --min_overlap 2 \
                #   --output $stat_result/combine_result_ITS_LSU.tsv \
                #   --combine_db ITS,LSU
             
             ##combine with taxonomy ranks
             cd $stat_result
             python  /home/microbiome/data_storage/SATA2/Fisher_test/MicroFisher/MicroFisher-Fungal-Profiling/python/run_merge_reports.py  \
                    --combine result_simulating_${num}species_${replicate}.short_read_min${length}_dbITS1_fisher_report.tsv \
                              result_simulating_${num}species_${replicate}.short_read_min${length}_dbITS2_fisher_report.tsv \
                    --mode raw
             for merged_file in  $(ls merged_output*) ; do
                 cp $merged_file  $classification_result_collection/result_simulating_${num}species_${replicate}.short_read_min${length}_ITS_${merged_file}
                 mv $merged_file ITS_${merged_file}
             done
             
             python  /home/microbiome/data_storage/SATA2/Fisher_test/MicroFisher/MicroFisher-Fungal-Profiling/python/run_merge_reports.py  \
                    --combine result_simulating_${num}species_${replicate}.short_read_min${length}_dbLSU_D1_fisher_new_report.tsv \
                              result_simulating_${num}species_${replicate}.short_read_min${length}_dbLSU_D2_fisher_new_report.tsv \
                    --mode raw
             for merged_file in  $(ls merged_output*) ; do
                 cp $merged_file  $classification_result_collection/result_simulating_${num}species_${replicate}.short_read_min${length}_LSU_${merged_file}
                 mv $merged_file LSU_${merged_file}
             done
                   
             python  /home/microbiome/data_storage/SATA2/Fisher_test/MicroFisher/MicroFisher-Fungal-Profiling/python/run_merge_reports.py  \
                    --combine result_simulating_${num}species_${replicate}.short_read_min${length}_dbLSU_D1_fisher_new_report.tsv \
                              result_simulating_${num}species_${replicate}.short_read_min${length}_dbLSU_D2_fisher_new_report.tsv \
                              result_simulating_${num}species_${replicate}.short_read_min${length}_dbITS1_fisher_report.tsv \
                              result_simulating_${num}species_${replicate}.short_read_min${length}_dbITS2_fisher_report.tsv \
                    --mode raw
             for merged_file in  $(ls merged_output*) ; do
                 cp $merged_file  $classification_result_collection/result_simulating_${num}species_${replicate}.short_read_min${length}_ITS_LSU_${merged_file}
                 mv $merged_file ITS_LSU_${merged_file}
             done
             
        done
    done
 done
 }                    
                        
                        
                        
                        
 
 
 function print () {
    arg1=$1
    arg2=$2
    
    echo $arg1
    echo $arg2
    echo $arg1$arg2
 }
 
 
 
