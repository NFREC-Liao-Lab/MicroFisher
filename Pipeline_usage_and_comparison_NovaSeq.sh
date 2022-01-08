#The scripts for an example pipeline usage and comparisons with other packages

#The example data collection
#Download the RefSeq 28S_fungal_sequences from NCBI
wget -c https://ftp.ncbi.nih.gov/blast/db/28S_fungal_sequences.tar.gz
tar -zxvf 28S_fungal_sequences.tar.gz

#Extract sample sequences from the 28S_fungal_sequences

#######################################################
#test using the refseq rRNA sequences
#######################################################
#!/bin/bash

conda activate metagenome

run_miniLength_test () {

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

for length in 70 80 90 100 110 120 130 140 150; do
   for num in 100; do
       for replicate in 1 2 3 4 5; do      
           result_dir=$Test_Fungi_RefSeq/hitlength_test/simulating_${num}species_${replicate}
           seq_dir=$Test_Fungi_RefSeq/hitlength_test/simulating_${num}species_${replicate}/splite_seq
           stat_result=$Test_Fungi_RefSeq/hitlength_test/stat_result/stat_result_simulating_${num}species_${length}_${replicate}
           #check if the result directory exist or not
           if [ ! -d "$result_dir/" ];then
              mkdir $result_dir
              mkdir $seq_dir
              mkdir $seq_dir/ITS
              mkdir $seq_dir/LSU
           else
              echo "$result_dir   exist"
           fi
           #check the simulating dataset exist or not, if not, then simulate the from RefSeq
           if [ ! -f "$result_dir/simulating_${num}species_${replicate}.short_read_R1.fastq.gz" ];then
                 #randomly select fungal species from the RefSeq database from simulating
                 cat $Test_Fungi_RefSeq/RefSeq_Fungi/28S_ITS_common_species_unique |shuf -n $num > $result_dir/species_list
                 cat $result_dir/species_list |while read taxa; do
                     genus=$(echo $taxa |cut -d" " -f 1)
                     species=$(echo $taxa |cut -d" " -f 2)
                     cat $Test_Fungi_RefSeq/RefSeq_Fungi/ITS_RefSeq_Fungi.fa  |seqkit seq -w 0 -g -u |grep -A 1 "${taxa}" > $seq_dir/ITS/${genus}_${species}.fasta
                     cat $Test_Fungi_RefSeq/RefSeq_Fungi/28S_fungal_sequences.fa  |seqkit seq -w 0 -g -u |grep -A 1 "${taxa}" > $seq_dir/LSU/${genus}_${species}.fasta
                 done
                 #generate the taxonomy of selected fungal species
                 for fasta in $(ls $seq_dir/ITS/*.fasta) ; do
                     for acc in $(cat $fasta |seqkit seq -i -n);do
                         accession=$acc
                         taxid=$(echo $acc |esummary -db nuccore  | xtract -pattern DocumentSummary -element TaxId)
                         linage=$(echo $taxid | taxonkit lineage --data-dir $DB_taxonomy |taxonkit reformat -P --data-dir $DB_taxonomy --delimiter ";" --format "{k};{p};{c};{o};{f};{g};{s}" | cut -f 1,3)
                         taxonomy=$accession'   '$linage
                         echo $taxonomy >> $result_dir/simulating_${num}species_${replicate}.taxonomy.txt
                     done
                 done
                 
                 for fasta in $(ls $seq_dir/LSU/*.fasta) ; do
                     for acc in $(cat $fasta |seqkit seq -i -n);do
                         accession=$acc
                         taxid=$(echo $acc |esummary -db nuccore  | xtract -pattern DocumentSummary -element TaxId)
                         linage=$(echo $taxid | taxonkit lineage --data-dir $DB_taxonomy |taxonkit reformat -P --data-dir $DB_taxonomy --delimiter ";" --format "{k};{p};{c};{o};{f};{g};{s}" | cut -f 1,3)
                         taxonomy=$accession'   '$linage
                         echo $taxonomy >> $result_dir/simulating_${num}species_${replicate}.taxonomy.txt
                     done
                 done
               #replace sequences id with fungal species name
               #  for sequence in $(ls $seq_dir/ITS/*.fasta); do
               #      seq_taxa=$(echo $sequence |cut -d"/" -f 12 |cut -d "." -f 1)
               #      cat $sequence |seqkit replace -p ".+" -r "${seq_taxa}" > $seq_dir/ITS/$seq_taxa.fa
               #      rm -f $seq_dir/ITS/$seq_taxa.fasta
               #  done
                 
               #  for sequence in $(ls $seq_dir/LSU/*.fasta); do
               #      seq_taxa=$(echo $sequence |cut -d"/" -f 12 |cut -d "." -f 1)
               #      cat $sequence |seqkit replace -p ".+" -r "${seq_taxa}" > $seq_dir/LSU/$seq_taxa.fa
               #      rm -f $seq_dir/LSU/$seq_taxa.fasta
                # done
                 #gernrate the simulating dataset using iss
                 time iss generate --cpus 24 --quiet  --compress \
                       --genomes $seq_dir/ITS/*.fasta \
                       --model novaseq \
                       --abundance exponential \
                       --n_reads 5M \
                       --output $result_dir/simulating_${num}species_${replicate}_ITS.short_read
           #cat $result_dir/simulating_${num}species_${length}_${replicate}_ITS.short_read_abundance.txt |sed 's/ITS/LSU/g' > $result_dir/simulating_${num}species_${length}_${replicate}_LSU.short_read_abundance.txt
                 time iss generate --cpus 24 --quiet  --compress \
                       --genomes $seq_dir/LSU/*.fasta \
                       --model novaseq \
                       --abundance exponential \
                       --n_reads 5M \
                       --output $result_dir/simulating_${num}species_${replicate}_LSU.short_read
            #combine the simulating dataset from ITS and LSU RefSeq database
            cat $result_dir/simulating_${num}species_${replicate}_ITS.short_read_R1.fastq.gz  $result_dir/simulating_${num}species_${replicate}_LSU.short_read_R1.fastq.gz > $result_dir/simulating_${num}species_${replicate}.short_read_R1.fastq.gz
            cat $result_dir/simulating_${num}species_${replicate}_ITS.short_read_R2.fastq.gz  $result_dir/simulating_${num}species_${replicate}_LSU.short_read_R2.fastq.gz > $result_dir/simulating_${num}species_${replicate}.short_read_R2.fastq.gz
            
            rm -f $result_dir/simulating_${num}species_${replicate}_ITS.short_read_R1.fastq.gz  $result_dir/simulating_${num}species_${replicate}_LSU.short_read_R1.fastq.gz
            rm -f $result_dir/simulating_${num}species_${replicate}_ITS.short_read_R2.fastq.gz  $result_dir/simulating_${num}species_${replicate}_LSU.short_read_R2.fastq.gz
            
            else
                echo "$result_dir/simulating_${num}species_${replicate}.short_read_R1.fastq.gz  exist"
            fi
            
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
                    
                    
                    #evaluate the classification performance 
                    rm -f $stat_result/simulating_${num}species_${length}_${replicate}_${DB_cat}_Ture_Positive.txt
                    rm -f $stat_result/simulating_${num}species_${length}_${replicate}_${DB_cat}_False_Positive.txt
                    rm -f $stat_result/simulating_${num}species_${length}_${replicate}_${DB_cat}_False_negative.txt
                    echo "###################################################################################################" >> $stat_result/simulating_${num}species_${length}_${replicate}_${DB_cat}_Ture_Positive.txt
                    echo "###################################################################################################" >> $stat_result/simulating_${num}species_${length}_${replicate}_${DB_cat}_False_Positive.txt
                    echo "###################################################################################################" >> $stat_result/simulating_${num}species_${length}_${replicate}_${DB_cat}_False_negative.txt
                    for level in S G F O C P ; do
                        if [ "$level" = "S" ]; then
                           TL=7
                        elif [ "$level" = "G" ]; then
                           TL=6
                        elif [ "$level" = "F" ]; then
                           TL=5  
                        elif [ "$level" = "O" ]; then
                           TL=4  
                        elif [ "$level" = "C" ]; then
                           TL=3  
                        elif [ "$level" = "P" ]; then
                           TL=2
                        fi
                        echo $TL
         
                        echo $level"######################################" >> $stat_result/simulating_${num}species_${length}_${replicate}_${DB_cat}_Ture_Positive.txt
                        echo $level"######################################" >> $stat_result/simulating_${num}species_${length}_${replicate}_${DB_cat}_False_Positive.txt
                        echo $level"######################################" >> $stat_result/simulating_${num}species_${length}_${replicate}_${DB_cat}_False_negative.txt
                        cat $stat_result/simulating_${num}species_${length}_${replicate}_${DB_cat}.short_read.kreprot.tsv |grep -w "${level}" |awk '{ for(i=1; i<=5; i++){ $i=""} ; print $0 }' |sort -u|while read line; do
                               var1=$(cat $result_dir/simulating_${num}species_${replicate}.taxonomy.txt |cut -d ";" -f $TL | grep "$line" )
                               if [ "$var1" != "" ]; then
                                   echo $line >> $stat_result/simulating_${num}species_${length}_${replicate}_${DB_cat}_Ture_Positive.txt
                               else 
                                   echo $line >> $stat_result/simulating_${num}species_${length}_${replicate}_${DB_cat}_False_Positive.txt
                               fi
                        done
         
                        cat $result_dir/simulating_${num}species_${replicate}.taxonomy.txt |cut -d ";" -f $TL |cut -d"_" -f 3 |sort -u |while read taxa; do
                            var2=$(cat $stat_result/simulating_${num}species_${length}_${replicate}_${DB_cat}.short_read.kreprot.tsv |grep -w "${level}" |grep "$taxa")
                              if [ "$var2" = "" ]; then
                                   echo $taxa >> $stat_result/simulating_${num}species_${length}_${replicate}_${DB_cat}_False_negative.txt
                              fi
                        done
                  done
                                  
                  #generate and combine the classification performance evaluation values
                  rm -f $stat_result/simulating_${num}species_${length}_${replicate}_${DB_cat}_*.stat.csv
                  rm -f $stat_result/simulating_${num}species_${length}_${replicate}_${DB_cat}.stat_combine.csv
                  for value in False_negative False_Positive Ture_Positive; do
                          file=$stat_result/simulating_${num}species_${length}_${replicate}_${DB_cat}_${value}.txt
                          stat_file=$stat_result/simulating_${num}species_${length}_${replicate}_${DB_cat}_${value}.stat.csv
                          stat_file_combine=$stat_result/simulating_${num}species_${length}_${replicate}_${DB_cat}.stat_combine.csv
                          
                          touch $stat_file_combine
                          FN_S=$[$(cat $file |grep "S###########" -A 99999999  |grep "G#########" -B 9999999|wc -l)-2]
                          FN_G=$[$(cat $file |grep "G###########" -A 99999999  |grep "F#########" -B 9999999|wc -l)-2]
                          FN_F=$[$(cat $file |grep "F###########" -A 99999999  |grep "O#########" -B 9999999|wc -l)-2]
                          FN_O=$[$(cat $file |grep "O###########" -A 99999999  |grep "C#########" -B 9999999|wc -l)-2]
                          FN_C=$[$(cat $file |grep "C###########" -A 99999999  |grep "P#########" -B 9999999|wc -l)-2]
                          FN_P=$[$(cat $file |grep "P###########" -A 99999999  |wc -l)-1]
                          
                          echo "level" " " "${value}" >> $stat_file
                          echo "Species" " " $FN_S >> $stat_file
                          echo "Genus" " " $FN_G >> $stat_file
                          echo "Family" " " $FN_F >> $stat_file
                          echo "Order" " " $FN_O >> $stat_file
                          echo "Class" " " $FN_C >> $stat_file
                          echo "Phylum" " " $FN_P >> $stat_file
                          
                          paste $stat_file_combine $stat_file > $stat_result/new
                          mv -f $stat_result/new $stat_file_combine
                          rm -f $stat_result/new     
                  done
                     
            done
       done             
   done
done

}


export -f run_miniLength_test
 
run_miniLength_test



rm -f DBs_separated/*
echo -ne "group_id\tspeciesNum\tDatabase\tlength\ttaxa_level\t\tFalse_negative\tFalse_Positive\tTure_Positive\n" >> DBs_separated/hitlength_DBs_separated.txt
for f in stat_result_simulating* ; do
    i=$(echo $f |cut -d"_" -f 3,4,5,6,7,8)
    group_id=$(echo $f |cut -d"_" -f 3,4,5,6,7,8)
    length=$(echo $f |cut -d"_" -f 5)
    speciesNum=$(echo $f |cut -d"_" -f 4)
    
    for DB in ITS1 ITS2 LsuD1 LsuD2; do
       for taxa_level in Species Genus Family Order Class Phylum; do
          if [ "$taxa_level" = "Species" ]; then
                     TL=2
          elif [ "$taxa_level" = "Genus" ]; then
                     TL=3
          elif [ "$taxa_level" = "Family" ]; then
                     TL=4  
          elif [ "$taxa_level" = "Order" ]; then
                     TL=5  
          elif [ "$taxa_level" = "Class" ]; then
                     TL=6  
          elif [ "$taxa_level" = "Phylum" ]; then
                     TL=7
          fi
          
          value=$(cat $f/${i}_${DB}.stat_combine.csv |awk NR==$TL |sed "s/$taxa_level//g")
          echo -ne "$group_id\t$speciesNum\t$DB\t$length\t$taxa_level $value\n" >> DBs_separated/hitlength_DBs_separated.txt
       done
    done
done







for f in stat_result_simulating* ; do
    i=$(echo $f |cut -d"_" -f 3,4,5,6,7,8)
    echo ${i} >> DBs_separated/${i}.final_combine.csv
    echo "ITS1" >> DBs_separated/${i}.final_combine.csv
    cat $f/${i}_ITS1.stat_combine.csv >> DBs_separated/${i}.final_combine.csv
    echo "ITS2" >> DBs_separated/${i}.final_combine.csv
    cat $f/${i}_ITS2.stat_combine.csv >> DBs_separated/${i}.final_combine.csv  
    echo "LsuD1" >> DBs_separated/${i}.final_combine.csv
    cat $f/${i}_LsuD1.stat_combine.csv >> DBs_separated/${i}.final_combine.csv  
    echo "LsuD2" >> DBs_separated/${i}.final_combine.csv
    cat $f/${i}_LsuD2.stat_combine.csv >> DBs_separated/${i}.final_combine.csv  
done




run_miniLength_combination_test () {

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


for length in 150 140 130 120 110 100 90 80 70; do
   for num in 100; do
       for replicate in 1 2 3 4 5; do      
           result_dir=$Test_Fungi_RefSeq/hitlength_test/simulating_${num}species_${replicate}
           seq_dir=$Test_Fungi_RefSeq/hitlength_test/simulating_${num}species_${replicate}/splite_seq
           stat_result=$Test_Fungi_RefSeq/hitlength_test/stat_result/stat_result_simulating_${num}species_${length}_${replicate}
           
           ITS1_report=$stat_result/simulating_${num}species_${length}_${replicate}_ITS1.short_read.reprot.tsv
           ITS2_report=$stat_result/simulating_${num}species_${length}_${replicate}_ITS2.short_read.reprot.tsv
           LusD1_report=$stat_result/simulating_${num}species_${length}_${replicate}_LsuD1.short_read.reprot.tsv
           LusD2_report=$stat_result/simulating_${num}species_${length}_${replicate}_LsuD2.short_read.reprot.tsv
           
           cat $ITS1_report |cut -f 2 > $stat_result/ITS1_taxID
           cat $ITS2_report |cut -f 2 > $stat_result/ITS2_taxID
           cat $LusD1_report |cut -f 2 > $stat_result/LusD1_taxID
           cat $LusD2_report |cut -f 2 > $stat_result/LusD2_taxID
           
           
           
           rm -f $stat_result/ITS_LSU_database_final_prediction
           cat $stat_result/ITS1_taxID | while read line ; do
               taxID=$line
               var1=$( cat $stat_result/LusD1_taxID | grep -w "${taxID}" )
               var2=$( cat $stat_result/LusD2_taxID | grep -w "${taxID}" )
               var3=$( cat $stat_result/ITS2_taxID | grep -w "${taxID}" )
               if [ "$var1" != "" -o "$var2" != "" -o "$var3" != "" ]; then
                  echo $taxID >> $stat_result/ITS_LSU_database_final_prediction
               fi
            done
            cat $stat_result/ITS2_taxID | while read line ; do
               taxID=$line
               var4=$( cat $stat_result/LusD1_taxID | grep -w "${taxID}" )
               var5=$( cat $stat_result/LusD2_taxID | grep -w "${taxID}" )
               if [ "$var4" != "" -o "$var5" != "" ]; then
                  echo $taxID >> $stat_result/ITS_LSU_database_final_prediction
               fi
            done
            cat $stat_result/LusD1_taxID | while read line ; do
               taxID=$line
               var6=$( cat $stat_result/LusD2_taxID | grep -w "${taxID}" )
               if [ "$var6" ]; then
                  echo $taxID >> $stat_result/ITS_LSU_database_final_prediction
               fi
            done
            
            
            rm -f $stat_result/ITS_database_final_prediction
           cat $stat_result/ITS1_taxID | while read line ; do
               taxID=$line
               var1=$( cat $stat_result/ITS2_taxID | grep -w "${taxID}" )
               if [ "$var1" != "" ]; then
                  echo $taxID >> $stat_result/ITS_database_final_prediction
               fi
            done
            
            
            rm -f $stat_result/LSU_database_final_prediction
            cat $stat_result/LusD1_taxID | while read line ; do
               taxID=$line
               var1=$( cat $stat_result/LusD2_taxID | grep -w "${taxID}" )
               if [ "$var1" != "" ]; then
                  echo $taxID >> $stat_result/LSU_database_final_prediction
               fi
            done
              
            rm -f $stat_result/ITS1_taxID $stat_result/ITS2_taxID $stat_result/LusD1_taxID $stat_result/LusD2_taxID
                   
            rm -f $stat_result/ITS_LSU_database_final_prediction.report.tsv
            cat $ITS1_report |head -n 1 >> $stat_result/ITS_LSU_database_final_prediction.report.tsv
            cat $stat_result/ITS_LSU_database_final_prediction|sort -u |while read line; do
                taxID=$line
                var1=$(cat $ITS1_report |awk -F "\t"  '{if ($2 == '$taxID') print }')
                var2=$(cat $ITS2_report |awk -F "\t"  '{if ($2 == '$taxID') print }')
                var3=$(cat $LusD2_report |awk -F "\t"  '{if ($2 == '$taxID') print }')
                var4=$(cat $LusD1_report |awk -F "\t"  '{if ($2 == '$taxID') print }')
                if [ "$var1" != "" ]; then
                   cat $ITS1_report |awk -F "\t"  '{if ($2 == '$taxID') print }' >> $stat_result/ITS_LSU_database_final_prediction.report.tsv
                elif [ "$var2" != "" ]; then
                   cat $ITS2_report |awk -F "\t"  '{if ($2 == '$taxID') print }' >> $stat_result/ITS_LSU_database_final_prediction.report.tsv
                elif [ "$var3" != "" ]; then
                   cat $LusD2_report |awk -F "\t"  '{if ($2 == '$taxID') print }' >> $stat_result/ITS_LSU_database_final_prediction.report.tsv
                elif [ "$var4" != "" ]; then
                   cat $LusD1_report |awk -F "\t"  '{if ($2 == '$taxID') print }' >> $stat_result/ITS_LSU_database_final_prediction.report.tsv
                fi
            done

            
            rm -f $stat_result/ITS_database_final_prediction.report.tsv
            cat $ITS1_report |head -n 1 >> $stat_result/ITS_database_final_prediction.report.tsv
            cat $stat_result/ITS_database_final_prediction|sort -u |while read line; do
                taxID=$line
                var1=$(cat $ITS1_report |awk -F "\t"  '{if ($2 == '$taxID') print }')
                if [ "$var1" != "" ]; then
                   cat $ITS1_report |awk -F "\t"  '{if ($2 == '$taxID') print }' >> $stat_result/ITS_database_final_prediction.report.tsv
                fi
            done
            
            rm -f $stat_result/LSU_database_final_prediction.report.tsv
            cat $ITS1_report |head -n 1 >> $stat_result/LSU_database_final_prediction.report.tsv
            cat $stat_result/LSU_database_final_prediction|sort -u |while read line; do
                taxID=$line
                var1=$(cat $LusD2_report |awk -F "\t"  '{if ($2 == '$taxID') print }')
                if [ "$var1" != "" ]; then
                   cat $LusD2_report |awk -F "\t"  '{if ($2 == '$taxID') print }' >> $stat_result/LSU_database_final_prediction.report.tsv
                fi
            done
            
            
            time centrifuge-kreport -x $DB_combination_fisher  \
                          $stat_result/ITS_LSU_database_final_prediction.report.tsv \
                            > $stat_result/ITS_LSU_database_final_prediction.kreprot.tsv            
            
            time centrifuge-kreport -x $DB_combination_fisher  \
                          $stat_result/ITS_database_final_prediction.report.tsv \
                            > $stat_result/ITS_database_final_prediction.kreprot.tsv                        
            
            time centrifuge-kreport -x $DB_combination_fisher  \
                          $stat_result/LSU_database_final_prediction.report.tsv \
                            > $stat_result/LSU_database_final_prediction.kreprot.tsv

    
       for combine in ITS_LSU ITS LSU;do
       
            rm -f $stat_result/simulating_${num}species_${length}_${replicate}_${combine}_combination_Ture_Positive.txt
            rm -f $stat_result/simulating_${num}species_${length}_${replicate}_${combine}_combination_False_Positive.txt
            rm -f $stat_result/simulating_${num}species_${length}_${replicate}_${combine}_combination_False_negative.txt
            echo "###################################################################################################" >> $stat_result/simulating_${num}species_${length}_${replicate}_${combine}_combination_Ture_Positive.txt
            echo "###################################################################################################" >> $stat_result/simulating_${num}species_${length}_${replicate}_${combine}_combination_False_Positive.txt
            echo "###################################################################################################" >> $stat_result/simulating_${num}species_${length}_${replicate}_${combine}_combination_False_negative.txt
            for level in S G F O C P ; do
                 if [ "$level" = "S" ]; then
                     TL=7
                 elif [ "$level" = "G" ]; then
                     TL=6
                 elif [ "$level" = "F" ]; then
                     TL=5  
                 elif [ "$level" = "O" ]; then
                     TL=4  
                 elif [ "$level" = "C" ]; then
                     TL=3  
                 elif [ "$level" = "P" ]; then
                     TL=2
                 fi
                 echo $TL
         
                 echo $level"######################################" >> $stat_result/simulating_${num}species_${length}_${replicate}_${combine}_combination_Ture_Positive.txt
                 echo $level"######################################" >> $stat_result/simulating_${num}species_${length}_${replicate}_${combine}_combination_False_Positive.txt
                 echo $level"######################################" >> $stat_result/simulating_${num}species_${length}_${replicate}_${combine}_combination_False_negative.txt
                 cat $stat_result/${combine}_database_final_prediction.kreprot.tsv |grep -w "${level}" |awk '{ for(i=1; i<=5; i++){ $i=""} ; print $0 }' |sort -u|while read line; do
                      var1=$(cat $result_dir/simulating_${num}species_${replicate}.taxonomy.txt |cut -d ";" -f $TL | grep "$line" )
                      if [ "$var1" != "" ]; then
                          echo $line >> $stat_result/simulating_${num}species_${length}_${replicate}_${combine}_combination_Ture_Positive.txt
                      else 
                          echo $line >> $stat_result/simulating_${num}species_${length}_${replicate}_${combine}_combination_False_Positive.txt
                      fi
                  done
         
                 cat $result_dir/simulating_${num}species_${replicate}.taxonomy.txt |cut -d ";" -f $TL |cut -d"_" -f 3 |sort -u |while read taxa; do
                      var2=$(cat $stat_result/${combine}_database_final_prediction.kreprot.tsv |grep -w "${level}" |grep "$taxa")
                      if [ "$var2" = "" ]; then
                          echo $taxa >> $stat_result/simulating_${num}species_${length}_${replicate}_${combine}_combination_False_negative.txt
                      fi
                 done
             done
                                  
             rm -f $stat_result/simulating_${num}species_${length}_${replicate}_${combine}_combination_*.final_stat.csv
             rm -f $stat_result/simulating_${num}species_${length}_${replicate}_${combine}_combination.final_stat_combine.csv
             for value in False_negative False_Positive Ture_Positive; do
                 file=$stat_result/simulating_${num}species_${length}_${replicate}_${combine}_combination_${value}.txt
                 stat_file=$stat_result/simulating_${num}species_${length}_${replicate}_${combine}_combination_${value}.final_stat.csv
                 stat_file_combine=$stat_result/simulating_${num}species_${length}_${replicate}_${combine}_combination.final_stat_combine.csv
                 touch $stat_file_combine
                 FN_S=$[$(cat $file |grep "S###########" -A 99999999  |grep "G#########" -B 9999999|wc -l)-2]
                 FN_G=$[$(cat $file |grep "G###########" -A 99999999  |grep "F#########" -B 9999999|wc -l)-2]
                 FN_F=$[$(cat $file |grep "F###########" -A 99999999  |grep "O#########" -B 9999999|wc -l)-2]
                 FN_O=$[$(cat $file |grep "O###########" -A 99999999  |grep "C#########" -B 9999999|wc -l)-2]
                 FN_C=$[$(cat $file |grep "C###########" -A 99999999  |grep "P#########" -B 9999999|wc -l)-2]
                 FN_P=$[$(cat $file |grep "P###########" -A 99999999  |wc -l)-1]
                          
                 echo "level" " " "${value}" >> $stat_file
                 echo "Species" " " $FN_S >> $stat_file
                 echo "Genus" " " $FN_G >> $stat_file
                 echo "Family" " " $FN_F >> $stat_file
                 echo "Order" " " $FN_O >> $stat_file
                 echo "Class" " " $FN_C >> $stat_file
                 echo "Phylum" " " $FN_P >> $stat_file
                          
                 paste $stat_file_combine $stat_file > $stat_result/new
                 mv -f $stat_result/new $stat_file_combine
                 rm -f $stat_result/new     
             done
          done
       done             
   done
done

}


export -f run_miniLength_combination_test
 
run_miniLength_combination_test



rm -f DBs_combine/*
echo -ne "group_id\tspeciesNum\tDatabase\tlength\ttaxa_level\t\tFalse_negative\tFalse_Positive\tTure_Positive\n" >> DBs_combine/hitlength_DBs_combine.txt
for f in stat_result_simulating* ; do
    i=$(echo $f |cut -d"_" -f 3,4,5,6,7,8)
    group_id=$(echo $f |cut -d"_" -f 3,4,5,6,7,8)
    length=$(echo $f |cut -d"_" -f 5)
    speciesNum=$(echo $f |cut -d"_" -f 4)
    
    for DB in ITS LSU ITS_LSU; do
       for taxa_level in Species Genus Family Order Class Phylum; do
          if [ "$taxa_level" = "Species" ]; then
                     TL=2
          elif [ "$taxa_level" = "Genus" ]; then
                     TL=3
          elif [ "$taxa_level" = "Family" ]; then
                     TL=4  
          elif [ "$taxa_level" = "Order" ]; then
                     TL=5  
          elif [ "$taxa_level" = "Class" ]; then
                     TL=6  
          elif [ "$taxa_level" = "Phylum" ]; then
                     TL=7
          fi
          
          value=$(cat $f/${i}_${DB}_combination.final_stat_combine.csv |awk NR==$TL |sed "s/$taxa_level//g")
          echo -ne "$group_id\t$speciesNum\t$DB\t$length\t$taxa_level $value\n" >> DBs_combine/hitlength_DBs_combine.txt
       done
    done
done











#######################################################################################################
#######################################################################################################



run_speciesNum_test () {

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
mkdir $Test_Fungi_RefSeq/SpeciesNum_test
mkdir $Test_Fungi_RefSeq/SpeciesNum_test/stat_result

for length in 140 ; do
   for num in  50 100 200; do
       for replicate in 1 2 3 4 5; do      
           result_dir=$Test_Fungi_RefSeq/SpeciesNum_test/simulating_${num}species_${replicate}
           seq_dir=$Test_Fungi_RefSeq/SpeciesNum_test/simulating_${num}species_${replicate}/splite_seq
           stat_result=$Test_Fungi_RefSeq/SpeciesNum_test/stat_result/stat_result_simulating_${num}species_${length}_${replicate}
           #check if the result directory exist or not
           if [ ! -d "$result_dir/" ];then
              mkdir $result_dir
              mkdir $seq_dir
              mkdir $seq_dir/ITS
              mkdir $seq_dir/LSU
           else
              echo "$result_dir   exist"
           fi
           #check the simulating dataset exist or not, if not, then simulate the from RefSeq
           if [ ! -f "$result_dir/simulating_${num}species_${replicate}.short_read_R1.fastq.gz" ];then
                 #randomly select fungal species from the RefSeq database from simulating
                 cat $Test_Fungi_RefSeq/RefSeq_Fungi/28S_ITS_common_species_unique |shuf -n $num > $result_dir/species_list
                 cat $result_dir/species_list |while read taxa; do
                     genus=$(echo $taxa |cut -d" " -f 1)
                     species=$(echo $taxa |cut -d" " -f 2)
                     cat $Test_Fungi_RefSeq/RefSeq_Fungi/ITS_RefSeq_Fungi.fa  |seqkit seq -w 0 -g -u |grep -A 1 "${taxa}" > $seq_dir/ITS/${genus}_${species}.fasta
                     cat $Test_Fungi_RefSeq/RefSeq_Fungi/28S_fungal_sequences.fa  |seqkit seq -w 0 -g -u |grep -A 1 "${taxa}" > $seq_dir/LSU/${genus}_${species}.fasta
                 done
                 #generate the taxonomy of selected fungal species
                 for fasta in $(ls $seq_dir/ITS/*.fasta) ; do
                     for acc in $(cat $fasta |seqkit seq -i -n);do
                         accession=$acc
                         taxid=$(echo $acc |esummary -db nuccore  | xtract -pattern DocumentSummary -element TaxId)
                         linage=$(echo $taxid | taxonkit lineage --data-dir $DB_taxonomy |taxonkit reformat -P --data-dir $DB_taxonomy --delimiter ";" --format "{k};{p};{c};{o};{f};{g};{s}" | cut -f 1,3)
                         taxonomy=$accession'   '$linage
                         echo $taxonomy >> $result_dir/simulating_${num}species_${replicate}.taxonomy.txt
                     done
                 done
                 
                 for fasta in $(ls $seq_dir/LSU/*.fasta) ; do
                     for acc in $(cat $fasta |seqkit seq -i -n);do
                         accession=$acc
                         taxid=$(echo $acc |esummary -db nuccore  | xtract -pattern DocumentSummary -element TaxId)
                         linage=$(echo $taxid | taxonkit lineage --data-dir $DB_taxonomy |taxonkit reformat -P --data-dir $DB_taxonomy --delimiter ";" --format "{k};{p};{c};{o};{f};{g};{s}" | cut -f 1,3)
                         taxonomy=$accession'   '$linage
                         echo $taxonomy >> $result_dir/simulating_${num}species_${replicate}.taxonomy.txt
                     done
                 done
               #replace sequences id with fungal species name
               #  for sequence in $(ls $seq_dir/ITS/*.fasta); do
               #      seq_taxa=$(echo $sequence |cut -d"/" -f 12 |cut -d "." -f 1)
               #      cat $sequence |seqkit replace -p ".+" -r "${seq_taxa}" > $seq_dir/ITS/$seq_taxa.fa
               #      rm -f $seq_dir/ITS/$seq_taxa.fasta
               #  done
                 
               #  for sequence in $(ls $seq_dir/LSU/*.fasta); do
               #      seq_taxa=$(echo $sequence |cut -d"/" -f 12 |cut -d "." -f 1)
               #      cat $sequence |seqkit replace -p ".+" -r "${seq_taxa}" > $seq_dir/LSU/$seq_taxa.fa
               #      rm -f $seq_dir/LSU/$seq_taxa.fasta
                # done
                 #gernrate the simulating dataset using iss
                 time iss generate --cpus 24 --quiet  --compress \
                       --genomes $seq_dir/ITS/*.fasta \
                       --model novaseq \
                       --abundance exponential \
                       --n_reads 5M \
                       --output $result_dir/simulating_${num}species_${replicate}_ITS.short_read
           #cat $result_dir/simulating_${num}species_${length}_${replicate}_ITS.short_read_abundance.txt |sed 's/ITS/LSU/g' > $result_dir/simulating_${num}species_${length}_${replicate}_LSU.short_read_abundance.txt
                 time iss generate --cpus 24 --quiet  --compress \
                       --genomes $seq_dir/LSU/*.fasta \
                       --model novaseq \
                       --abundance exponential \
                       --n_reads 5M \
                       --output $result_dir/simulating_${num}species_${replicate}_LSU.short_read
            #combine the simulating dataset from ITS and LSU RefSeq database
            cat $result_dir/simulating_${num}species_${replicate}_ITS.short_read_R1.fastq.gz  $result_dir/simulating_${num}species_${replicate}_LSU.short_read_R1.fastq.gz > $result_dir/simulating_${num}species_${replicate}.short_read_R1.fastq.gz
            cat $result_dir/simulating_${num}species_${replicate}_ITS.short_read_R2.fastq.gz  $result_dir/simulating_${num}species_${replicate}_LSU.short_read_R2.fastq.gz > $result_dir/simulating_${num}species_${replicate}.short_read_R2.fastq.gz
            
            rm -f $result_dir/simulating_${num}species_${replicate}_ITS.short_read_R1.fastq.gz  $result_dir/simulating_${num}species_${replicate}_LSU.short_read_R1.fastq.gz
            rm -f $result_dir/simulating_${num}species_${replicate}_ITS.short_read_R2.fastq.gz  $result_dir/simulating_${num}species_${replicate}_LSU.short_read_R2.fastq.gz
            
            else
                echo "$result_dir/simulating_${num}species_${replicate}.short_read_R1.fastq.gz  exist"
            fi
            
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
                    
                    
                    #evaluate the classification performance 
                    rm -f $stat_result/simulating_${num}species_${length}_${replicate}_${DB_cat}_Ture_Positive.txt
                    rm -f $stat_result/simulating_${num}species_${length}_${replicate}_${DB_cat}_False_Positive.txt
                    rm -f $stat_result/simulating_${num}species_${length}_${replicate}_${DB_cat}_False_negative.txt
                    echo "###################################################################################################" >> $stat_result/simulating_${num}species_${length}_${replicate}_${DB_cat}_Ture_Positive.txt
                    echo "###################################################################################################" >> $stat_result/simulating_${num}species_${length}_${replicate}_${DB_cat}_False_Positive.txt
                    echo "###################################################################################################" >> $stat_result/simulating_${num}species_${length}_${replicate}_${DB_cat}_False_negative.txt
                    for level in S G F O C P ; do
                        if [ "$level" = "S" ]; then
                           TL=7
                        elif [ "$level" = "G" ]; then
                           TL=6
                        elif [ "$level" = "F" ]; then
                           TL=5  
                        elif [ "$level" = "O" ]; then
                           TL=4  
                        elif [ "$level" = "C" ]; then
                           TL=3  
                        elif [ "$level" = "P" ]; then
                           TL=2
                        fi
                        echo $TL
         
                        echo $level"######################################" >> $stat_result/simulating_${num}species_${length}_${replicate}_${DB_cat}_Ture_Positive.txt
                        echo $level"######################################" >> $stat_result/simulating_${num}species_${length}_${replicate}_${DB_cat}_False_Positive.txt
                        echo $level"######################################" >> $stat_result/simulating_${num}species_${length}_${replicate}_${DB_cat}_False_negative.txt
                        cat $stat_result/simulating_${num}species_${length}_${replicate}_${DB_cat}.short_read.kreprot.tsv |grep -w "${level}" |awk '{ for(i=1; i<=5; i++){ $i=""} ; print $0 }' |sort -u|while read line; do
                               var1=$(cat $result_dir/simulating_${num}species_${replicate}.taxonomy.txt |cut -d ";" -f $TL | grep "$line" )
                               if [ "$var1" != "" ]; then
                                   echo $line >> $stat_result/simulating_${num}species_${length}_${replicate}_${DB_cat}_Ture_Positive.txt
                               else 
                                   echo $line >> $stat_result/simulating_${num}species_${length}_${replicate}_${DB_cat}_False_Positive.txt
                               fi
                        done
         
                        cat $result_dir/simulating_${num}species_${replicate}.taxonomy.txt |cut -d ";" -f $TL |cut -d"_" -f 3 |sort -u |while read taxa; do
                            var2=$(cat $stat_result/simulating_${num}species_${length}_${replicate}_${DB_cat}.short_read.kreprot.tsv |grep -w "${level}" |grep "$taxa")
                              if [ "$var2" = "" ]; then
                                   echo $taxa >> $stat_result/simulating_${num}species_${length}_${replicate}_${DB_cat}_False_negative.txt
                              fi
                        done
                  done
                                  
                  #generate and combine the classification performance evaluation values
                  rm -f $stat_result/simulating_${num}species_${length}_${replicate}_${DB_cat}_*.stat.csv
                  rm -f $stat_result/simulating_${num}species_${length}_${replicate}_${DB_cat}.stat_combine.csv
                  for value in False_negative False_Positive Ture_Positive; do
                          file=$stat_result/simulating_${num}species_${length}_${replicate}_${DB_cat}_${value}.txt
                          stat_file=$stat_result/simulating_${num}species_${length}_${replicate}_${DB_cat}_${value}.stat.csv
                          stat_file_combine=$stat_result/simulating_${num}species_${length}_${replicate}_${DB_cat}.stat_combine.csv
                          
                          touch $stat_file_combine
                          FN_S=$[$(cat $file |grep "S###########" -A 99999999  |grep "G#########" -B 9999999|wc -l)-2]
                          FN_G=$[$(cat $file |grep "G###########" -A 99999999  |grep "F#########" -B 9999999|wc -l)-2]
                          FN_F=$[$(cat $file |grep "F###########" -A 99999999  |grep "O#########" -B 9999999|wc -l)-2]
                          FN_O=$[$(cat $file |grep "O###########" -A 99999999  |grep "C#########" -B 9999999|wc -l)-2]
                          FN_C=$[$(cat $file |grep "C###########" -A 99999999  |grep "P#########" -B 9999999|wc -l)-2]
                          FN_P=$[$(cat $file |grep "P###########" -A 99999999  |wc -l)-1]
                          
                          echo "level" " " "${value}" >> $stat_file
                          echo "Species" " " $FN_S >> $stat_file
                          echo "Genus" " " $FN_G >> $stat_file
                          echo "Family" " " $FN_F >> $stat_file
                          echo "Order" " " $FN_O >> $stat_file
                          echo "Class" " " $FN_C >> $stat_file
                          echo "Phylum" " " $FN_P >> $stat_file
                          
                          paste $stat_file_combine $stat_file > $stat_result/new
                          mv -f $stat_result/new $stat_file_combine
                          rm -f $stat_result/new     
                  done
                     
            done
       done             
   done
done

}


export -f run_speciesNum_test
 
run_speciesNum_test




rm -f DBs_separated/*
echo -ne "group_id\tspeciesNum\tDatabase\tlength\ttaxa_level\t\tFalse_negative\tFalse_Positive\tTure_Positive\n" >> DBs_separated/hitlength_DBs_separated.txt
for f in stat_result_simulating* ; do
    i=$(echo $f |cut -d"_" -f 3,4,5,6,7,8)
    group_id=$(echo $f |cut -d"_" -f 3,4,5,6,7,8)
    length=$(echo $f |cut -d"_" -f 5)
    speciesNum=$(echo $f |cut -d"_" -f 4| sed 's/species//')
    
    for DB in ITS1 ITS2 LsuD1 LsuD2; do
       for taxa_level in Species Genus Family Order Class Phylum; do
          if [ "$taxa_level" = "Species" ]; then
                     TL=2
          elif [ "$taxa_level" = "Genus" ]; then
                     TL=3
          elif [ "$taxa_level" = "Family" ]; then
                     TL=4  
          elif [ "$taxa_level" = "Order" ]; then
                     TL=5  
          elif [ "$taxa_level" = "Class" ]; then
                     TL=6  
          elif [ "$taxa_level" = "Phylum" ]; then
                     TL=7
          fi
          
          value=$(cat $f/${i}_${DB}.stat_combine.csv |awk NR==$TL |sed "s/$taxa_level//g")
          echo -ne "$group_id\t$speciesNum\t$DB\t$length\t$taxa_level $value\n" >> DBs_separated/hitlength_DBs_separated.txt
       done
    done
done




#########################################################################################################
#########################################################################################################

run_speciesNum_combination_test () {

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


for length in 140; do
   for num in 50 100 200 ; do
       for replicate in 1 2 3 4 5; do      
           result_dir=$Test_Fungi_RefSeq/SpeciesNum_test/simulating_${num}species_${replicate}
           seq_dir=$Test_Fungi_RefSeq/SpeciesNum_test/simulating_${num}species_${replicate}/splite_seq
           stat_result=$Test_Fungi_RefSeq/SpeciesNum_test/stat_result/stat_result_simulating_${num}species_${length}_${replicate}
           
           ITS1_report=$stat_result/simulating_${num}species_${length}_${replicate}_ITS1.short_read.reprot.tsv
           ITS2_report=$stat_result/simulating_${num}species_${length}_${replicate}_ITS2.short_read.reprot.tsv
           LusD1_report=$stat_result/simulating_${num}species_${length}_${replicate}_LsuD1.short_read.reprot.tsv
           LusD2_report=$stat_result/simulating_${num}species_${length}_${replicate}_LsuD2.short_read.reprot.tsv
           
           cat $ITS1_report |cut -f 2 > $stat_result/ITS1_taxID
           cat $ITS2_report |cut -f 2 > $stat_result/ITS2_taxID
           cat $LusD1_report |cut -f 2 > $stat_result/LusD1_taxID
           cat $LusD2_report |cut -f 2 > $stat_result/LusD2_taxID
           
           
           
           rm -f $stat_result/all_database_final_prediction
           cat $stat_result/ITS1_taxID | while read line ; do
               taxID=$line
               var1=$( cat $stat_result/LusD1_taxID | grep -w "${taxID}" )
               var2=$( cat $stat_result/LusD2_taxID | grep -w "${taxID}" )
               var3=$( cat $stat_result/ITS2_taxID | grep -w "${taxID}" )
               if [ "$var1" != "" -o "$var2" != "" -o "$var3" != "" ]; then
                  echo $taxID >> $stat_result/all_database_final_prediction
               fi
            done
            cat $stat_result/ITS2_taxID | while read line ; do
               taxID=$line
               var4=$( cat $stat_result/LusD1_taxID | grep -w "${taxID}" )
               var5=$( cat $stat_result/LusD2_taxID | grep -w "${taxID}" )
               if [ "$var4" != "" -o "$var5" != "" ]; then
                  echo $taxID >> $stat_result/all_database_final_prediction
               fi
            done
            cat $stat_result/LusD1_taxID | while read line ; do
               taxID=$line
               var6=$( cat $stat_result/LusD2_taxID | grep -w "${taxID}" )
               if [ "$var6" ]; then
                  echo $taxID >> $stat_result/all_database_final_prediction
               fi
            done
            
            
            
            
            rm -f $stat_result/ITS1_taxID $stat_result/ITS2_taxID $stat_result/LusD1_taxID $stat_result/LusD2_taxID
                   
            rm -f $stat_result/all_database_final_prediction.report.tsv
            cat $ITS1_report |head -n 1 >> $stat_result/all_database_final_prediction.report.tsv
            cat $stat_result/all_database_final_prediction|sort -u |while read line; do
                taxID=$line
                var1=$(cat $ITS1_report |awk -F "\t"  '{if ($2 == '$taxID') print }')
                var2=$(cat $ITS2_report |awk -F "\t"  '{if ($2 == '$taxID') print }')
                var3=$(cat $LusD1_report |awk -F "\t"  '{if ($2 == '$taxID') print }')
                var4=$(cat $LusD2_report |awk -F "\t"  '{if ($2 == '$taxID') print }')
                if [ "$var1" != "" ]; then
                   cat $ITS1_report |awk -F "\t"  '{if ($2 == '$taxID') print }' >> $stat_result/all_database_final_prediction.report.tsv
                elif [ "$var2" != "" ]; then
                   cat $ITS2_report |awk -F "\t"  '{if ($2 == '$taxID') print }' >> $stat_result/all_database_final_prediction.report.tsv
                elif [ "$var3" != "" ]; then
                   cat $LusD2_report |awk -F "\t"  '{if ($2 == '$taxID') print }' >> $stat_result/all_database_final_prediction.report.tsv
                elif [ "$var4" != "" ]; then
                   cat $LusD1_report |awk -F "\t"  '{if ($2 == '$taxID') print }' >> $stat_result/all_database_final_prediction.report.tsv
                fi
            done

            
            time centrifuge-kreport -x $DB_combination_fisher  \
                          $stat_result/all_database_final_prediction.report.tsv \
                            > $stat_result/all_database_final_prediction.kreprot.tsv

            rm -f $stat_result/simulating_${num}species_${length}_${replicate}_combination_Ture_Positive.txt
            rm -f $stat_result/simulating_${num}species_${length}_${replicate}_combination_False_Positive.txt
            rm -f $stat_result/simulating_${num}species_${length}_${replicate}_combination_False_negative.txt
            echo "###################################################################################################" >> $stat_result/simulating_${num}species_${length}_${replicate}_combination_Ture_Positive.txt
            echo "###################################################################################################" >> $stat_result/simulating_${num}species_${length}_${replicate}_combination_False_Positive.txt
            echo "###################################################################################################" >> $stat_result/simulating_${num}species_${length}_${replicate}_combination_False_negative.txt
            for level in S G F O C P ; do
                 if [ "$level" = "S" ]; then
                     TL=7
                 elif [ "$level" = "G" ]; then
                     TL=6
                 elif [ "$level" = "F" ]; then
                     TL=5  
                 elif [ "$level" = "O" ]; then
                     TL=4  
                 elif [ "$level" = "C" ]; then
                     TL=3  
                 elif [ "$level" = "P" ]; then
                     TL=2
                 fi
                 echo $TL
         
                 echo $level"######################################" >> $stat_result/simulating_${num}species_${length}_${replicate}_combination_Ture_Positive.txt
                 echo $level"######################################" >> $stat_result/simulating_${num}species_${length}_${replicate}_combination_False_Positive.txt
                 echo $level"######################################" >> $stat_result/simulating_${num}species_${length}_${replicate}_combination_False_negative.txt
                 cat $stat_result/all_database_final_prediction.kreprot.tsv |grep -w "${level}" |awk '{ for(i=1; i<=5; i++){ $i=""} ; print $0 }' |sort -u|while read line; do
                      var1=$(cat $result_dir/simulating_${num}species_${replicate}.taxonomy.txt |cut -d ";" -f $TL | grep "$line" )
                      if [ "$var1" != "" ]; then
                          echo $line >> $stat_result/simulating_${num}species_${length}_${replicate}_combination_Ture_Positive.txt
                      else 
                          echo $line >> $stat_result/simulating_${num}species_${length}_${replicate}_combination_False_Positive.txt
                      fi
                  done
         
                 cat $result_dir/simulating_${num}species_${replicate}.taxonomy.txt |cut -d ";" -f $TL |cut -d"_" -f 3 |sort -u |while read taxa; do
                      var2=$(cat $stat_result/all_database_final_prediction.kreprot.tsv |grep -w "${level}" |grep "$taxa")
                      if [ "$var2" = "" ]; then
                          echo $taxa >> $stat_result/simulating_${num}species_${length}_${replicate}_combination_False_negative.txt
                      fi
                 done
             done
                                  
             rm -f $stat_result/simulating_${num}species_${length}_${replicate}_combination_*.final_stat.csv
             rm -f $stat_result/simulating_${num}species_${length}_${replicate}_combination.final_stat_combine.csv
             for value in False_negative False_Positive Ture_Positive; do
                 file=$stat_result/simulating_${num}species_${length}_${replicate}_combination_${value}.txt
                 stat_file=$stat_result/simulating_${num}species_${length}_${replicate}_combination_${value}.final_stat.csv
                 stat_file_combine=$stat_result/simulating_${num}species_${length}_${replicate}_combination.final_stat_combine.csv
                 touch $stat_file_combine
                 FN_S=$[$(cat $file |grep "S###########" -A 99999999  |grep "G#########" -B 9999999|wc -l)-2]
                 FN_G=$[$(cat $file |grep "G###########" -A 99999999  |grep "F#########" -B 9999999|wc -l)-2]
                 FN_F=$[$(cat $file |grep "F###########" -A 99999999  |grep "O#########" -B 9999999|wc -l)-2]
                 FN_O=$[$(cat $file |grep "O###########" -A 99999999  |grep "C#########" -B 9999999|wc -l)-2]
                 FN_C=$[$(cat $file |grep "C###########" -A 99999999  |grep "P#########" -B 9999999|wc -l)-2]
                 FN_P=$[$(cat $file |grep "P###########" -A 99999999  |wc -l)-1]
                          
                 echo "level" " " "${value}" >> $stat_file
                 echo "Species" " " $FN_S >> $stat_file
                 echo "Genus" " " $FN_G >> $stat_file
                 echo "Family" " " $FN_F >> $stat_file
                 echo "Order" " " $FN_O >> $stat_file
                 echo "Class" " " $FN_C >> $stat_file
                 echo "Phylum" " " $FN_P >> $stat_file
                          
                 paste $stat_file_combine $stat_file > $stat_result/new
                 mv -f $stat_result/new $stat_file_combine
                 rm -f $stat_result/new     
             done
       done             
   done
done

}


export -f run_speciesNum_combination_test
 
run_speciesNum_combination_test







rm -f DBs_combine/*
echo -ne "group_id\tspeciesNum\tDatabase\tlength\ttaxa_level\t\tFalse_negative\tFalse_Positive\tTure_Positive\n" >> DBs_combine/hitlength_DBs_combine.txt
for f in stat_result_simulating* ; do
    i=$(echo $f |cut -d"_" -f 3,4,5,6,7,8)
    group_id=$(echo $f |cut -d"_" -f 3,4,5,6,7,8)
    length=$(echo $f |cut -d"_" -f 5)
    speciesNum=$(echo $f |cut -d"_" -f 4)
    
    for DB in ITS LSU ITS_LSU; do
       for taxa_level in Species Genus Family Order Class Phylum; do
          if [ "$taxa_level" = "Species" ]; then
                     TL=2
          elif [ "$taxa_level" = "Genus" ]; then
                     TL=3
          elif [ "$taxa_level" = "Family" ]; then
                     TL=4  
          elif [ "$taxa_level" = "Order" ]; then
                     TL=5  
          elif [ "$taxa_level" = "Class" ]; then
                     TL=6  
          elif [ "$taxa_level" = "Phylum" ]; then
                     TL=7
          fi
          
          value=$(cat $f/${i}_${DB}_combination.final_stat_combine.csv |awk NR==$TL |sed "s/$taxa_level//g")
          echo -ne "$group_id\t$speciesNum\t$DB\t$length\t$taxa_level $value\n" >> DBs_combine/hitlength_DBs_combine.txt
       done
    done
done


