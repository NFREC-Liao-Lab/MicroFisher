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
ITS1_fisher=$ITS_DBs/ITS1_fisher
ITS2_fisher=$ITS_DBs/ITS2_fisher
LsuD1_fisher=$LSU_D1D2_DBs_new/LSU_D1_fisher_new
LsuD2_fisher=$LSU_D1D2_DBs_new/LSU_D2_fisher_new
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
                     cat $Test_Fungi_RefSeq/RefSeq_Fungi/ITS_RefSeq_Fungi.fa  |seqkit seq -w 0 -g -u |grep -A 1 "${taxa}" >> $seq_dir/ITS/${genus}_${species}.fasta
                     cat $Test_Fungi_RefSeq/RefSeq_Fungi/28S_fungal_sequences.fa  |seqkit seq -w 0 -g -u |grep -A 1 "${taxa}" >> $seq_dir/LSU/${genus}_${species}.fasta
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
                       --model MiSeq \
                       --abundance exponential \
                       --n_reads 40M \
                       --output $result_dir/simulating_${num}species_${replicate}_ITS.short_read
           #cat $result_dir/simulating_${num}species_${length}_${replicate}_ITS.short_read_abundance.txt |sed 's/ITS/LSU/g' > $result_dir/simulating_${num}species_${length}_${replicate}_LSU.short_read_abundance.txt
                 time iss generate --cpus 24 --quiet  --compress \
                       --genomes $seq_dir/LSU/*.fasta \
                       --model MiSeq \
                       --abundance exponential \
                       --n_reads 40M \
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
                   time centrifuge -p 24 -k 1  --min-hitlen $length -x $DB  \
                          -1  $result_dir/simulating_${num}species_${replicate}.short_read_R1.fastq.gz  \
                          -2  $result_dir/simulating_${num}species_${replicate}.short_read_R2.fastq.gz  \
                          -S  $stat_result/simulating_${num}species_${length}_${replicate}_${DB_cat}.short_read.result \
                          --report-file $stat_result/simulating_${num}species_${length}_${replicate}_${DB_cat}.short_read.reprot.tsv 
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
         
                        cat $stat_result/simulating_${num}species_${replicate}.taxonomy.txt |cut -d ";" -f $TL |cut -d"_" -f 3 |sort -u |while read taxa; do
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




