# This is the scripts for generation of mock communities
##########################################################################################################
#The example data collection
#Download the RefSeq from NCBI
###########################################################################################################
wget -c https://ftp.ncbi.nih.gov/blast/db/28S_fungal_sequences.tar.gz
tar -zxvf 28S_fungal_sequences.tar.gz
wget -c https://ftp.ncbi.nih.gov/blast/db/ITS_RefSeq_Fungi.tar.gz
tar -zxvf ITS_RefSeq_Fungi.tar.gz

#Extract sample sequences from the 28S_fungal_sequences
#This step will extract the sequences of species that have both ITS and 28S rRNA sequences

#######################################################
#generate the mock community
#######################################################
#!/bin/bash

conda activate metagenome




mock_community_generation () {
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
mkdir $Test_Fungi_RefSeq/mock_community
mkdir $Test_Fungi_RefSeq/mock_community/classification_result

 #set the number of fungal species in the simulating dataset
   for num in 50 100 200; do
      #set the number of replicates
       for replicate in 1 2 3 4 5; do      
           result_dir=$Test_Fungi_RefSeq/mock_community/simulating_${num}species_${replicate}
           seq_dir=$Test_Fungi_RefSeq/mock_community/simulating_${num}species_${replicate}/splite_seq
           
          #check if the result directory exist or not
           if [ ! -d "$result_dir/" ];then
              mkdir $result_dir
              mkdir $seq_dir
           else
              echo "$result_dir   exist"
           fi
          #check the simulating dataset exist or not, if not, then extract the sequences from RefSeq
           if [ ! -f "$result_dir/simulating_${num}species_${replicate}.short_read_R1.fastq.gz" ];then
                #randomly select fungal species from the RefSeq database
                 cat $Test_Fungi_RefSeq/RefSeq_Fungi/28S_ITS_common_species_unique |shuf -n $num > $result_dir/species_list
                #extract the sequence of selected species
                 cat $result_dir/species_list |while read taxa; do
                     genus=$(echo $taxa |cut -d" " -f 1)
                     species=$(echo $taxa |cut -d" " -f 2)
                     cat $Test_Fungi_RefSeq/RefSeq_Fungi/ITS_RefSeq_Fungi.fa  |seqkit seq -w 0 -g -u |grep -A 1 "${taxa}" > $seq_dir/${genus}_${species}.fasta
                     cat $Test_Fungi_RefSeq/RefSeq_Fungi/28S_fungal_sequences.fa  |seqkit seq -w 0 -g -u |grep -A 1 "${taxa}" >> $seq_dir/${genus}_${species}.fasta
                 done
                #generate the taxonomy of selected fungal species
                 for fasta in $(ls $seq_dir/*.fasta) ; do
                     for acc in $(cat $fasta |seqkit seq -i -n);do
                         accession=$acc
                         taxid=$(echo $acc |esummary -db nuccore  | xtract -pattern DocumentSummary -element TaxId)
                         linage=$(echo $taxid | taxonkit lineage --data-dir $DB_taxonomy |taxonkit reformat -P --data-dir $DB_taxonomy --delimiter ";" --format "{k};{p};{c};{o};{f};{g};{s}" | cut -f 1,3)
                         taxonomy=$accession'   '$linage
                         echo $taxonomy >> $result_dir/simulating_${num}species_${replicate}.taxonomy.txt
                     done
                 done
                 

                 #gernrate the simulating dataset using iss
                 time iss generate --cpus 24 --quiet  --compress \
                       --draft $seq_dir/*.fasta \
                       --model novaseq \
                       --abundance exponential \
                       --n_reads 5M \
                       --output $result_dir/simulating_${num}species_${replicate}.short_read
                 #cat $result_dir/simulating_${num}species_${length}_${replicate}_ITS.short_read_abundance.txt |sed 's/ITS/LSU/g' > $result_dir/simulating_${num}species_${length}_${replicate}_LSU.short_read_abundance.txt
                 
            else
                echo "$result_dir/simulating_${num}species_${replicate}.short_read_R1.fastq.gz  exist"
            fi
        done
    done
}


