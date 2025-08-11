# This is the scripts for generation of mock communities
##########################################################################################################
#The example data collection
#Download the RefSeq from NCBI
###########################################################################################################
wget -c https://ftp.ncbi.nih.gov/blast/db/28S_fungal_sequences.tar.gz
tar -zxvf 28S_fungal_sequences.tar.gz
wget -c https://ftp.ncbi.nih.gov/blast/db/ITS_RefSeq_Fungi.tar.gz
tar -zxvf ITS_RefSeq_Fungi.tar.gz
##Fungal genomes were downloaded from JGI


#This step will extract the sequences of species that have both ITS and 28S rRNA sequences, and Fungal genome list


cat FunDB_genomes_info_FungalTraits_2022_06_01.tsv | while read line 
  do  
  echo ${line}
  species=$(echo $line |awk '{print $3, $4}')
  ITS_28S=$(cat 28S_ITS_common_species_unique | grep "$species")
  echo $ITS_28S
  if [ -z $ITS_28S ];then
      echo "${species} is not common found in genome and RefSeq"
  else
      echo "{$species} is common found in genome and RefSeq"
      echo $line >> 28S_ITS_jgiGenome_common_species
  fi

 done


cat 28S_ITS_jgiGenome_common_species | while read line 
do
  echo $line
  taxa=$(echo $line |awk '{print $3, $4}')
  genus=$(echo $taxa |cut -d" " -f 1)
  species=$(echo $taxa |cut -d" " -f 2)
  genome_folder=$(echo $line |awk '{print $2}')
  target_folder="/home/microbiome/data_storage/SATA3/Fisher_test/Test_Fungi_RefSeq/RefSeq_Fungi/Genome_RefSeq_folder"
  file_folder="/home/microbiome/data_storage/SATA3/Fungi_genome/fungal_genomes"
  genome_file=$(echo $file_folder/$genome_folder/${genome_folder}_AssemblyScaffolds.fasta.gz)
  if [ -f $genome_file ]; then
     cp $genome_file $target_folder/${genus}_${species}_AssemblyScaffolds.fasta.gz
     gunzip $target_folder/${genus}_${species}_AssemblyScaffolds.fasta.gz

     cat ITS_RefSeq_Fungi.fa  |seqkit seq -w 0 -g -u |grep -A 1 "${taxa}" > $target_folder/${genus}_${species}_RefSeq_ITS.fasta
     for i in {1..2000}; do cat $target_folder/${genus}_${species}_RefSeq_ITS.fasta; done > $target_folder/${genus}_${species}_RefSeq_ITS_500x.fasta
     awk 'BEGIN{RS=">"; ORS=""} NR==1 {print ">"$0; next} {seq=""; n=split($0,a,"\n"); for(i=2;i<=n;i++) seq=seq a[i] "\n"; if(!seen[seq]++){print ">"$0} else {print seq}}' $target_folder/${genus}_${species}_RefSeq_ITS_500x.fasta > $target_folder/${genus}_${species}_RefSeq_ITS_500x_1.fasta

     cat 28S_fungal_sequences.fa  |seqkit seq -w 0 -g -u |grep -A 1 "${taxa}" > $target_folder/${genus}_${species}_RefSeq_28S.fasta
     for i in {1..2000}; do cat $target_folder/${genus}_${species}_RefSeq_28S.fasta; done > $target_folder/${genus}_${species}_RefSeq_28S_500x.fasta
     awk 'BEGIN{RS=">"; ORS=""} NR==1 {print ">"$0; next} {seq=""; n=split($0,a,"\n"); for(i=2;i<=n;i++) seq=seq a[i] "\n"; if(!seen[seq]++){print ">"$0} else {print seq}}' $target_folder/${genus}_${species}_RefSeq_28S_500x.fasta > $target_folder/${genus}_${species}_RefSeq_28S_500x_1.fasta

     cat $target_folder/${genus}_${species}_AssemblyScaffolds.fasta  $target_folder/${genus}_${species}_RefSeq_ITS_500x_1.fasta $target_folder/${genus}_${species}_RefSeq_28S_500x_1.fasta > $target_folder/${genus}_${species}_genome_combined.fasta
     
     #rename the seq
     filename=${genus}_${species}_genome_combined
     awk -v fn=$filename '
    BEGIN {count = 1}
    /^>/ {print ">" fn "_" count++; next}
    {print}
    ' $target_folder/${genus}_${species}_genome_combined.fasta > $target_folder/${genus}_${species}_genome_combined_renamed.fasta


     echo $line >> 28S_ITS_jgiGenome_common_species_update
  fi
done

#######################################################
#generate the mock community
#######################################################
#!/bin/bash

conda activate metagenome

for replicate in 1 2 3 4 5 6 7 8 9 10; do
  species_num=20
  result_dir=/home/microbiome/data_storage/SATA3/Fisher_test/Test_Fungi_RefSeq/mock_community_Genome
  seq_dir=/home/microbiome/data_storage/SATA3/Fisher_test/Test_Fungi_RefSeq/mock_community_Genome/splite_seq

#generate the fungal species list
  cat 28S_ITS_jgiGenome_common_species_update |shuf -n $species_num > $result_dir/species_list_replicate_$replicate

  cat $result_dir/species_list_replicate_$replicate | while read line 
  do
    echo $line
    taxa=$(echo $line |awk '{print $3, $4}')
    genus=$(echo $taxa |cut -d" " -f 1)
    species=$(echo $taxa |cut -d" " -f 2)
    target_folder="/home/microbiome/data_storage/SATA3/Fisher_test/Test_Fungi_RefSeq/RefSeq_Fungi/Genome_RefSeq_folder"

    echo "$target_folder/${genus}_${species}_genome_combined_renamed.fasta" >> $result_dir/file_list_replicate_$replicate
  done
  sed ':a;N;$!ba;s/\n/ /g' -i  $result_dir/file_list_replicate_$replicate

  genomes=$(cat $result_dir/file_list_replicate_$replicate)
  time iss generate --cpus 24 --quiet  --compress \
                       --draft $genomes \
                       --model novaseq \
                       --abundance exponential \
                       --n_reads 2000M \
                       --output $result_dir/simulating_${species_num}species_${replicate}.short_read

done














