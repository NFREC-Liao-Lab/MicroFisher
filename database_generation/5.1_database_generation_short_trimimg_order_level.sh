#!/bin/bash -l
#SBATCH --job-name=LSU_database_generation_preprocess_trimming      # Job name
#SBATCH --mail-type=all            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=wanghaihua@ufl.edu       # Where to send mail	
#SBATCH --ntasks=1                      # Number of MPI ranks
#SBATCH --cpus-per-task=1               # Number of cores per MPI rank 
#SBATCH --nodes=1                       # Number of nodes
#SBATCH --ntasks-per-node=1             # How many tasks on each node
#SBATCH --ntasks-per-socket=1           # How many tasks on each CPU or socket
#SBATCH --mem-per-cpu=3G             # Memory per core
#SBATCH --time=30-00:00:00                 # Time limit hrs:min:sec
#SBATCH --output=LSU_database_generation_preprocess_trimming_%j.log     # Standard output and error log

pwd; hostname; date


CPU="18"


#module load faSomeRecords
module load ucsc/20210803
module load clustalo
module load rstudio/1.4
module load seqkit/0.10.2
module load clustalw/2.1
#module load mega/7.0.26
module load gcc/5.2.0
module load parallel/20150122


#This is the code for generation of mulit- short database D1D2 region
work_dir="/home/wanghaihua/all_data/database_generation/fungal_lsu_database/result_dir"
cd $work_dir
# SILVA_138.1_LSUParc_tax_silva_LSU_DNA_filter_final_uniq_assigned.fasta 
# SILVA_138.1_LSUParc_tax_silva_LSU_DNA_filter_final_uniq_assigned_nmdump_removed.fasta 
# accession2taxid_fungi_assigned 
# taxid2lineage_fungi 


######################################################################
#trim the seqences in database into 3 short fragments

#1. extract taxonomy names in Order level
grep -r ";o__" taxid2lineage_fungi   > taxid2lineage_fungi_0

while read line
do
echo $line  |cut -d ";" -f 4 |cut -d "_" -f 3  >> order_level_ids_0
done < taxid2lineage_fungi_0

#count the non-order taxonomy number
cat order_level_ids_0  | wc -l

cat order_level_ids_0 |sort -u > order_level_ids
rm order_level_ids_0
rm taxid2lineage_fungi_0

mkdir $work_dir/sep_process
mkdir $work_dir/short_seqs
mkdir $work_dir/seqLogo
mkdir $work_dir/aligned_aln


###extract sequence at order level
mkdir $work_dir/sep_process/seq_less_10_Order
for Order in $(cat order_level_ids)
do

  ###setup the work file named with Order name
  mkdir $work_dir/sep_process/${Order}_Order
  
  ###extract the order level's taxid
  grep -r "${Order}" taxid2lineage_fungi |awk '{print $1}' > $work_dir/sep_process/${Order}_Order/${Order}_taxid
  
  ###extract the order level's accession number 
  for taxid in $(cat $work_dir/sep_process/${Order}_Order/${Order}_taxid)
  do
     grep -w "${taxid}" accession2taxid_fungi_assigned | awk '{print $1}'|cut -d ":" -f 2 |sort -u >>  $work_dir/sep_process/${Order}_Order/${Order}_accession_num
  done
  
  ###extract the order level's sequences
     faSomeRecords conbined_remove_low_quality.fasta \
                                          $work_dir/sep_process/${Order}_Order/${Order}_accession_num \
                                          $work_dir/sep_process/${Order}_Order/${Order}.fasta
                                          
  #combine the order's sequences which sequence number less than 10
  if [[ "$(cat $work_dir/sep_process/${Order}_Order/${Order}.fasta |grep -c ">" )" -le 10 ]] ;
  then
  cat $work_dir/sep_process/${Order}_Order/${Order}.fasta >> $work_dir/sep_process/seq_less_10_Order/seq_less_10.fasta
  rm -rf $work_dir/sep_process/${Order}_Order
  fi
done
  
  
#multiple sequences alignment

  MSA () {
  work_dir="/home/wanghaihua/all_data/database_generation/fungal_lsu_database/result_dir"
  i=$1
  Order=$(basename "$i" _Order)
  
  if [ ! -f "$work_dir/sep_process/${Order}_Order/${Order}.aln" ]; then
   #remove the dumplications
   seqkit rmdup --by-name $work_dir/sep_process/${Order}_Order/${Order}.fasta > $work_dir/sep_process/${Order}_Order/${Order}_duprm.fasta

   #clustalo --force  -i $work_dir/sep_process/${Order}_Order/${Order}.fasta -o $work_dir/sep_process/${Order}_Order/${Order}.aln
     clustalw -infile=$work_dir/sep_process/${Order}_Order/${Order}_duprm.fasta \
              -align \
              -outfile=$work_dir/sep_process/${Order}_Order/${Order}.aln \
              -maxseqlen=100000000 -OUTPUT=FASTA -QUIET
     
     cp $work_dir/sep_process/${Order}_Order/${Order}.aln  $work_dir/aligned_aln
     cp $work_dir/sep_process/${Order}_Order/${Order}.aln $work_dir/sep_process/${Order}_Order/order_seq_aligned.aln
     
  else 
       echo "${Order} has been aligned"
  fi        
  
  }
  
export -f MSA
time parallel -j 20 MSA ::: $(ls $work_dir/sep_process)
  
  


printf "Order_name","sequence_num","algned_seq_length","database_1_start","database_1_end","database_1_mean_Bits","database_1_length","database_2_start","database_2_end","database_2_mean_Bits","database_2_length" >> $work_dir/seq_info_in_order_level.txt

printf "file","format  type","num_seqs","sum_len","min_len","avg_len","max_len""\n" >> $work_dir/seq_stat_in_order_level.txt


for i in $(ls $work_dir/sep_process)
do
  Order=$(basename "$i" _Order)
  
  ####write the sequence information
  printf "\n"${Order} >> $work_dir/seq_info_in_order_level.txt
  printf "," >> $work_dir/seq_info_in_order_level.txt
  # count the number of sequences in each order
  printf $(cat $work_dir/sep_process/${Order}_Order/${Order}.fasta |grep -c ">" ) >> $work_dir/seq_info_in_order_level.txt
  printf "," >> $work_dir/seq_info_in_order_level.txt

  ###processing
  #1. seqlogo
     cd $work_dir/sep_process/${Order}_Order
     Rscript  /home/wanghaihua/scripts_file/LSU_database_generation/6_seqlogo_plot.R
     mv -f sequence_logo.pdf $work_dir/seqLogo/${Order}_seqLogo.pdf
     printf $(cat aligned_length.txt |awk 'NR==2{print $2}') >> $work_dir/seq_info_in_order_level.txt
     printf "," >> $work_dir/seq_info_in_order_level.txt
  
  #2. position calculation
     Rscript /home/wanghaihua/scripts_file/LSU_database_generation/7_sliding_window_comparation.r
     if [[ "$(cat cut_start_1.txt)" = "" ]] || [[ "$(cat cut_start_2.txt)" = "" ]] ;
     then
     Rscript /home/wanghaihua/scripts_file/LSU_database_generation/8_sliding_window_comparation.r
     fi
  #3. write the selected position and information
     cat cut_start_1.txt |awk 'NR==2{print $2}' > cut_start_1
     printf "$(cat cut_start_1.txt |awk 'NR==2{print $2}')", >> $work_dir/seq_info_in_order_level.txt
     cat cut_end_1.txt |awk 'NR==2{print $2}' > cut_end_1
     printf "$(cat cut_end_1.txt |awk 'NR==2{print $2}')", >> $work_dir/seq_info_in_order_level.txt
     printf "$(cat Bits_mean_1.txt |awk 'NR==2{print $2}')", >> $work_dir/seq_info_in_order_level.txt
     printf "$(cat length_D1.txt |awk 'NR==2{print $2}')", >> $work_dir/seq_info_in_order_level.txt
    
    
     cat cut_start_2.txt |awk 'NR==2{print $2}' > cut_start_2
     printf "$(cat cut_start_2.txt |awk 'NR==2{print $2}')", >> $work_dir/seq_info_in_order_level.txt
     cat cut_end_2.txt |awk 'NR==2{print $2}' > cut_end_2
     printf "$(cat cut_end_2.txt |awk 'NR==2{print $2}')", >> $work_dir/seq_info_in_order_level.txt
     printf "$(cat Bits_mean_2.txt |awk 'NR==2{print $2}')", >> $work_dir/seq_info_in_order_level.txt
     printf "$(cat length_D2.txt |awk 'NR==2{print $2}')", >> $work_dir/seq_info_in_order_level.txt
  
     rm -f cut_start_1.txt cut_end_1.txt  cut_start_2.txt cut_end_2.txt
     
   #4. get the short sequences  
     cat order_seq_aligned.aln | seqkit subseq -r $(cat cut_start_1):$(cat cut_end_1) > database_1.aln
     cat order_seq_aligned.aln | seqkit subseq -r $(cat cut_start_2):$(cat cut_end_2) > database_2.aln

     cat database_1.aln | seqkit seq -g -l |seqkit seq -m 100 > ${Order}_short_1.fasta
     cat database_2.aln | seqkit seq -g -l |seqkit seq -m 100 > ${Order}_short_2.fasta

  
     cat ${Order}_short_1.fasta > $work_dir/short_seqs/${Order}_short_1.fasta
     cat ${Order}_short_2.fasta > $work_dir/short_seqs/${Order}_short_2.fasta

  #5. get the sequence stat

  seqkit stat ${Order}_short_1.fasta |awk 'NR==2' >> $work_dir/seq_stat_in_order_level.txt
  seqkit stat ${Order}_short_2.fasta |awk 'NR==2' >> $work_dir/seq_stat_in_order_level.txt

  cd $work_dir

done