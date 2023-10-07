#!/bin/bash -l
#SBATCH --job-name=LSU_database_generation_preprocess      # Job name
#SBATCH --mail-type=all            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=wanghaihua@ufl.edu       # Where to send mail	
#SBATCH --ntasks=100                      # Number of MPI ranks
#SBATCH --cpus-per-task=9               # Number of cores per MPI rank 
#SBATCH --nodes=1                       # Number of nodes
#SBATCH --ntasks-per-node=1             # How many tasks on each node
#SBATCH --ntasks-per-socket=1           # How many tasks on each CPU or socket
#SBATCH --mem-per-cpu=7G             # Memory per core
#SBATCH --time=30-00:00:00                 # Time limit hrs:min:sec
#SBATCH --output=LSU_database_generation_preprocess_%j.log     # Standard output and error log

pwd; hostname; date


CPU="18"


module load metaxa2
module load itsx/1.1b
#module load seqkit
module load clustalo
module load conda
conda activate RNASeq

#This is the code for generation of mulit- short database D1D2 region

apt-get install ncbi-blast+-legacy
###############################################
      
###############################################
#splite the big file
seqkit split2 -s 100000 fungi_nucleotide.fasta --out-dir splited_seq --threads 24


#pre-process the sequences with primers for ITS


extract_with_primer () {
ITS1_fwd="TCCGTAGGTGAACCTGCGG"
ITS1_rev="CCGCAGGTTCACCTACGGA"
ITS4_rev="TCCTCCGCTTATTGATATGC"
ITS4_fwd="GCATATCAATAAGCGGAGGA"
ITS2_rev="GCTGCGTTCTTCATCGATGC"
ITS3_fwd="GCATCGATGAAGAACGCAGC"

ITS_extract=/home/microbiome/data_storage/SATA2/fungal_nucleotide/ITS_regions/ITS_extract_by_primers

i=$1
base=$(basename "$i" .fasta)

seqkit seq --upper-case -g -w 0 $i > $base.organized.fasta

grep -B 1  "${ITS1_fwd}"  $base.organized.fasta >> $base.organized_ITS.fasta
grep -B 1  "${ITS1_rev}"  $base.organized.fasta >> $base.organized_ITS.fasta
grep -B 1  "${ITS4_rev}"  $base.organized.fasta >> $base.organized_ITS.fasta
grep -B 1  "${ITS4_fwd}"  $base.organized.fasta >> $base.organized_ITS.fasta
grep -B 1  "${ITS2_rev}"  $base.organized.fasta >> $base.organized_ITS.fasta
grep -B 1  "${ITS3_fwd}"  $base.organized.fasta >> $base.organized_ITS.fasta

seqkit rmdup -n  $base.organized_ITS.fasta -w 0 > $base.organized_ITS_rmdup.fasta

}

 
export -f extract_with_primer

time parallel -j 24 --eta --load 99% --noswap  extract_with_primer ::: $(ls *.fasta)





#remove sequence > 100kb and remove whole genome shotgun sequences
remove_long_seq () {
i=$1
base=$(basename "$i" .organized_ITS_rmdup.fasta)
seqkit grep -n -r -v -p   "whole genome shotgun"  $i -w 0 > $base.organized_ITS_rmdup_rmgenome.fasta
seqkit seq --remove-gaps --max-len 5000 --min-len 150 -w 0 $base.organized_ITS_rmdup_rmgenome.fasta > $base.organized_ITS_rmdup_rmgenome_rmlongSeq.fasta
}

export -f remove_long_seq

time parallel -j 24 --eta --load 99% --noswap  remove_long_seq ::: $(ls *.organized_ITS_rmdup.fasta)






#Remove the non-LSU region by using the ITSx
#https://github.com/ncbi/ITSx/blob/master/ITSx%20User's%20Guide.pdf

extract_rRNA () {

i=$1
base=$(basename $i .organized_ITS_rmdup_rmlongSeq.fasta)
ITSx -i $i \
     -o $base.rRNA_ITS \
     -t Fungi \
     --save_regions ITS1,5.8S,ITS2 \
     --only_full F \
     --selection_priority sum \
     --complement T \
     --truncate T \
     --not_found T \
     --preserve T \
     --cpu 24 --multi_thread T

}

export -f extract_rRNA

time parallel -j 12 --eta --load 99% --noswap  extract_rRNA ::: $(ls *.organized_ITS_rmdup_rmlongSeq.fasta)







##############################
#code for ITS regin selectin
ITS1-ITS4:600-700
ITS1-ITS3：200-350
ITS3-ITS4: 350-400
###############################

get_ITS_region () {
ITS1_fwd_start="TCCGTAGGTGAACCTGCGG"
ITS2_rev_mid="GCTGCGTTCTTCATCGATGC"
ITS3_fwd_mid="GCATCGATGAAGAACGCAGC"
ITS4_rev_end="TCCTCCGCTTATTGATATGC"
ITS4_fwd_end="GCATATCAATAAGCGGAGGA"

i=$1
base=$(basename $i .organized_ITS_rmdup_rmgenome_rmlongSeq.fasta)
while read line
do
   var1=$(echo $line | grep ">" )
   var2=$(echo $line | grep "${ITS1_fwd_start}" )
   var3=$(echo $line | grep "${ITS3_fwd_mid}" )
   var4=$(echo $line | grep "${ITS4_fwd_end}" )
   
if [[ "$var1" != "" ]];
then
    seq_name="$line"
    acc=$(echo $seq_name | cut -d" " -f 1)
    #echo $line >> $result_file
else
   if [[ "$var2" != "" ]];
   then
      if [[ "$var3" != "" ]];
      then
         if [[ "$var4" != "" ]];
         then
             echo -ne "$acc\n$line\n" |seqkit amplicon -F $ITS1_fwd_start -R $ITS2_rev_mid  >> $base.ITS1.fasta
             echo -ne "$acc\n$line\n" |seqkit amplicon -F $ITS3_fwd_mid -R $ITS4_rev_end  >> $base.ITS2.fasta
         else
             echo -ne "$acc\n$line\n" |seqkit amplicon -F $ITS1_fwd_start -R $ITS2_rev_mid  >> $base.ITS1.fasta
             echo -ne "$acc\n$line\n" |seqkit amplicon -F $ITS3_fwd_mid -r 1:450 >> $base.ITS2.fasta
         fi
      else
         echo -ne "$acc\n$line\n" |seqkit amplicon -F $ITS1_fwd_start -r 1:450 >> $base.ITS1.fasta
      fi
   else
      if [[ "$var3" != "" ]];
      then
         if [[ "$var4" != "" ]];
         then
             echo -ne "$acc\n$line\n" |seqkit amplicon -R $ITS2_rev_mid -r -450:-1 >> $base.ITS1.fasta
             echo -ne "$acc\n$line\n" |seqkit amplicon -F $ITS3_fwd_mid -R $ITS4_rev_end  >> $base.ITS2.fasta
         else
             echo -ne "$acc\n$line\n" |seqkit amplicon -R $ITS2_rev_mid -r -450:-1 >> $base.ITS1.fasta
             echo -ne "$acc\n$line\n" |seqkit amplicon -F $ITS3_fwd_mid -r 1:450 >> $base.ITS2.fasta
         fi
      else
         if [[ "$var4" != "" ]];
         then
            echo -ne "$acc\n$line\n" |seqkit amplicon -R $ITS4_rev_end -r -450:-1 >> $base.ITS2.fasta
         fi
      fi
    fi
    
fi

done < $i

cat $base.ITS1.fasta |seqkit seq -m 120 -M 450 -w 0 --remove-gaps > $base.ITS1.modified.fasta
cat $base.ITS2.fasta |seqkit seq -m 120 -M 450 -w 0 --remove-gaps > $base.ITS2.modified.fasta
}

export -f get_ITS_region

time parallel -j 24 --eta --load 99% --noswap  get_ITS_region ::: $(ls *.organized_ITS_rmdup_rmgenome_rmlongSeq.fasta)















get_ITS_region_2 () {
ITS1_fwd_start="TCCGTAGGTGAACCTGCGG"
ITS2_rev_mid="GCTGCGTTCTTCATCGATGC"
ITS3_fwd_mid="GCATCGATGAAGAACGCAGC"
ITS4_rev_end="TCCTCCGCTTATTGATATGC"
ITS4_fwd_end="GCATATCAATAAGCGGAGGA"

i=$1
base=$(basename $i .organized_ITS_rmdup_rmgenome_rmlongSeq.fasta)
while read line
do
   var1=$(echo $line | grep ">" )
   var2=$(echo $line | grep "${ITS1_fwd_start}" )
   var3=$(echo $line | grep "${ITS3_fwd_mid}" )
   var4=$(echo $line | grep "${ITS4_fwd_end}" )
   
if [[ "$var1" != "" ]];
then
    seq_name="$line"
    acc=$(echo $seq_name | cut -d" " -f 1)
    #echo $line >> $result_file
else
   if [[ "$var2" != "" ]];
   then
      if [[ "$var3" != "" ]];
      then
         if [[ "$var4" != "" ]];
         then
             echo -ne "$acc\n$line\n" |seqkit amplicon -F $ITS1_fwd_start -R $ITS2_rev_mid  >> $base.ITS1_2.fasta
             echo -ne "$acc\n$line\n" |seqkit amplicon -F $ITS3_fwd_mid -R $ITS4_rev_end  >> $base.ITS2_2.fasta
         else
             echo -ne "$acc\n$line\n" |seqkit amplicon -F $ITS1_fwd_start -R $ITS2_rev_mid  >> $base.ITS1_2.fasta
         fi
      fi
   else
      if [[ "$var3" != "" ]];
      then
         if [[ "$var4" != "" ]];
         then
             echo -ne "$acc\n$line\n" |seqkit amplicon -F $ITS3_fwd_mid -R $ITS4_rev_end  >> $base.ITS2_2.fasta
         fi
      fi
    fi
    
fi

done < $i

cat $base.ITS1_2.fasta |seqkit seq -m 120 -M 450 -w 0 --remove-gaps > $base.ITS1_2.modified.fasta
cat $base.ITS2_2.fasta |seqkit seq -m 120 -M 450 -w 0 --remove-gaps > $base.ITS2_2.modified.fasta
}

export -f get_ITS_region_2

time parallel -j 24 --eta --load 99% --noswap  get_ITS_region_2 ::: $(ls *.organized_ITS_rmdup_rmgenome_rmlongSeq.fasta)










###############################################
#filter the LSU region by using one of the D1D2 primers, meanwhile, trim the sequence based on the primers




##############################################################################################
#extract LSU D1 AND D2 region D1:175;D2:250
##############################################################################################
extract_rRNA_lsu () {

i=$1
base=$(basename $i .fasta)

     metaxa2 -i $i \
        -o $base.rRNA.lsu \
        -f fasta \
        -t e \
        -g lsu \
        --complement T \
        --truncate T \
        --mode a \
        --selection_priority sum \
        --not_found T \
        --preserve T \
        --cpu 18 --multi_thread T 
        

}

export -f extract_rRNA_lsu

nohup time parallel -j 2 --eta --load 90% --noswap  extract_rRNA_lsu ::: $(ls *.fasta) &







extract_LSU_D1D2 () {

forward_LR0R_primer="ACCCGCTGAACTTAAGC"
forward_CTB6_primer="GCATATCAATAAGCGGAGG"

reverse_LF340_1_primer="TACTTGTGCGCTATCGG"
reverse_LF340_2_primer="TACTTGTTCGCTATCGG"
forward_LF340_1_primer="CCGATAGCGCACAAGTA"
forward_LF340_2_primer="CCGATAGCGAACAAGTA"

forward_LR3R_primer="GTCTTGAAACACGGACC"
reverse_LR3R_primer="GGTCCGTGTTTCAAGAC"

reverse_LR3_primer="CCGTGTTTCAAGACGGG"
forward_LR3_primer="CCCGTCTTGAAACACGG"

i=$1
base=$(basename $i .fasta)

cat $i |seqkit seq --only-id -w 0 -g  -u > $base.singline.fasta
rm -f $base.D1.result.fasta $base.D2.result.fasta
while read line
do
   var1=$(echo $line | grep ">" )
   var2=$(echo $line | grep "${forward_CTB6_primer}" )
   var3=$(echo $line | grep "${forward_LF340_1_primer}" )
   var4=$(echo $line | grep "${forward_LF340_2_primer}" )
   var5=$(echo $line | grep "${forward_LR3R_primer}" )   
   
if [ "$var1" != "" ];
then
    seq_name="$line"
    #echo $line >> $result_file
else
    if [ "$var2" != "" ];
    then
        if [ "$var3" != "" ];
        then
            if [ "$var5" != "" ];
            then
                echo -ne "$seq_name\n$line\n" |seqkit amplicon --quiet -F $forward_CTB6_primer -R $reverse_LF340_1_primer -r 70:-70 >> $base.D1.result.fasta 
                echo -ne "$seq_name\n$line\n" |seqkit amplicon --quiet -F $forward_LF340_1_primer -R $reverse_LR3R_primer -r 90:-80 >> $base.D2.result.fasta 
            else
                echo -ne "$seq_name\n$line\n" |seqkit amplicon --quiet -F $forward_CTB6_primer -R $reverse_LF340_1_primer -r 70:-70 >> $base.D1.result.fasta
                echo -ne "$seq_name\n$line\n" |seqkit amplicon --quiet -F $forward_LF340_1_primer -r 90:340 >> $base.D2.result.fasta
            fi
        elif [ "$var4" != "" ];
        then
            if [ "$var5" != "" ];
            then
                echo -ne "$seq_name\n$line\n" |seqkit amplicon --quiet -F $forward_CTB6_primer -R $reverse_LF340_2_primer -r 70:-70 >> $base.D1.result.fasta
                echo -ne "$seq_name\n$line\n" |seqkit amplicon --quiet -F $forward_LF340_2_primer -R $reverse_LR3R_primer -r 90:-80 >> $base.D2.result.fasta
            else
                echo -ne "$seq_name\n$line\n" |seqkit amplicon --quiet -F $forward_CTB6_primer -R $reverse_LF340_2_primer -r 70:-70 >> $base.D1.result.fasta
                echo -ne "$seq_name\n$line\n" |seqkit amplicon --quiet -F $forward_LF340_2_primer -r 90:340 >> $base.D2.result.fasta
            fi
        else
           var6=$(echo $line | grep "${forward_LR0R_primer}" )
           if [ "$var6" != "" ];
           then
              echo -ne "$seq_name\n$line\n" |seqkit amplicon --quiet -F $forward_CTB6_primer -r 70:245 >> $base.D1.result.fasta
           fi
        fi
    else
       if [ "$var3" != "" ];
       then
           if [ "$var5" != "" ];
           then
               echo -ne "$seq_name\n$line\n" |seqkit amplicon --quiet -F $forward_LF340_1_primer -r -245:-70 -f >> $base.D1.result.fasta
               echo -ne "$seq_name\n$line\n" |seqkit amplicon --quiet -F $forward_LF340_1_primer -R $reverse_LR3R_primer -r 90:-80 >> $base.D2.result.fasta
           else
               SeqLength=$(echo -ne "$seq_name\n$line\n" |seqkit fx2tab -l| awk '{print $3}')
               if  [ "$SeqLength" -lt 700 ];
               then
                   echo -ne "$seq_name\n$line\n" |seqkit amplicon --quiet -F $forward_LF340_1_primer -r -245:-70 -f >> $base.D1.result.fasta
                   echo -ne "$seq_name\n$line\n" |seqkit amplicon --quiet -F $forward_LF340_1_primer -r 90:340 >> $base.D2.result.fasta
                fi
            fi
       elif  [ "$var4" != "" ];
       then
           if [ "$var5" != "" ];
           then
               echo -ne "$seq_name\n$line\n" |seqkit amplicon --quiet -F $forward_LF340_2_primer -r -245:-70 -f >> $base.D1.result.fasta
               echo -ne "$seq_name\n$line\n" |seqkit amplicon --quiet -F $forward_LF340_2_primer -R $reverse_LR3R_primer -r 90:-80 >> $base.D2.result.fasta
           else
               SeqLength=$(echo -ne "$seq_name\n$line\n" |seqkit fx2tab -l| awk '{print $3}')
               if  [ "$SeqLength" -lt 700 ];
               then
                   echo -ne "$seq_name\n$line\n" |seqkit amplicon --quiet -F $forward_LF340_2_primer -r -245:-70 -f >> $base.D1.result.fasta
                   echo -ne "$seq_name\n$line\n" |seqkit amplicon --quiet -F $forward_LF340_2_primer -r 90:340 >> $base.D2.result.fasta
                fi
            fi
        else
           if [ "$var5" != "" ];
           then
              var7=$(echo $line | grep "${forward_LR3_primer}" )
              if [ "$var7" != "" ];
              then
                  echo -ne "$seq_name\n$line\n" |seqkit amplicon --quiet -F $forward_LR3R_primer -r -330:-80 >> $base.D2.result.fasta
              fi
           fi
        fi
    fi
fi

done < $base.singline.fasta

}


export -f extract_LSU_D1D2

nohup time parallel -j 24 --eta --load 99% --noswap  extract_LSU_D1D2 ::: $(cat file_list) &



cat *.D1.result.fasta |seqkit seq -m 120 -M 200 |seqkit rmdup -n > fungi_LSU_D1.fasta
cat *.D2.result.fasta |seqkit seq -m 120 -M 300 |seqkit rmdup -n > fungi_LSU_D2.fasta
               





















