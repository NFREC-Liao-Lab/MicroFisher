#!/bin/bash -l
#SBATCH --job-name=LSU_database_generation_preprocess_taxonomy      # Job name
#SBATCH --mail-type=all            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=wanghaihua@ufl.edu       # Where to send mail	
#SBATCH --ntasks=1000000000000000                      # Number of MPI ranks
#SBATCH --cpus-per-task=1               # Number of cores per MPI rank 
#SBATCH --nodes=1                       # Number of nodes
#SBATCH --ntasks-per-node=1             # How many tasks on each node
#SBATCH --ntasks-per-socket=1           # How many tasks on each CPU or socket
#SBATCH --mem-per-cpu=7G             # Memory per core
#SBATCH --time=1-00:00:00                 # Time limit hrs:min:sec
#SBATCH --output=LSU_database_generation_preprocess_taxonomy_%j.log     # Standard output and error log

pwd; hostname; date


CPU="18"
#module load entrez-direct
#module load ncbitax2lin
#module load taxonkit 
module load ucsc/20210803
module load seqkit/0.10.2
module load conda 
conda activate RNASeq


# conda install -c bioconda entrez-direct
# pip install -U ncbitax2lin /conda install -c bioconda taxonkit 

#This is the code for generation of mulit- short database D1D2 region


###############################################
#prepare the sequence taxonomy
#link the accession number to taxid and taxname

#1. download database
# mkdir taxdump
# wget -c https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
# mv taxdump.tar.gz taxdump
# tar -zxvf taxdump/taxdump.tar.gz -C taxdump

#2. extract the accession number
grep -r ">" SILVA_138.1_LSUParc_tax_silva_LSU_DNA_filter_final_uniq.fasta|cut -d">" -f 2  > accession.number

wc -l accession.number > number_of_seq


#3. assign the accession number to taxid
touch accession2taxid
for i in $(cat accession.number)
do
esummary -db nuccore -id $i | xtract -pattern DocumentSummary -element Caption,TaxId >> accession2taxid
done


#list unsigned taxids
for i in $(cat accession.number)
do
  var1=$(echo | grep -r  "$i" accession2taxid)
  if [[ "$var1" = "" ]];
  then
    echo $i >> taxid_unassigned
    esummary -db nuccore -id $i | xtract -pattern DocumentSummary -element Caption,TaxId >> taxid_unassigned
   fi
done




# cat accession | epost -db nuccore | esummary -db nuccore | xtract -pattern DocumentSummary -element Caption,TaxId > accession2taxid


#4. extract the taxid
awk '{print $2}' accession2taxid |grep -v -w "0" |sort -u > taxids
wc -l taxids

#5. link taxid to lineage
taxonkit lineage --data-dir taxdump taxids > lineage.txt

#if only fungi included,shoud extract the only fungal lineage
grep -r "Fungi" >> lineage_fungi.txt

#6. formate the lineage
#taxonkit reformat -P --data-dir ../taxdump lineage.txt \
#--delimiter ";" \
#--fill-miss-rank \
#--format "{k};{p};{c};{o};{f};{g};{s}" \
#--miss-rank-repl "0" \
#--miss-rank-repl-prefix "unclassified " \
#| cut -f 1,3 > taxid2lineage

taxonkit reformat -P  lineage_fungi.txt \
--data-dir taxdump \
--delimiter ";" \
--format "{k};{p};{c};{o};{f};{g};{s}" \
| cut -f 1,3  > taxid2lineage_fungi



#7. refine the lineage (remove the lineage without "Order" level )
 # sed "/o__;f__;/d" taxid2lineage_fungi > taxid2lineage_assigned

#8. extract assigned taxid
awk '{print $1}' taxid2lineage_fungi > taxid_assigned_fungi

#9. extract assigned accession2taxid
touch accession2taxid_fungi_assigned

for i in $(cat taxid_assigned_fungi)
do
grep -r "$i" accession2taxid >> accession2taxid_fungi_assigned
done

awk '{print $1}' accession2taxid_fungi_assigned |cut -d":" -f 2 |sort -u > accession_assigned_fungi


########################################################################
#extract fasta with new accession ids

faSomeRecords SILVA_138.1_LSUParc_tax_silva_LSU_DNA_filter_final_uniq.fasta \
                                       accession_assigned_fungi \
                                       SILVA_138.1_LSUParc_tax_silva_LSU_DNA_filter_final_uniq_assigned.fasta


seqkit seq SILVA_138.1_LSUParc_tax_silva_LSU_DNA_filter_final_uniq_assigned.fasta -w 0 > SILVA_138.1_LSUParc_tax_silva_LSU_DNA_filter_final_uniq_assigned_nmdump_removed.fasta 

mv SILVA_138.1_LSUParc_tax_silva_LSU_DNA_filter_final_uniq_assigned.fasta ../result_dir
mv SILVA_138.1_LSUParc_tax_silva_LSU_DNA_filter_final_uniq_assigned_nmdump_removed.fasta ../result_dir
mv accession2taxid_fungi_assigned ../result_dir
mv taxid2lineage_fungi ../result_dir


















########################################################################################################################################
#precess the fungal accession file
#########################################################################################################################################
wget -c https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
gunzip nucl_gb.accession2taxid.gz

awk '{print $2,$3}' nucl_gb.accession2taxid > accession2taxid

download the fungi taxonomy ids：fungi_taxonomy_id



#######################################################
###ITS accession
##
split -l 3000 ITS_accession ITS_accession_splited


acc2taxid () {

i=$1

for acc in $(cat $i)
do
grep -w "$acc" accession2taxid >> $i.accession2taxid
done

}

export -f acc2taxid

time parallel -j 24 --eta --load 99% --noswap  acc2taxid ::: $(ls ITS_accession_splited*) 



###########################################################
###LSU accession
###

split -l 2000 fungal_LSU_acc fungal_LSU_acc_splited

acc2taxid () {

i=$1

for acc in $(cat $i)
do
grep -w "$acc" accession2taxid >> $i.accession2taxid_lsu
done

}

export -f acc2taxid

time parallel -j 24 --eta --load 99% --noswap  acc2taxid ::: $(ls fungal_LSU_acc_splited*) 





