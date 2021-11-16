#The scripts for an example pipeline usage and comparisons with other packages

#The example data collection
#Download the RefSeq 28S_fungal_sequences from NCBI
wget -c https://ftp.ncbi.nih.gov/blast/db/28S_fungal_sequences.tar.gz
tar -zxvf 28S_fungal_sequences.tar.gz

#Extract sample sequences from the 28S_fungal_sequences

for num in 10 50 100 ;do
   for i in $(seq 1 5);do
       mkdir $Test_Fungi_RefSeq/simulating_$num""seqs.${i}
       mkdir $Test_Fungi_RefSeq/simulating_$num""seqs.${i}/splite_seq
       
       seed=$[RANDOM%100+10]
       cat $Test_Fungi_RefSeq/RefSeq_Fungi/28S_fungal_sequences.fa  |seqkit sample -n $num -s $seed > $Test_Fungi_RefSeq/simulating_$num""seqs.${i}/28S_fungal_sequences_sample_$num""seqs.${i}.fasta
       cat $Test_Fungi_RefSeq/simulating_$num""seqs.${i}/28S_fungal_sequences_sample_$num""seqs.${i}.fasta |seqkit split2 -s 1 --out-dir $Test_Fungi_RefSeq/simulating_$num""seqs.${i}/splite_seq
       
       for seqs in $(ls $Test_Fungi_RefSeq/simulating_$num""seqs.${i}/splite_seq/*.fasta); do 
           id=$(grep ">" $seqs)
           genus=$(echo $id |cut -d" "  -f 2)
           species=$(echo $id |cut -d" " -f 3)
           name=$genus"_"$species
           mv $seqs $Test_Fungi_RefSeq/simulating_$num""seqs.${i}/splite_seq/$name.fasta
       done
       
       for acc in $(cat $Test_Fungi_RefSeq/simulating_$num""seqs.${i}/28S_fungal_sequences_sample_$num""seqs.${i}.fasta |seqkit seq -i -n);do
           accession=$acc
           taxid=$(echo $acc |esummary -db nuccore  | xtract -pattern DocumentSummary -element TaxId)
           linage=$(echo $taxid | taxonkit lineage --data-dir $DB_taxonomy |taxonkit reformat -P --data-dir $DB_taxonomy --delimiter ";" --format "{k};{p};{c};{o};{f};{g};{s}" | cut -f 1,3)
           species=$accession'   '$linage
           echo $species >> $Test_Fungi_RefSeq/simulating_$num""seqs.${i}/28S_fungal_sequences_sample_$num""seqs.${i}.taxonomy.txt
       done
           
       
       time iss generate --cpus 24 --quiet  --compress \
            --draft $Test_Fungi_RefSeq/simulating_$num""seqs.${i}/splite_seq/*.fasta \
            --model NovaSeq \
            --n_reads 20M \
            --output $Test_Fungi_RefSeq/simulating_$num""seqs.${i}/simulating_$num""seqs.${i}.short_read
  
       time centrifuge  -x $LSU_D1D2_DBs_new/LSU_D1_fisher_new \
           -1  $Test_Fungi_RefSeq/simulating_$num""seqs.${i}/simulating_$num""seqs.${i}.short_read_R1.fastq.gz  \
           -2  $Test_Fungi_RefSeq/simulating_$num""seqs.${i}/simulating_$num""seqs.${i}.short_read_R2.fastq.gz  \
           -S $Test_Fungi_RefSeq/simulating_$num""seqs.${i}/simulating_$num""seqs.${i}.LSU_DB1.short_read.result \
           --report-file $Test_Fungi_RefSeq/simulating_$num""seqs.${i}/simulating_$num""seqs.${i}.LSU_DB1.novaseq_reads.reprot.tsv -p 24 -k 1  --min-hitlen 140
        time centrifuge-kreport -x $LSU_D1D2_DBs_new/LSU_D1_fisher_new \
                   $Test_Fungi_RefSeq/simulating_$num""seqs.${i}/simulating_$num""seqs.${i}.LSU_DB1.novaseq_reads.reprot.tsv \
                   > $Test_Fungi_RefSeq/simulating_$num""seqs.${i}/simulating_$num""seqs.${i}.LSU_DB1.novaseq_reads.kreprot.tsv
   
                     
       time centrifuge  -x $LSU_D1D2_DBs_new/LSU_D2_fisher_new \
           -1  $Test_Fungi_RefSeq/simulating_$num""seqs.${i}/simulating_$num""seqs.${i}.short_read_R1.fastq.gz  \
           -2  $Test_Fungi_RefSeq/simulating_$num""seqs.${i}/simulating_$num""seqs.${i}.short_read_R2.fastq.gz  \
           -S $Test_Fungi_RefSeq/simulating_$num""seqs.${i}/simulating_$num""seqs.${i}.LSU_DB2.short_read.result \
           --report-file $Test_Fungi_RefSeq/simulating_$num""seqs.${i}/simulating_$num""seqs.${i}.LSU_DB2.novaseq_reads.reprot.tsv -p 24 -k 1  --min-hitlen 140
        time centrifuge-kreport -x $LSU_D1D2_DBs_new/LSU_D1_fisher_new \
                   $Test_Fungi_RefSeq/simulating_$num""seqs.${i}/simulating_$num""seqs.${i}.LSU_DB2.novaseq_reads.reprot.tsv \
                   > $Test_Fungi_RefSeq/simulating_$num""seqs.${i}/simulating_$num""seqs.${i}.LSU_DB2.novaseq_reads.kreprot.tsv

done
