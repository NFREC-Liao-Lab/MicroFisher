#!/bin/bash -l
#SBATCH --job-name=LSU_database_generation_preprocess_vsearch      # Job name
#SBATCH --mail-type=all            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=wanghaihua@ufl.edu       # Where to send mail	
#SBATCH --ntasks=10                      # Number of MPI ranks
#SBATCH --cpus-per-task=9               # Number of cores per MPI rank 
#SBATCH --nodes=1                       # Number of nodes
#SBATCH --ntasks-per-node=1             # How many tasks on each node
#SBATCH --ntasks-per-socket=1           # How many tasks on each CPU or socket
#SBATCH --mem-per-cpu=7G             # Memory per core
#SBATCH --time=30-00:00:00                 # Time limit hrs:min:sec
#SBATCH --output=LSU_database_generation_preprocess_vsearch_%j.log     # Standard output and error log

pwd; hostname; date


CPU="18"


#vsearch==v2.3.4
module load vsearch
###########################################################
##ITS
#################################################################

vsearch --derep_fulllength fungi_rRNA.ITS1.fasta \
        -output fungi_rRNA.ITS1_uniq.fasta \
        --id 0.99999 \
        --centroids \
        --xsize \
        --threads 24 \
        --minseqlength 120 \
        --fasta_width 0
       
        
vsearch --derep_fulllength fungi_rRNA.ITS2.fasta \
        -output fungi_rRNA.ITS2_uniq.fasta \
        --id 0.99999 \
        --centroids \
        --xsize \
        --threads 24 \
        --minseqlength 120 \
        --maxseqlength 300 \
        --fasta_width 0       
#vsearch --cluster_fast SILVA_138.1_LSUParc_tax_silva_LSU_DNA_filter_final.fasta --id 0.99999 --centroids  SILVA_138.1_LSUParc_tax_silva_LSU_DNA_filter_final_uniq.fasta





###https://drive5.com/usearch/manual/cmd_uchime2_ref.html
#remove Chimera 
#download usearch :https://drive5.com/downloads/usearch11.0.667_i86linux32.gz
#downlon
#uncompress and rename to usearch
chmod +x /home/sunny/tools/usearch
export PATH=/home/sunny/tools:$PATH

#usage
/home/sunny/tools/usearch  -uchime2_ref fungi_rRNA.ITS1_uniq.fasta \
                           -db /home/sunny/tools/uchime_reference_dataset_28.06.2017/ITS1_ITS2_datasets/uchime_reference_dataset_ITS1_28.06.2017.fasta \
                           -uchimeout fungi_rRNA.ITS1_uniq_chimera.txt \
                           -strand plus \
                           -notmatched fungi_rRNA.ITS1_uniq_chimera_filter.fa \
                           -mode high_confidence


/home/sunny/tools/usearch  -uchime2_ref fungi_rRNA.ITS2_uniq.fasta \
                           -db /home/sunny/tools/uchime_reference_dataset_28.06.2017/ITS1_ITS2_datasets/uchime_reference_dataset_ITS2_28.06.2017.fasta \
                           -uchimeout fungi_rRNA.ITS2_uniq_chimera.txt \
                           -strand plus \
                           -notmatched fungi_rRNA.ITS2_uniq_chimera_filter.fa \
                           -mode high_confidence



#############################################################################################
###LSU
###############################################################################################

vsearch --derep_fulllength fungi_LSU_D1.fasta \
        -output fungi_LSU_D1_uniq.fasta \
        --id 1 \
        --centroids \
        --xsize \
        --threads 24 \
        --minseqlength 120 \
        --fasta_width 0
       
        
vsearch --derep_fulllength fungi_LSU_D2.fasta \
        -output fungi_LSU_D2_uniq.fasta \
        --id 1 \
        --centroids \
        --xsize \
        --threads 24 \
        --minseqlength 120 \
        --maxseqlength 300 \
        --fasta_width 0       