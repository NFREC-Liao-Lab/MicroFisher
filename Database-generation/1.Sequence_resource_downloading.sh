#this is for sequence resource downloding 
#!/bin/bash -l
#SBATCH --job-name=LSU_database_generation_download      # Job name
#SBATCH --mail-type=all            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=wanghaihua@ufl.edu       # Where to send mail	
#SBATCH --ntasks=1                      # Number of MPI ranks
#SBATCH --cpus-per-task=9               # Number of cores per MPI rank 
#SBATCH --nodes=1                       # Number of nodes
#SBATCH --ntasks-per-node=1             # How many tasks on each node
#SBATCH --ntasks-per-socket=1           # How many tasks on each CPU or socket
#SBATCH --mem-per-cpu=7G             # Memory per core
#SBATCH --time=30-00:00:00                 # Time limit hrs:min:sec
#SBATCH --output=LSU_database_generation_download_%j.log     # Standard output and error log

pwd; hostname; date


CPU="18"

#This is the code for generation of mulit- short database D1D2 region

###############################################################

#1.download the fungal nucleotide sequences from NCBI (INSDC and Refseq)
workdir=
mkdir $workdir/fungal_nucleotide

#download ncbi accession numbers(taxid4751:fungi)
esearch -db "nucleotide" -query "txid4751[Organism]"|efetch -format fasta > $workdir/fungal_nucleotide/fungi_nucleotide 
esearch -db "nucleotide" -query "txid4751[Organism]"|efetch -format acc > $workdir/fungal_nucleotide/fungal_nucleotide.acc








#Download the raw LSU database #


#The database was downloaded from SILVA collected un-filtered database 
#As well, the sequences could be collected from International Nucleotide Sequence Database Collaboration (INSDC) databases

#workdir="/home/wanghaihua/all_data/database_generation"
#mkdir $work_dir/raw_database
#cd $work_dir/raw_database

#wget -c https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_LSUParc_tax_silva.fasta.gz $work_dir/raw_database
#gunzip $work_dir/raw_database/SILVA_138.1_LSUParc_tax_silva.fasta.gz

#RNA to DNA
#seqkit seq --threads 10 --rna2dna SILVA_138.1_LSUParc_tax_silva.fasta > SILVA_138.1_LSUParc_tax_silva_DNA.fasta

