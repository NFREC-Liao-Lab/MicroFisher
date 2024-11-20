# The script for build of contrifuge database index.

#!/bin/bash

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

#http://www.ccb.jhu.edu/software/centrifuge/manual.shtml
centrifuge-build --conversion-table $DB_taxonomy/accession2taxid --taxonomy-tree $DB_taxonomy/nodes.dmp --name-table $DB_taxonomy/names.dmp $ITS_DBs/fungi_rRNA.ITS1_uniq_chimera_filter.fa $ITS_DBs/ITS1_fisher
centrifuge-build --conversion-table $DB_taxonomy/accession2taxid --taxonomy-tree $DB_taxonomy/nodes.dmp --name-table $DB_taxonomy/names.dmp $ITS_DBs/fungi_rRNA.ITS2_uniq_chimera_filter.fa $ITS_DBs/ITS2_fisher

centrifuge-build --conversion-table $LSU_D1D2_DBs_new/accession2taxid --taxonomy-tree $DB_taxonomy/nodes.dmp --name-table $DB_taxonomy/names.dmp $LSU_D1D2_DBs_new/fungi_LSU_D1_uniq.fasta $LSU_D1D2_DBs_new/LSU_D1_fisher_new
centrifuge-build --conversion-table $LSU_D1D2_DBs_new/accession2taxid --taxonomy-tree $DB_taxonomy/nodes.dmp --name-table $DB_taxonomy/names.dmp $LSU_D1D2_DBs_new/fungi_LSU_D2_uniq.fasta $LSU_D1D2_DBs_new/LSU_D2_fisher_new

centrifuge-build --conversion-table $DB_combination/accession2taxid_combination --taxonomy-tree $DB_taxonomy/nodes.dmp --name-table $DB_taxonomy/names.dmp $DB_combination/combination_rmdup.fa $DB_combination_fisher







########################################################################
centrifuge-build --conversion-table $DB_taxonomy/accession2taxid --taxonomy-tree $DB_taxonomy/nodes.dmp --name-table $DB_taxonomy/names.dmp $ITS_DBs/fungi_rRNA.ITS1_uniq_chimera_filter_subsample.fa $ITS_DBs/ITS1_fisher
centrifuge-build --conversion-table $DB_taxonomy/accession2taxid --taxonomy-tree $DB_taxonomy/nodes.dmp --name-table $DB_taxonomy/names.dmp $ITS_DBs/fungi_rRNA.ITS2_uniq_chimera_filter_subsample.fa $ITS_DBs/ITS2_fisher

centrifuge-build --conversion-table $LSU_D1D2_DBs_new/accession2taxid --taxonomy-tree $DB_taxonomy/nodes.dmp --name-table $DB_taxonomy/names.dmp $LSU_D1D2_DBs_new/fungi_LSU_D1_uniq_subsample.fasta $LSU_D1D2_DBs_new/LSU_D1_fisher_new
centrifuge-build --conversion-table $LSU_D1D2_DBs_new/accession2taxid --taxonomy-tree $DB_taxonomy/nodes.dmp --name-table $DB_taxonomy/names.dmp $LSU_D1D2_DBs_new/fungi_LSU_D2_uniq_subsample.fasta $LSU_D1D2_DBs_new/LSU_D2_fisher_new

centrifuge-build --conversion-table $DB_combination/accession2taxid_combination --taxonomy-tree $DB_taxonomy/nodes.dmp --name-table $DB_taxonomy/names.dmp $DB_combination/combination_subsample_rmdup.fa $DB_combination_fisher
