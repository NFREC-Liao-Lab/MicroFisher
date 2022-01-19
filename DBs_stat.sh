workdir=/home/microbiome/data_storage/SATA2/Fisher_test/short_DBs/DBs_stat
taxdump_dir=/home/microbiome/data_storage/SATA2/Fisher_test/short_DBs/DB_taxonomy

#obtain the taxonomy lineage for each of taxid

taxonkit lineage --data-dir $taxdump_dir all_taxid |taxonkit reformat -F -P  | cut -f 1,3  > all_taxid_taxonomy


rm -f *.2lineage *.without_lineage
DBs_stat () {
file=$(basename  $1 _ACCESSION)
   cat $1 |while read line; do
        acc2taxid=$(cat accession2taxid_combination |grep -w "${line}" )
        if [ "$acc2taxid" != "" ]; then
           taxid=$(echo $acc2taxid |cut -d" " -f 2)
           lineage=$(cat all_taxid_taxonomy |grep -w $taxid)
           if [ "$lineage" != "" ]; then
              acc2taxid2lineage=$acc2taxid" "$lineage
              echo $acc2taxid2lineage >> $file.2lineage
           else
              echo $line >> $file.without_lineage
           fi
        else
           echo $line >> $file.without_lineage
        fi
    done

}           

export -f DBs_stat


time parallel -j 4 --eta --load 99% --noswap  DBs_stat ::: $(ls *_ACCESSION)










esummary -db nuccore -id $i | xtract -pattern DocumentSummary -element Caption,TaxId









