from ete3 import NCBITaxa
ncbi = NCBITaxa()
# ncbi.update_taxonomy_database()


DESIRED_RANKS = ["family", "genus", "species"]
FULL_RANKS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]


def get_desired_taxa_ranks(taxid, desired_ranks=None):

    if desired_ranks is None:
        desired_ranks = FULL_RANKS
    # match = re.match(PATTERN_Genus_species, desc)
    # if match:
    #     tax_g_s = f"{match.group(1)} {match.group(2)}"
    # taxid = ncbi.get_name_translator([tax_g_s])[tax_g_s][0]
    lineage = ncbi.get_lineage(taxid)
    names = ncbi.get_taxid_translator(lineage)
    lineage2ranks = ncbi.get_rank(names)
    ranks2lineage = dict((rank, taxid)
                         for (taxid, rank) in lineage2ranks.items())
    ranks_dict = {'{}'.format(rank): ranks2lineage.get(
                  rank, '<not present>') for rank in desired_ranks}
    # ranks_filtered = dict(filter(lambda kv: kv[1]!="<not present>", ranks_dict.items()))
    ranks_filtered = {k: v for k, v in ranks_dict.items()
                      if v != "<not present>"}

    id_to_name = ncbi.get_taxid_translator(ranks_filtered.values())
    des_rank_names = {}
    for k, v in ranks_filtered.items():
        des_rank_names[k] = (v, id_to_name[v])
    # des_rank_names = list(ncbi.get_taxid_translator(ranks_dict{"order"}).values())
    # des_rank_names = list(ncbi.get_taxid_translator(ranks_dict.values()).values())
    # s = f"{match.group(1)}_{match.group(2)}"
    # s = f"O_{des_rank_names['order']}_F_{des_rank_names['family']}"
    return des_rank_names
