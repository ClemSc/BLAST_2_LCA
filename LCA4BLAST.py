import argparse
import gzip
import numpy as np
import pandas as pd

def get_lineage(taxid, nodes_df, result={},
                ranks=['class', 'order','family', 'genus', 'species']):
    """
    Takes an NCBI taxid and a dataframe loaded from nodes.dmp and return the lineage
    of the taxid limited to the list of provided ranks.
    Results are a dictionnary
    """
    node = nodes_df.loc[taxid]
    if node['rank'] in ranks:
        result[node['rank']] = taxid
    if node['rank'] == 'kingdom' or node['parent_taxid'] == 1:
        return result
    else:
        return get_lineage(node['parent_taxid'],
                           nodes_df,
                           result=result,
                           ranks=ranks)

def get_lca(blast_out_df, perc_hits_threshold = 90, species_level=True):
    result_df = pd.DataFrame()
    for query in blast_out_df['qseqid'].unique():
        df = blast_out_df[blast_out_df['qseqid'] == query].copy()
        for rank in ranks:
            df_val_count = df[rank].value_counts()
            df_val_perc = (df_val_count * 100) / df.shape[0]
            if df_val_perc.iloc[0] >= perc_hits_threshold:
                accepted_taxid = df_val_perc.index[0]
                df = df[df[rank] == accepted_taxid]
                df['LCA_rank'] = rank
                df['hits_num'] = df.shape[0]
                df = df.sort_values(by=['evalue'], ascending=True).iloc[0, :].copy()
                if not species_level:
                    # For low sensitivity LCA, species level is not accepted. If any,
                    # the LCA is set to the genus.
                    df['species'] = np.NaN
                    df['LCA_rank'] = df['LCA_rank'].replace('species', 'genus')
                result_df = pd.concat([result_df, df], axis=1)
                break
            else:
                accepted_taxid = False
                df[rank] = np.NaN
    return result_df.T

def get_scientific_name(lca_df, names_df, ranks):
    for rank in ranks:
        lca_df[rank] = lca_df[[rank]].merge(names_df, how='left',
                                left_on=rank, right_on='taxid')['name'].values
    return lca_df

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
            prog='LCA4BLAST',
            description='Perform a Last Common Ancestor (LCA) estimation from Blastn search results.',
            epilog='Feel free to report issues in the GitHub, contributions welcome.')

    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-t', '--nodes')
    parser.add_argument('-n', '--names')
    parser.add_argument('-f', '--fields')
    parser.add_argument('-l', '--length')
    parser.add_argument('-L', '--low_pident')
    parser.add_argument('-H', '--high_pident')
    parser.add_argument('-p', '--p_hits')
    args = parser.parse_args()

    print('Loading NCBI taxonomy')
    nodes_df = pd.read_csv(args.nodes,
                    compression='infer',
                    header=None, sep='\t',
                    usecols=[0,2,4],
                    names=["taxid", "parent_taxid", "rank"],
                    index_col=0)

    names_df = pd.read_csv(args.names,
                    compression='infer',
                    header=None, sep='\t',
                    usecols=[0,2,6],
                    names=["taxid", "name", "type"])
    names_df = names_df[names_df['type'] == 'scientific name'][['taxid', 'name']]

    # Retrieving parameters or applying default value.
    if args.length:
        length = args.length
    else:
        length = 350

    if args.low_pident:
        ls_threshold = args.low_pident
    else:
        ls_threshold = 80

    if args.high_pident:
        hs_threshold = args.high_pident
    else:
        hs_threshold = 95

    if args.p_hits:
        perc_hits_threshold = args.p_hits
    else:
        perc_hits_threshold = 90


    print('Loading Blast output file')
    blast_out_df = pd.read_csv(args.input,
                    compression='infer',
                    header=None, sep='\t')

    blast_fields = "qseqid saccver pident qcovs length evalue bitscore staxid"
    blast_out_df.columns = blast_fields.split(" ")

    print("Apply minimum thresholds then drop duplicate hit-subjects")
    blast_out_df = blast_out_df[blast_out_df['length'] >= length].copy()
    blast_out_df = blast_out_df[blast_out_df['pident'] >= ls_threshold].copy()
    blast_out_df = blast_out_df.sort_values(by=['evalue'], ascending=True)
    blast_out_df = blast_out_df.drop_duplicates(subset=['saccver'], keep='first')

    ranks =  ['species', 'genus', 'family', 'order', 'class']

    print(f'Building lineages for each unique OTU hit using following ranks:\n{ranks}')
    set_taxids = blast_out_df["staxid"].unique()
    all_lineages = list()
    for taxid in set_taxids:
        lineage = get_lineage(taxid, nodes_df, ranks=ranks)
        lineage["staxid"] = taxid
        all_lineages.append(lineage.copy())
    set_lineages_df = pd.DataFrame.from_dict(all_lineages)
    blast_out_df = blast_out_df.merge(set_lineages_df, how='left', on="staxid")

    print('Perform high sensitivity LCA')
    hs_df = blast_out_df[blast_out_df['pident'] >= hs_threshold]
    hs_df = get_lca(hs_df, perc_hits_threshold=perc_hits_threshold, species_level=True)

    print('Perform low sensitivity LCA')
    ls_df = blast_out_df[~blast_out_df['qseqid'].isin(hs_df['qseqid'].unique())]
    ls_df = get_lca(ls_df, perc_hits_threshold=perc_hits_threshold, species_level=False)

    lca_df = pd.concat([hs_df, ls_df], axis=0)

    print('Associate scientific names')
    lca_df = get_scientific_name(lca_df, names_df, ranks)

    print(f"Writing results to {args.output}, under TSV file format")
    lca_df.to_csv(args.output, sep='\t', index=False)
