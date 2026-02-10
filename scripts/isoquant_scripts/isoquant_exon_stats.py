import os
import sys
import pandas as pd
from collections import defaultdict




def calculate_reads_per_exon(tbl):

    df = pd.read_table(tbl, sep='\t', header=None, usecols=[0, 1, 3, 4, 7], 
                   names=["read_id","chr","tx_id", "gene_id", "exons"])
    # df.set_index("read_id", inplace=True)
    print(df.head())

    # Initialize a dictionary to store exon IDs and corresponding read IDs
    exon_dict = defaultdict(list)

    # Iterate through each row in the DataFrame
    for index, row in df.iterrows():
        read_id = row['read_id']
        chr_id = row['chr']
        exons = row['exons'].split(',')
        
        # Concatenate chromosome ID with each exon block to form unique exon IDs
        for exon in exons:
            exon_id = f"{chr_id}_{exon}"
            exon_dict[exon_id].append(read_id)

    exon_count = pd.DataFrame({'exon_id': list(exon_dict.keys()), 
                               'read_count': [len(read_ids) for read_ids in exon_dict.values()]})
    print(exon_count.head())

    # Create a DataFrame with exon_id to list of read_ids
    # exon_to_reads_df = pd.DataFrame({'exon_id': list(exon_dict.keys()), 
    #                                  'read_ids': list(exon_dict.values())})
    # print(exon_to_reads_df.head())

    # get a table maps an exon_id -> a read_id 
    # exon_df = pd.DataFrame([(exon, read_id) for exon, read_ids in exon_dict.items() ], columns=['exon_id', 'read_id'])
    # print(exon_df.head())

    return exon_count


def main(tbl, output):
    # tbl_name = os.path.basename(tbl).split(".")[-1]
    # output_file = f'{tbl_name}_readsPerExon'

    df = calculate_reads_per_exon(tbl)
    df.to_csv(output, sep='\t', header=None, index=None)

    print(f"stats created: {output}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script <read_assignments.tsv> <output>")
        sys.exit(1)

    tbl = sys.argv[1]
    output = sys.argv[2]
    if not os.path.isfile(tbl,):
        print(f"Error: {tbl} is not a valid file.")
        sys.exit(1)

    main(tbl, output)

