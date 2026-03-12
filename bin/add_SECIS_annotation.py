import pandas as pd 
import argparse

def read_gff(gff):
    df = pd.read_csv(gff, sep='\t', names=['chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])
    df = df[["chr", "start", "end"]]
    df = df.set_index('chr')
    return df

def read_df(df):
    df = pd.read_csv(df)
    df = df[['transcript_name', 'og_score_start', 'og_score_end', 'og_score_score', 'og_score_length', 're_score_start', 're_score_end', 're_score_score', 're_score_length', 'TGA_site_score']]
    df = df.set_index('transcript_name')
    return df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract the coding regions of each transcript and input to coding_potential script")
    # Add arguments for input and output directories
    parser.add_argument('--ORFsearch', type=str, required=True, help="Path to ORFsearch file")
    parser.add_argument('--all_secis', type=str, required=True, help="Path to all secis elements file")  
    parser.add_argument('--filtered_secis', type=str, required=True, help="Path to filtered secis elements file")   
    parser.add_argument('--output', type=str, required=True, help="Path to output table")
    args = parser.parse_args()
    
    all_secis = read_gff(args.all_secis)
    filtered_secis = read_gff(args.filtered_secis)
    df = read_df(args.ORFsearch)
    
    df1 = df.join(all_secis, how='left', rsuffix='_all_secis')
    df2 = df1.join(filtered_secis, how='left', rsuffix='_filtered_secis')
    df2.to_csv(args.output, index=True)
    
 