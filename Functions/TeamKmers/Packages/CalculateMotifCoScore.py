import pandas as pd

def main(new_folder, CS_DF, K1_OutPutDF, K2_OutPutDF):
    merged1 = CS_DF.reset_index().merge(K1_OutPutDF, left_on="Kmer1", right_on="KmerSeq", how="left").drop(columns=["KmerSeq"])
    merged2 = merged1.merge(K2_OutPutDF, left_on="Kmer2", right_on="KmerSeq", how="left").drop(columns=["KmerSeq"])
    merged2["MCo Score"] = merged2["TP Cosine Similarity"] * merged2["Score_x"] * merged2["Score_y"]
    print(merged2.shape, merged2.head())
    MCo_DF = merged2[["MotifInfo_x",
                      "MotifInfo_y",
                      "MCo Score",
                      "Kmer1_Count",
                      "Kmer2_Count",
                      "2Kmer_CO-Count"]] \
        .rename(columns={"MotifInfo_x": "Motif1",
                         "MotifInfo_y": "Motif2",
                         "Kmer1_Count": "Motif1_Count",
                         "Kmer2_Count": "Motif2_Count",
                         "2Kmer_CO-Count": "2Motifs_CO-Count"})
    print(MCo_DF.shape, MCo_DF.head())

    MCo_DF["Motif1_Count"] = MCo_DF["Motif1_Count"].astype(int)
    MCo_DF["Motif2_Count"] = MCo_DF["Motif2_Count"].astype(int)
    MCo_DF["2Motifs_CO-Count"] = MCo_DF["2Motifs_CO-Count"].astype(int)
    
    MCo_DF = MCo_DF.groupby(["Motif1", "Motif2"], as_index=False).sum()
    MCo_DF = MCo_DF.sort_values(by="MCo Score", ascending=False)
    MCo_DF.to_csv(f"{new_folder}/Outputs/MotifCoocurrenceScore.tsv", sep="\t", index=None)


    return MCo_DF