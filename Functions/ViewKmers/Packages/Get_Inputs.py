import os
import pandas

def get_KmerList(new_folder, TopKmers):
    KmerList_Folder = os.path.join(new_folder, "KmerList")
    KmerList_File = [file for file in sorted(os.listdir(KmerList_Folder)) if file.endswith(".kmers") or file.endswith(".tsv")][0]
    with open(os.path.join(KmerList_Folder, KmerList_File)) as f:
        Kmers = [line.split("\t")[0] for line in f.readlines()[1:] if line != ""]
        Kmers = Kmers[:round(len(Kmers)*TopKmers)]
    Kmers_Counts = len(Kmers)
    print(f"Selected {Kmers_Counts} Kmers (Top {TopKmers} of all input Kmers)")
    return Kmers


def main(new_folder):
    TPGeneXKmerDF, TNGeneXKmerDF = get_GeneXKmerDF(new_folder)
    return TPGeneXKmerDF, TNGeneXKmerDF

