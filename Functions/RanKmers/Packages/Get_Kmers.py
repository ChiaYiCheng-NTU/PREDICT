import os

def get_KmerList(new_folder):
    KmerList_Folder = os.path.join(new_folder, "KmerList")
    KmerList_File = [file for file in sorted(os.listdir(KmerList_Folder)) if file.endswith(".kmer") or file.endswith(".tsv")][0]
    with open(os.path.join(KmerList_Folder, KmerList_File)) as f:
        Kmers = []
        for line in f:
            line = line.strip()
            if not line or "Kmer" in line:
                continue
            first_col = line.split("\t")[0]
            Kmers.append(first_col.split("_")[1] if "_" in first_col else first_col)
    Kmers_Counts = len(Kmers)
    print(f"Ranking {Kmers_Counts} Kmers...")
    return Kmers

def WriteKmerLists(new_folder, Kmers):
    for fold in range(5):
        OutputPath = f"{new_folder}/Kmer/Train/Kmers_{fold+1:02d}.txt"
        with open(OutputPath, "w") as OutputFile:
            OutputFile.write("Kmer\tPoscounts\tNegcounts\tP-value\n")
            for kmer in Kmers:
                OutputFile.write(f"{kmer}\n")

def main(new_folder):
    Kmers = get_KmerList(new_folder)
    print(Kmers)
    WriteKmerLists(new_folder, Kmers)
