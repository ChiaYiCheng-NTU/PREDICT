import os

def Complement(kmer):
    complement = str.maketrans("ACGT", "TGCA")  # DNA complement mapping
    return kmer.translate(complement)

def get_Seqences(new_folder):
    GeneFasPath = sorted([f"{new_folder}/GeneList/" + gene for gene in sorted(os.listdir(f"{new_folder}/GeneList")) if gene.endswith(".fa")])[0]
    print(GeneFasPath)
    with open(GeneFasPath, "r") as GeneFasFile:
        lines = GeneFasFile.readlines()
        Names = [line.strip(">").strip().split(" ")[0] for line in lines[0::2]]
        Seqs = [line.strip().upper() for line in lines[1::2]]
        RC_Seqs = [Complement(seq) for seq in Seqs]
        GeneSequences = {Name: (Seq, RC_Seq) for Name, Seq, RC_Seq in zip(Names, Seqs, RC_Seqs)}
    return GeneSequences

def get_KmerList(new_folder):
    KmerList_Folder = os.path.join(new_folder, "KmerList")
    KmerList_File = [file for file in sorted(os.listdir(KmerList_Folder)) if file.endswith(".kmer") or file.endswith(".tsv")][0]
    try:
        with open(os.path.join(KmerList_Folder, KmerList_File)) as f:
            Kmers = [line.split("\t")[0].split("_")[1] for line in f.readlines()[1:] if line != ""]
    except IndexError:
        with open(os.path.join(KmerList_Folder, KmerList_File)) as f:
            Kmers = [line.strip() for line in f.readlines() if line != ""]
    return Kmers

def main(new_folder):
    SeqencesDict = get_Seqences(new_folder)
    KmerList = get_KmerList(new_folder)
    return SeqencesDict, KmerList