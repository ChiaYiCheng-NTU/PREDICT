import os



def get_KmerList(new_folder, TopKmers):
    KmerList_Folder = os.path.join(new_folder, "KmerList")
    KmerList_File = [file for file in sorted(os.listdir(KmerList_Folder)) if file.endswith(".kmer") or file.endswith(".tsv")][0]
    with open(os.path.join(KmerList_Folder, KmerList_File)) as f:
        lines = f.readlines()
        if "Kmer" in lines[0] or "kmer" in lines[0]:
            Kmers = [line.split("\t")[0].strip() for line in lines[1:] if line != ""]
        else :
            Kmers = [line.split("\t")[0].strip() for line in lines if line != ""]
        Kmers = Kmers[:round(len(Kmers)*TopKmers)]
    Kmers_Counts = len(Kmers)
    print(f"Selected {Kmers_Counts} Kmers (Top {TopKmers}/1 of all input Kmers)")
    return Kmers

def get_MotifList(new_folder):
    def MotifListIsFasta(Motif):
        if ">" in Motif[0]:
            New_Motif = []
            for index, line in enumerate(Motif):
                if index%2 == 0:
                    New_Motif.append([Motif[index+1], line.replace(">", "")])
            return New_Motif
        else:
            return Motif
    MotifList_Folder = os.path.join(new_folder, "MotifList")
    MotifList_File = [file for file in sorted(os.listdir(MotifList_Folder)) if file.endswith(".motif") or file.endswith(".fasta")][0]
    with open(os.path.join(MotifList_Folder, MotifList_File)) as f:
        Motif = [line.strip() for line in f.readlines() if line != ""]
    
    return MotifListIsFasta(Motif)

def main(new_folder, TopKmers=0.1):
    Kmers_List = get_KmerList(new_folder, TopKmers)
    Motif_List = get_MotifList(new_folder)
    return Kmers_List, Motif_List
