import pandas as pd

def Generate_KmerID(Kmer):
        mapping_dic = {"A": "1", "T": "2", "C": "3", "G": "4"}
        if Kmer.startswith("nt_"):
            Kmer = Kmer.replace("nt_", "")
            Kmer_ID = f"nt_K{str(len(Kmer)).zfill(2)}_"
        elif Kmer.startswith("t_"):
            Kmer = Kmer.replace("t_", "")
            Kmer_ID = f"t_K{str(len(Kmer)).zfill(2)}_"
        for mer in Kmer:
            Kmer_ID += mapping_dic[mer]
        return Kmer_ID



def main(new_folder, OutputDF):
    DFindex = OutputDF.index
    K1_K1ID_List = [(KmerPair[0], Generate_KmerID(KmerPair[0])) for KmerPair in DFindex]
    K2_K2ID_List = [(KmerPair[1], Generate_KmerID(KmerPair[1])) for KmerPair in DFindex]
    K1OutputDF = pd.DataFrame(K1_K1ID_List, index=None, columns=["Kmer", "Kmer_ID"])
    K2OutputDF = pd.DataFrame(K2_K2ID_List, index=None, columns=["Kmer", "Kmer_ID"])
    
    K1OutputDF.to_csv(f"{new_folder}/Outputs/Kmer2Motif_inputs/K1/K1.kmer", index=None, sep="\t")
    K2OutputDF.to_csv(f"{new_folder}/Outputs/Kmer2Motif_inputs/K2/K2.kmer", index=None, sep="\t")

    



