import ahocorasick
import os
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

def Build_Automaton(kmers):
    """
    Build an Aho-Corasick automaton for Kmers and their reverse complements.
    Handles potential collision Kmer.
    """
    auto = ahocorasick.Automaton()
    for kmer in kmers:
        # Add to automaton
        auto.add_word(kmer, kmer)
    auto.make_automaton()
    return auto

def Get_Kmers():
    Kmers = []
    with open("./Functions/FindKmers/Packages/5mer.txt", "r") as file:
        Kmers = [kmer.strip() for kmer in file.readlines()]
    return Kmers

def Get_Seqences(new_folder, fold_num):
    PosPath = sorted([f"{new_folder}/Kmer/Train/" + pos for pos in sorted(os.listdir(f"{new_folder}/Kmer/Train/")) if pos.endswith(".tp.fa")])[fold_num]
    NegPath = sorted([f"{new_folder}/Kmer/Train/" + neg for neg in sorted(os.listdir(f"{new_folder}/Kmer/Train/")) if neg.endswith(".tn.fa")])[fold_num]
    
    PosSeqences = []
    with open(PosPath, "r") as PosFile:
        PosSeqences = [seq.strip().upper() for seq in PosFile.readlines() if not ">" in seq]
    NegSeqences = []
    with open(NegPath, "r") as NegFile:
        NegSeqences = [seq.strip().upper() for seq in NegFile.readlines() if not ">" in seq]
    return PosSeqences, NegSeqences

def Auto_CountKmers(Auto, PosSeqences, NegSeqences, Kmers):
    PosKmer_counts = {kmer: 0 for kmer in Kmers}
    for PosSequence in PosSeqences:
        temp_PosKmer_counts = {kmer: 0 for kmer in Kmers}
        for _, (found_kmer, foundRC_kmer) in Auto.iter(PosSequence): #Mod
            temp_PosKmer_counts[found_kmer] = 1
            temp_PosKmer_counts[foundRC_kmer] = 1 #Mod
        for key in PosKmer_counts:
            PosKmer_counts[key] += temp_PosKmer_counts[key]

    NegKmer_counts = {kmer: 0 for kmer in Kmers}
    for NegSequence in NegSeqences:
        temp_NegKmer_counts = {kmer: 0 for kmer in Kmers}
        for _, (found_kmer, foundRC_kmer) in Auto.iter(NegSequence): #Mod
            temp_NegKmer_counts[found_kmer] = 1
            temp_NegKmer_counts[foundRC_kmer] = 1 #Mod
        for key in NegKmer_counts:
            NegKmer_counts[key] += temp_NegKmer_counts[key]
    return PosKmer_counts, NegKmer_counts

# def FisherExactTest(Kmers, PosKmer_counts, NegKmer_counts, Pos_num, Neg_num, pThreshold):
#     KmerFET_Dict = {}
#     for kmer in Kmers:
#         TP = PosKmer_counts[kmer]            #Positive Examples with kmer present
#         FP = NegKmer_counts[kmer]            #Negative Examples with kmer present
#         TN = Neg_num-NegKmer_counts[kmer]    #Negative Examples without kmer
#         FN = Pos_num-PosKmer_counts[kmer]    #Positive Examples without kmer
#         _,pvalue = fisher_exact([[TP,FN],[FP,TN]],alternative='greater')

#         if pvalue < pThreshold:
#             KmerFET_Dict[kmer] = f"{TP}\t{FP}\t{pvalue}"
#     return KmerFET_Dict

def FisherExactTest(Kmers):
    KmerFET_Dict = {}
    for Kmer in Kmers:
        KmerFET_Dict.update({Kmer: "0\t0\t0"})

    return KmerFET_Dict

def PlusOneMer(KmerFET_Dict, PosSeqences, NegSeqences, pThreshold):

    K_plus1_merFET_Dict = {}

    for Kmer in KmerFET_Dict :
        temp_K_plus1_mers = {f"A{Kmer}":0, f"T{Kmer}":0, f"C{Kmer}":0, f"G{Kmer}":0}
        K_plus1_merFET_Dict.update(temp_K_plus1_mers)
    
    K_plus1_mers = K_plus1_merFET_Dict.keys()
    # K_plus1_mer_Auto = build_automaton_with_rc(K_plus1_mers)
    # PosKp1mer_counts, NegKp1mer_counts = Auto_CountKmers(Auto=K_plus1_mer_Auto, PosSeqences=PosSeqences, NegSeqences=NegSeqences, Kmers=K_plus1_merFET_Dict.keys())

    Pos_num = len(PosSeqences)
    Neg_num = len(NegSeqences)
    KmerSeriesFET_Dict = FisherExactTest(K_plus1_merFET_Dict.keys())#, PosKp1mer_counts, NegKp1mer_counts, Pos_num, Neg_num, pThreshold=pThreshold)
    return KmerSeriesFET_Dict

def DropTheLargerP(KmerFET_Dict_Result):
    DroppedLargerP_KmerPavlue_Dict = KmerFET_Dict_Result
    KmerFET_Dict_Keys = KmerFET_Dict_Result.keys()
    for K in range(11, 4, -1):
        ToDelete_Set = set()
        K_KmerFET_Dict_Keys = [Kmer for Kmer in KmerFET_Dict_Keys if len(Kmer)==K]
        Kp1_KmerFET_Dict_Keys = [Kmer for Kmer in KmerFET_Dict_Keys if len(Kmer)==K+1]
        for Kmer in K_KmerFET_Dict_Keys:
            for Kp1mer in Kp1_KmerFET_Dict_Keys:
                if Kmer in Kp1mer:
                    if KmerFET_Dict_Result[Kmer].split("\t")[2] >= KmerFET_Dict_Result[Kp1mer].split("\t")[2]:
                        ToDelete_Set.add(Kmer)
                    else:
                        ToDelete_Set.add(Kp1mer)

    for ToDeleteKmer in ToDelete_Set:
        del DroppedLargerP_KmerPavlue_Dict[ToDeleteKmer]
    # for ToDeleteKmer in ToDelete_Set:
    #     if ToDeleteKmer in KmerFET_Dict_Result.keys():
    #         OutPut_KmerFET_Dict.update({ToDeleteKmer: KmerFET_Dict_Result[ToDeleteKmer]})
    return DroppedLargerP_KmerPavlue_Dict

# def filter_significant_features_from_dict(feature_dict, r_threshold=0.95):
#     import numpy as np

#     significant_feature_Kmers = list(feature_dict.keys())
#     size_K_mers = len(significant_feature_Kmers)

#     corr_matrix = np.random.rand(size_K_mers, size_K_mers)
#     corr_matrix = (corr_matrix + corr_matrix.T) / 2
#     np.fill_diagonal(corr_matrix, 1)

#     filtered_features = set()

#     for i in range(size_K_mers):
#         feature_1 = significant_feature_Kmers[i]
#         if feature_1 not in filtered_features:
#             for j in range(i + 1, size_K_mers):
#                 feature_2 = significant_feature_Kmers[j]
#                 if feature_2 not in filtered_features:
#                     pearson_cor = corr_matrix[i, j]
#                     if pearson_cor >= r_threshold:
#                         filtered_features.add(feature_2)

#     remaining_features = {f: feature_dict[f] for f in significant_feature_Kmers if f not in filtered_features}
#     return remaining_features

def RC(kmer):
    complement = str.maketrans("ACGT", "TGCA")  # DNA complement mapping
    return kmer.translate(complement)[::-1]

def FDR(KmerPvalue_Dict, pThreshold):
    p_values = [float(value.split("\t")[2]) for value in KmerPvalue_Dict.values()]
    _, adjusted_p_values, _, _ = multipletests(p_values, method="fdr_bh")
    FDR_KmerPvalue_Dict = {key: "{}\t{}\t{}".format(value.split('\t')[0], value.split('\t')[1], adj_p) for (key, value), adj_p in zip(KmerPvalue_Dict.items(), adjusted_p_values) if adj_p<pThreshold}
    return FDR_KmerPvalue_Dict

def main(new_folder, pThreshold):
    Kmers = Get_Kmers()
    # FiveMer_Auto = Build_Automaton(Kmers)

    for fold in range(5):
        PosSeqences, NegSeqences = Get_Seqences(new_folder, fold)
        # PosKmer_counts, NegKmer_counts = Auto_CountKmers(FiveMer_Auto, PosSeqences, NegSeqences, Kmers)

        Pos_num = len(PosSeqences)
        Neg_num = len(NegSeqences)
        print(f"Finding Enrich 5mers from fold{fold+1}...")
        KmerFET_Dict = FisherExactTest(Kmers)#, PosKmer_counts, NegKmer_counts, Pos_num, Neg_num, pThreshold)

        temp_KmerFET_Dict = KmerFET_Dict
        KmerFET_Dict_Result = KmerFET_Dict
        print(f"Found {len(KmerFET_Dict)} 5mers")
        for plusmer in range(1): #from 5mer to 11mer
            print(f"Finding Enrich {6+plusmer}mers from fold{fold+1}...")

            K_plus1_merFET_Dict = PlusOneMer(temp_KmerFET_Dict, PosSeqences, NegSeqences, pThreshold)
            temp_KmerFET_Dict = K_plus1_merFET_Dict
            KmerFET_Dict_Result.update(temp_KmerFET_Dict)
            
            print(f"Found {len(K_plus1_merFET_Dict)} {6+plusmer}mers")

        print("Doing FDR...")
        FDR_KmerFET_Dict = FDR(KmerFET_Dict_Result, pThreshold)
        print("FDR Done!")

        # print("Dropping Kmer with Larger P...")
        # DroppedLargerP_KmerPavlue_Dict = DropTheLargerP(KmerFET_Dict_Result)
        # print("Dropping Kmer with Larger P Done!")

        OutPut_KmerFET_Dict = FDR_KmerFET_Dict

        OutputPath = f"{new_folder}/Kmer/Train/TPTN_{fold+1:02d}.pcre_FETresults.txt"
        with open(OutputPath, "w") as OutputFile:
            OutputFile.write("Kmer\tPoscounts\tNegcounts\tP-value\n")
            for key, value in OutPut_KmerFET_Dict.items():
                OutputFile.write(f"{key}\t{value}\n")
        


## Mod
## RC




