import ahocorasick
import os
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from itertools import product
import multiprocessing

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

def Generate_Kmers(k=5):
    kmers = [''.join(p) for p in product("ATCG", repeat=k)]
    return kmers

def Get_Seqences(new_folder, fold_num):
    PosPath = sorted([f"{new_folder}/Kmer/Train/" + pos for pos in sorted(os.listdir(f"{new_folder}/Kmer/Train/")) if pos.endswith(".tp.fa")])[fold_num]
    NegPath = sorted([f"{new_folder}/Kmer/Train/" + neg for neg in sorted(os.listdir(f"{new_folder}/Kmer/Train/")) if neg.endswith(".tn.fa")])[fold_num]
    
    PosSequences = []
    with open(PosPath, "r") as PosFile:
        PosSequences = [seq.strip().upper() for seq in PosFile.readlines() if not ">" in seq]
    NegSequences = []
    with open(NegPath, "r") as NegFile:
        NegSequences = [seq.strip().upper() for seq in NegFile.readlines() if not ">" in seq]
    return PosSequences, NegSequences

def Auto_CountKmers(Auto, PosSequences, NegSequences, Kmers):
    PosKmer_counts = {kmer: [0 for PosSeqNums in PosSequences] for kmer in Kmers}
    for PosSeqNum, PosSequence in enumerate(PosSequences):
        rc_PosSequence = RC(PosSequence)
        for _, found_kmer in Auto.iter(PosSequence):
            PosKmer_counts[found_kmer][PosSeqNum] += 1
        for _, found_kmer in Auto.iter(rc_PosSequence):
            PosKmer_counts[found_kmer][PosSeqNum] += 1

    NegKmer_counts = {kmer: [0 for NegSeqNums in NegSequences] for kmer in Kmers}
    for NegSeqNum, NegSequence in enumerate(NegSequences):
        rc_NegSequence = RC(NegSequence)
        for _, found_kmer in Auto.iter(NegSequence):
            NegKmer_counts[found_kmer][NegSeqNum] += 1
        for _, found_kmer in Auto.iter(rc_NegSequence):
            NegKmer_counts[found_kmer][NegSeqNum] += 1

    return PosKmer_counts, NegKmer_counts

def WilcoxonRankSum(Kmers, PosKmer_counts, NegKmer_counts):
    KmerWRS_Dict = {}
    from scipy.stats import ranksums
    for kmer in Kmers:
        Pos_kmer_list = PosKmer_counts[kmer]
        Neg_kmer_list = NegKmer_counts[kmer]
        _, pvalue = ranksums(Pos_kmer_list, Neg_kmer_list)

        PosTotalCount = sum(Pos_kmer_list)
        NegTotalCount = sum(Neg_kmer_list)
        KmerWRS_Dict[kmer] = f"{PosTotalCount}\t{NegTotalCount}\t{pvalue}"
    return KmerWRS_Dict


def FindExistedKmers(sequences, k):
    kmers = set()
    for sequence in sequences:
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i + k]
            kmers.add(kmer)
    KmerList = list(kmers)
    outputKmers = sorted(KmerList, key=atcg_sort_key)
    return outputKmers

def DropTheLargerP(KmerWRS_Dict_Result):
    DroppedLargerP_KmerPavlue_Dict = KmerWRS_Dict_Result
    KmerWRS_Dict_Keys = KmerWRS_Dict_Result.keys()
    for K in range(11, 4, -1):
        ToDelete_Set = set()
        K_KmerWRS_Dict_Keys = [Kmer for Kmer in KmerWRS_Dict_Keys if len(Kmer)==K]
        Kp1_KmerWRS_Dict_Keys = [Kmer for Kmer in KmerWRS_Dict_Keys if len(Kmer)==K+1]
        for Kmer in K_KmerWRS_Dict_Keys:
            for Kp1mer in Kp1_KmerWRS_Dict_Keys:
                if Kmer in Kp1mer:
                    if KmerWRS_Dict_Result[Kmer].split("\t")[2] >= KmerWRS_Dict_Result[Kp1mer].split("\t")[2]:
                        ToDelete_Set.add(Kmer)
                    else:
                        ToDelete_Set.add(Kp1mer)

    for ToDeleteKmer in ToDelete_Set:
        del DroppedLargerP_KmerPavlue_Dict[ToDeleteKmer]

    return DroppedLargerP_KmerPavlue_Dict


def RC(kmer):
    complement = str.maketrans("ACGT", "TGCA")  # DNA complement mapping
    return kmer.translate(complement)[::-1]

def FDR(KmerPvalue_Dict, pThreshold):
    p_values = [float(value.split("\t")[2]) for value in KmerPvalue_Dict.values()]
    _, adjusted_p_values, _, _ = multipletests(p_values, method="fdr_bh")
    FDR_KmerPvalue_Dict = {key: "{}\t{}\t{}".format(value.split('\t')[0], value.split('\t')[1], adj_p) for (key, value), adj_p in zip(KmerPvalue_Dict.items(), adjusted_p_values) if adj_p<pThreshold}
    return FDR_KmerPvalue_Dict

def process_large_kmers(kmer_subset, PosSequences, NegSequences):
    KMer_Auto = Build_Automaton(kmer_subset)
    PosKmer_counts, NegKmer_counts = Auto_CountKmers(KMer_Auto, PosSequences, NegSequences, kmer_subset)
    temp_KmerWRS_Dict = WilcoxonRankSum(kmer_subset, PosKmer_counts, NegKmer_counts)
    return temp_KmerWRS_Dict

def atcg_sort_key(kmer):
    atcg_order = {'A': 0, 'T': 1, 'C': 2, 'G': 3}
    return [atcg_order[base] for base in kmer if base in atcg_order]

def main(new_folder, pThreshold):
    Fivemers = Generate_Kmers(k=5)
    FiveMer_Auto = Build_Automaton(Fivemers)

    for fold in range(5):
        PosSequences, NegSequences = Get_Seqences(new_folder, fold)
        PosKmer_counts, NegKmer_counts = Auto_CountKmers(FiveMer_Auto, PosSequences, NegSequences, Fivemers)

        print(f"Finding Enrich 5mers from fold{fold+1}...")
        KmerWRS_Dict = WilcoxonRankSum(Fivemers, PosKmer_counts, NegKmer_counts)

        temp_KmerWRS_Dict = KmerWRS_Dict.copy()
        KmerWRS_Dict_Result = KmerWRS_Dict.copy()
        print(f"Found {len(KmerWRS_Dict)} 5mers")

        for plusmer in range(6): #from 5mer to 11mer
            K = 6+plusmer
            print(f"Finding Enrich {K}mers from fold{fold+1}...")
            
            if K <= 8:
                Kmers = Generate_Kmers(k=K)
                print(f"Found {len(Kmers)} {K}mers")
                KMer_Auto = Build_Automaton(Kmers)
                PosKmer_counts, NegKmer_counts = Auto_CountKmers(KMer_Auto, PosSequences, NegSequences, Kmers)
                temp_KmerWRS_Dict = WilcoxonRankSum(Kmers, PosKmer_counts, NegKmer_counts)
                KmerWRS_Dict_Result.update(temp_KmerWRS_Dict)
            
            elif K > 8:
                Kmers = FindExistedKmers(PosSequences, k=K)
                print(f"Found {len(Kmers)} {K}mers")
                #MutiProcess
                num_processes = multiprocessing.cpu_count()
                chunk_size = len(Kmers) // num_processes
                kmer_chunks = [Kmers[i:i + chunk_size] for i in range(0, len(Kmers), chunk_size)]
                with multiprocessing.Pool(num_processes) as pool:
                    results = pool.starmap(process_large_kmers, [(chunk, PosSequences, NegSequences) for chunk in kmer_chunks])
                
                temp_Integrate_result_Dict = {}
                for result_Dict in results:
                    temp_Integrate_result_Dict.update(result_Dict)
                temp_Integrate_result_Dict = {Kmer: temp_Integrate_result_Dict[Kmer] for Kmer in Kmers}

                KmerWRS_Dict_Result.update(temp_Integrate_result_Dict)

                

        print("Doing FDR...")
        FDR_KmerWRS_Dict = FDR(KmerWRS_Dict_Result, pThreshold)
        print("FDR Done!")

        print("Dropping Kmer with Larger P...")
        DroppedLargerP_KmerPavlue_Dict = DropTheLargerP(FDR_KmerWRS_Dict)
        print("Dropping Kmer with Larger P Done!")

        OutPut_KmerWRS_Dict = DroppedLargerP_KmerPavlue_Dict.copy()

        

        OutputPath = f"{new_folder}/Kmer/Train/TPTN_{fold+1:02d}.pcre_WRSresults.txt"
        with open(OutputPath, "w") as OutputFile:
            OutputFile.write("Kmer\tPoscounts\tNegcounts\tP-value\n")
            for key, value in OutPut_KmerWRS_Dict.items():
                OutputFile.write(f"{key}\t{value}\n")




