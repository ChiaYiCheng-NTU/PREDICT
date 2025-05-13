import numpy as np
import pandas as pd
import random
from scipy.stats import wilcoxon
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.preprocessing import MinMaxScaler

# def get_Pairs(GeneXKmerDF):
#     Pairs = {}
#     cols = GeneXKmerDF.columns
#     for K1Index, K1 in enumerate(cols):
#         K1_counts = GeneXKmerDF.loc[:,K1]
#         for K2 in cols[K1Index:]:
#             K2_counts = GeneXKmerDF.loc[:,K2]
#             K1K2_IntersctionMin = tuple(np.minimum(K1_counts, K2_counts))
#             Pairs[f"{K1}/{K2}"] = K1K2_IntersctionMin
#     return Pairs ## Shape: {K1/K2: (G1CountMin, G2CountMin, ...)}
def RC(kmer):
    complement = str.maketrans("ACGT", "TGCA")  # DNA complement mapping
    return kmer.translate(complement)[::-1]

def get_Pairs(GeneXKmerDF, Data, TPCSThreshold=0.75): ##Fake one, this is to get Single cols in DFs
    Pairs = {}
    TPCSvalue_Dict = {}
    TPKmerTotalcount_DF = pd.DataFrame(columns=("Kmer1_Count", "Kmer2_Count", "2Kmer_CO-Count"))
    cols = GeneXKmerDF.columns
    for K1Index, K1 in enumerate(cols):
        K1_counts = GeneXKmerDF.loc[:,K1]
        Scaled_K1_counts = MinMaxScaler().fit_transform(np.array([K1_counts]).T).T
        for K2 in cols[K1Index+1:]:
            K2_counts = GeneXKmerDF.loc[:,K2]
            Scaled_K2_counts = MinMaxScaler().fit_transform(np.array([K2_counts]).T).T
            if K1.split("_")[0] != K2.split("_")[0] and K1.split("_")[1] == RC(K2.split("_")[1]):
                continue
            else:
                if Data == "TP":
                    CS = cosine_similarity(Scaled_K1_counts, Scaled_K2_counts)[0][0]
                    if CS > TPCSThreshold:
                        FlattenScaled_K1_counts = Scaled_K1_counts.flatten()
                        FlattenScaled_K2_counts = Scaled_K2_counts.flatten()
                        K1K2_IntersctionMin = pd.Series((FlattenScaled_K1_counts + FlattenScaled_K2_counts) / 2)
                        Pairs[f"{K1}/{K2}"] = K1K2_IntersctionMin
                        TPCSvalue_Dict[f"{K1}/{K2}"] = CS
                        ToAppendSeries = pd.Series({
                            'Kmer1_Count': sum(K1_counts),
                            'Kmer2_Count': sum(K2_counts),
                            '2Kmer_CO-Count': sum(pd.Series((K1_counts + K2_counts) / 2))
                        }, name=(K1, K2))
                        TPKmerTotalcount_DF = pd.concat([TPKmerTotalcount_DF, ToAppendSeries.to_frame().T])
                elif Data == "TN":
                    FlattenScaled_K1_counts = Scaled_K1_counts.flatten()
                    FlattenScaled_K2_counts = Scaled_K2_counts.flatten()
                    K1K2_IntersctionMin = pd.Series((FlattenScaled_K1_counts + FlattenScaled_K2_counts) / 2)
                    Pairs[f"{K1}/{K2}"] = K1K2_IntersctionMin
                else:
                    raise ValueError("Please choose correct Data type (TP or TN).")
    return Pairs, TPCSvalue_Dict, TPKmerTotalcount_DF ## Shape: {K1/K2: (G1CountMin, G2CountMin, ...)}

def Wilcoxon_TPTNPairs(TPPairs, TNPairs):
    def OverSampling(TPPair, TNPairm, seed=42):
        num_TP = len(TPPair)
        num_TN = len(TNPair)
        random.seed(seed)
        if num_TP > num_TN:
            TN_Resampled = random.choices(TNPair, k=num_TP)
            TP_Resampled = TPPair
        elif num_TP < num_TN:
            TP_Resampled = random.choices(TPPair, k=num_TN)
            TN_Resampled = TNPair
        else:
            TN_Resampled = TNPair
            TP_Resampled = TPPair
        return TP_Resampled, TN_Resampled
    
    KmerPairsPvalue_Dict = {pair: 1 for pair in TPPairs.keys()}
    for Pair in TPPairs:
        TPPair = TPPairs[Pair]
        TNPair = TNPairs[Pair]
        TP_Resampled, TN_Resampled = OverSampling(TPPair, TNPair)
        _, pvalue = wilcoxon(TP_Resampled, TN_Resampled)
        KmerPairsPvalue_Dict[Pair] = pvalue
    return KmerPairsPvalue_Dict ## Shape: {K1/K2: Pvalue, K1/K3: Pvalue, ..., K2/K3: Pvalue, ...}

# def Wilcoxon_TPTNPairs(Pairs): ##Fake one, this is to Wilcoxon two cols in one DF.
#     KmerPairsPvalue_Dict = {pair: 1 for pair in Pairs.keys()}
#     for Pair in Pairs:
#         Pair1 = Pairs[Pair][0]
#         Pair2 = Pairs[Pair][1]
#         _, pvalue = ranksums(Pair1, Pair2)
#         if Pair.split("/")[0] == Pair.split("/")[1]:
#             print(pvalue)
#         KmerPairsPvalue_Dict[Pair] = pvalue
#     return KmerPairsPvalue_Dict ## Shape: {K1/K2: Pvalue, K1/K3: Pvalue, ..., K2/K3: Pvalue, ...}

def CosineSimilarity_TPTNPairs(TPPairs, TNPairs):
    def OverSampling(TPPair, TNPairm, seed=42):
        num_TP = len(TPPair)
        num_TN = len(TNPair)
        random.seed(seed)
        if num_TP > num_TN:
            TN_Resampled = [random.choices(TNPair, k=num_TP)]
            TP_Resampled = [TPPair]
        elif num_TP < num_TN:
            TP_Resampled = [random.choices(TPPair, k=num_TN)]
            TN_Resampled = [TNPair]
        else:
            TN_Resampled = [TNPair]
            TP_Resampled = [TPPair]
        return TP_Resampled, TN_Resampled

    KmerCSvalue_Dict = {pair: 0 for pair in TPPairs.keys()}
    for Pair in TPPairs:
        TPPair = TPPairs[Pair]
        TNPair = TNPairs[Pair]
        TP_Resampled, TN_Resampled = OverSampling(TPPair, TNPair)
        CS_Value = cosine_similarity(TP_Resampled, TN_Resampled)
        KmerCSvalue_Dict[Pair] = CS_Value[0][0]
    return KmerCSvalue_Dict

def Make_DF(new_folder, KmerPairsPvalue_Dict, TPGeneXKmerDF, DF_Type):
    cols = TPGeneXKmerDF.columns
    rows = TPGeneXKmerDF.columns ## Yes, columns. BC TPGeneXKmerDF.index are genes, here we want is KmerXKmer Table.
    DF = pd.DataFrame(np.nan, index=rows, columns=cols)
    for pair in KmerPairsPvalue_Dict:
        K1 = pair.split("/")[0]
        K2 = pair.split("/")[1]
        DF.loc[K1, K2] = KmerPairsPvalue_Dict[pair]
    DF.to_csv(f"{new_folder}/Outputs/KmerXKmer{DF_Type}Table.tsv", sep="\t")
    return DF

def get_Pavlues(PvalueDF, pThreshold, Comparison):
    cols = PvalueDF.columns
    Pvalues = []
    for K1Index, K1 in enumerate(cols):
        for K2 in cols[K1Index+1:]:
            Pvalues.append(PvalueDF.loc[K1, K2])
    if Comparison == "Greater":
        Filtered_Pvalues = [p for p in Pvalues if p > pThreshold]
    elif Comparison == "Less":
        Filtered_Pvalues = [p for p in Pvalues if p < pThreshold]
    else:
        raise ValueError("arg: Comparison must be ('Greater') or 'Less'")
    return Filtered_Pvalues

def main(new_folder, TPGeneXKmerDF, TNGeneXKmerDF, TPCSThreshold, TPTNCSThreshold):
    TPPairs, TPCSvalue_Dict, TPKmerTotalcount_DF = get_Pairs(TPGeneXKmerDF, "TP", TPCSThreshold)
    print(f"Got {len(TPPairs)} TPPairs!")
    print("Applying Wilcoxon on each pair...")
    # inTPKmerPairsPvalue_Dict = Wilcoxon_TPTNPairs(TPPairs)
    # TPPvalueDF = Make_DF(new_folder, inTPKmerPairsPvalue_Dict, TPGeneXKmerDF, "inTPPvalue")
    # Filtered_Pvalues = get_Pavlues(TPPvalueDF, pThreshold, Comparison)
    TNPairs, _, _ = get_Pairs(TNGeneXKmerDF, "TN")
    KmerPairsPvalue_Dict = Wilcoxon_TPTNPairs(TPPairs, TNPairs)
    KmerCSvalue_Dict = CosineSimilarity_TPTNPairs(TPPairs, TNPairs)
    # print(f"Got {len(TNPairs)} TNPairs!")
    # print("Applying CosineSimilarity on each pair...")
    # KmerCSvalue_Dict = CosineSimilarity_TPTNPairs(TPPairs, TNPairs)
    TPCSDF = Make_DF(new_folder, TPCSvalue_Dict, TPGeneXKmerDF, "TPCosineSimilarity")
    Pvalue_DF = Make_DF(new_folder, KmerPairsPvalue_Dict, TPGeneXKmerDF, "Pvalue")
    TPTNCS_DF = Make_DF(new_folder, KmerCSvalue_Dict, TPGeneXKmerDF, "TPTNCosineSimilarity")

    StackedTPCSDF = TPCSDF.stack()
    StackedTPCSDF.name = "TP Cosine Similarity"
    StackedTPTNCSDF = TPTNCS_DF.stack()
    StackedTPTNCSDF.name = "TPTN Cosine Similarity"
    Qualified_Pairs = StackedTPTNCSDF[StackedTPTNCSDF < TPTNCSThreshold]
    OutputDF = pd.concat([Qualified_Pairs, StackedTPCSDF], axis=1, join='inner')\
        .round(4)\
        .sort_values(by="TPTN Cosine Similarity", ascending=True)
    OutputDF = pd.concat([OutputDF, TPKmerTotalcount_DF], axis=1, join='inner')
    OutputDF = OutputDF.rename_axis(index=["Kmer1", "Kmer2"])
    OutputDF.to_csv(f"{new_folder}/Outputs/CoocurredKmerPairs.tsv", sep="\t")
    

    return Pvalue_DF, OutputDF
    

