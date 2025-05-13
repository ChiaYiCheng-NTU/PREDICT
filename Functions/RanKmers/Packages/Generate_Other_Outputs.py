import ahocorasick
import glob
import pandas as pd
import numpy as np
import shutil

def RC(kmer):
    complement = str.maketrans("ACGT", "TGCA")  # DNA complement mapping
    return kmer.translate(complement)[::-1]

def GrabExistedResults(new_folder, BestCopy):
    ML_Train_Result_Path = f"{new_folder}/Kmer/ML_Train_Result"
    UsefulOutputs_Path = f"{new_folder}/UsefulOutputs"
    shutil.copy(f"{ML_Train_Result_Path}/ALLCopies_mean_score.txt", f"{UsefulOutputs_Path}/ALLCopies_mean_score.txt")
    shutil.copy(f"{ML_Train_Result_Path}/AllCopies_FeatureImportance.txt", f"{UsefulOutputs_Path}/AllCopies_FeatureImportance.txt")
    shutil.copy(f"{ML_Train_Result_Path}/BestModel_Copy{BestCopy}.joblib", f"{UsefulOutputs_Path}/BestModel_Copy{BestCopy}.joblib")


def Generate_GeneTable(new_folder, BestCopy):
    def GetFeatures():
        UsedTrainDF_Path = f"{new_folder}/Kmer/Train/Train_X_{int(BestCopy):02}_df.csv"
        with open(UsedTrainDF_Path, "r") as f:
            Features_List = [line.strip() for line in f.readline().split("\t") if line != ""]
        deRepeat_Features_List = list({}.fromkeys(Features_List).keys())
        return deRepeat_Features_List ## format: [Kmer01, Kmer02, Kmer03, ...]
    
    def GetTPGenes(BestCopy):
        ## Combining TPTrainFastas and TPTestFastas of Fold{BestCopy} == Whole TPFastas
        TPTrainFastas = glob.glob(f"{new_folder}/Kmer/Train/*{int(BestCopy):02}.tp.fa")[0]
        TPTestFastas = glob.glob(f"{new_folder}/Kmer/Test/*{int(BestCopy):02}.tp.fa")[0]
        with open(TPTrainFastas)as f:
            TPlines = [line.strip() for line in f.readlines()]
        with open(TPTestFastas)as f:
            TPlines.extend([line.strip() for line in f.readlines()])

        TPFastas_Dict = {}
        for index, TPline in enumerate(TPlines):
            if index % 2 == 0:
                TPID = TPline.replace(">", "").split(" ")[0]
                TPFastas_Dict[TPID] = [TPlines[index+1], RC(TPlines[index+1])]
        return TPFastas_Dict ## format: {GeneID01: [NonTemplateSequence01, TemplateSequence01],
                             ##          GeneID02: [NonTemplateSequence02, TemplateSequence02], ...}
        
    def GetTNGenes(BestCopy):
        ## Combining TNTrainFastas and TNTestFastas of Fold{BestCopy} == Whole TNFastas
        TNTrainFastas = glob.glob(f"{new_folder}/Kmer/Train/*{int(BestCopy):02}.tn.fa")[0]
        TNTestFastas = glob.glob(f"{new_folder}/Kmer/Test/*{int(BestCopy):02}.tn.fa")[0]
        with open(TNTrainFastas)as f:
            TNlines = [line.strip() for line in f.readlines()]
        with open(TNTestFastas)as f:
            TNlines.extend([line.strip() for line in f.readlines()])

        TNFastas_Dict = {}
        for index, TNline in enumerate(TNlines):
            if index % 2 == 0:
                TNID = TNline.replace(">", "").split(" ")[0]
                TNFastas_Dict[TNID] = [TNlines[index+1], RC(TNlines[index+1])]
        return TNFastas_Dict ## format: {GeneID01: [NonTemplateSequence01, TemplateSequence01],
                             ##          GeneID02: [NonTemplateSequence02, TemplateSequence02], ...}

    def build_automaton(kmers):
        ## Build an Aho-Corasick automaton for Kmers and their reverse complements.
        ## Handles potential collision Kmer.
        kmers = list(set([line.strip("nt_").strip("t_") for line in kmers if line != ""]))
        auto = ahocorasick.Automaton()
        for kmer in kmers:
            auto.add_word(kmer, kmer)
        auto.make_automaton()
        return auto

    def KmerCount_inSeq(DoubleStrand, Kmer_list, Auto):
        seq = DoubleStrand[0]
        rc_seq = DoubleStrand[1]
        KmerCount_inSeq_Dict = {kmer:0 for kmer in Kmer_list if "nt_" in kmer}
        KmerCount_inRCSeq_Dict = {kmer:0 for kmer in Kmer_list if not "nt_" in kmer}
        for _, kmer in Auto.iter(seq):
            kmer = "nt_" + kmer
            KmerCount_inSeq_Dict[kmer] += 1
        for _, kmer in Auto.iter(rc_seq):
            kmer = "t_" + kmer
            KmerCount_inRCSeq_Dict[kmer] += 1

        KmerCount_inSeq_list = list(KmerCount_inSeq_Dict.values()) + (list(KmerCount_inRCSeq_Dict.values()))
        return KmerCount_inSeq_list

    def MakeDF(ToMakeDF_List, rowNames, GeneType):
        colName = ToMakeDF_List[0]
        data = ToMakeDF_List[1:]
        GeneXKmer_DF = pd.DataFrame(data, columns=colName, index=rowNames)
        GeneXKmer_DF.to_csv(f"{new_folder}/UsefulOutputs/{GeneType}GeneXKmer_Table.tsv", sep="\t")
        return GeneXKmer_DF
    
    Features_List = GetFeatures()
    TPFastas_Dict = GetTPGenes(BestCopy)
    TNFastas_Dict = GetTNGenes(BestCopy)
    Auto = build_automaton(Features_List)
    ToMakeTPDF_List = [Features_List.copy()]
    ToMakeTNDF_List = [Features_List.copy()]

    for TPID in TPFastas_Dict:
        ToMakeTPDF_List.append(KmerCount_inSeq(TPFastas_Dict[TPID], Features_List, Auto))
    TPGeneXKmer_DF = MakeDF(ToMakeTPDF_List, TPFastas_Dict.keys(), "TP")

    for TPID in TNFastas_Dict:
        ToMakeTNDF_List.append(KmerCount_inSeq(TNFastas_Dict[TPID], Features_List, Auto))
    TNGeneXKmer_DF = MakeDF(ToMakeTNDF_List, TNFastas_Dict.keys(), "TN")

    return TPGeneXKmer_DF, TNGeneXKmer_DF

def Generate_KmerXKmerTable(new_folder, GeneXKmer_DF=pd.DataFrame):
    kmer_matrix = GeneXKmer_DF.to_numpy()
    T_matrix = kmer_matrix.T
    matrix = kmer_matrix.astype(bool)
    T_matrix = T_matrix.astype(int)
    matrix = matrix.astype(int)
    result_matrix = np.dot(T_matrix, matrix).astype(int)
    KmerXKmer_DF = pd.DataFrame(result_matrix, index=GeneXKmer_DF.columns, columns=GeneXKmer_DF.columns)
    KmerXKmer_DF.to_csv(f"{new_folder}/UsefulOutputs/KmerXKmer_Table.tsv", sep="\t")

    return KmerXKmer_DF

def Generate_KmerShortTable(new_folder):
    UsefulOutputs_Path = f"{new_folder}/UsefulOutputs"
    with open(f"{UsefulOutputs_Path}/AllCopies_FeatureImportance.txt", "r") as f:
        Kmer_KmerID_List = ["\t".join(line.split("\t")[3:5]) + "\n" for line in f.readlines()[1:]]

    with open(f"{UsefulOutputs_Path}/Kmer_ShortTable.tsv", "w") as f:
        f.write("Kmer\tKmer_ID\n")
        f.writelines(Kmer_KmerID_List)

def main(new_folder, BestCopy):
    GrabExistedResults(new_folder, BestCopy)
    TPGeneXKmer_DF, TNGeneXKmer_DF = Generate_GeneTable(new_folder, BestCopy)
    KmerXKmer_DF = Generate_KmerXKmerTable(new_folder, TPGeneXKmer_DF)
    Generate_KmerShortTable(new_folder)









