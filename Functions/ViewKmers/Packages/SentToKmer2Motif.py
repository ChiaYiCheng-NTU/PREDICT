import os
import shutil
import sys
sys.path.append(f'{os.path.dirname(os.path.realpath(__file__))}/../../Kmer2Motif')
import Kmer2Motif

def GetSrcDesPath(new_folder):
    KmerFolder = f"{new_folder}/KmerList/"
    MotifsFolder = f"{new_folder}/MotifList/"
    KmerSrcPath = KmerFolder + [file for file in os.listdir(KmerFolder) if file.endswith(".kmer") or file.endswith(".tsv")][0]
    MotifSrcPath = MotifsFolder + [file for file in os.listdir(MotifsFolder) if file.endswith(".motif") or file.endswith(".fasta")][0]
    
    KmerDesPath = f"{new_folder}/../../InputData/Kmer2Motif/KmerList"
    MotifDesPath = f"{new_folder}/../../InputData/Kmer2Motif/MotifList"
    return KmerSrcPath, MotifSrcPath, KmerDesPath, MotifDesPath

def Make_K2M_Dict(K2M_OutPutDF):
    keys = K2M_OutPutDF["KmerSeq"]
    values = K2M_OutPutDF["MotifInfo"] if K2M_OutPutDF["MotifInfo"][0] != None else K2M_OutPutDF["MotifSeq"]
    K2M_Dict = {key: [] for key in keys}
    for key, value in zip(keys, values):
        K2M_Dict[key].append(value)
    return K2M_Dict

def main(new_folder):
    KmerSrcPath, MotifSrcPath, KmerDesPath, MotifDesPath = GetSrcDesPath(new_folder)
    _, K2M_OutPutDF = Kmer2Motif.main(KmerSrcPath, MotifSrcPath, TopKmers = 1, KeepTopMotifs = 0.5, ScoreCutOff = 0.9)

    K2M_Dict = Make_K2M_Dict(K2M_OutPutDF)
    return K2M_Dict

