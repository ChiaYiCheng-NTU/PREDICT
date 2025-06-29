import os
import shutil
import sys
sys.path.append(f'{os.path.dirname(os.path.realpath(__file__))}/../../Kmer2Motif')
import Kmer2Motif

def GetSrcDesPath(new_folder):
    K1SrcPath = f"{new_folder}/Outputs/Kmer2Motif_inputs/K1/K1.kmer"
    K2SrcPath = f"{new_folder}/Outputs/Kmer2Motif_inputs/K2/K2.kmer"
    MotifsFolder = f"{new_folder}/MotifList/"
    MotifSrcPath = MotifsFolder + [file for file in os.listdir(MotifsFolder) if file.endswith(".motif") or file.endswith(".fasta")][0]
    
    KmerDesPath = f"{new_folder}/../../InputData/Kmer2Motif/KmerList"
    MotifDesPath = f"{new_folder}/../../InputData/Kmer2Motif/MotifList"
    return K1SrcPath, K2SrcPath, MotifSrcPath, KmerDesPath, MotifDesPath

def main(new_folder):
    K1SrcPath, K2SrcPath, MotifSrcPath, KmerDesPath, MotifDesPath = GetSrcDesPath(new_folder)
    shutil.copy(K1SrcPath, KmerDesPath)
    shutil.copy(MotifSrcPath, MotifDesPath)
    K1_K2Mfolder, K1_OutPutDF = Kmer2Motif.main(TopKmers = 1, KeepTopMotifs = 0.5, ScoreCutOff = 0.9)

    shutil.copy(K2SrcPath, KmerDesPath)
    shutil.copy(MotifSrcPath,MotifDesPath)
    K2_K2Mfolder, K2_OutPutDF = Kmer2Motif.main(TopKmers = 1, KeepTopMotifs = 0.5, ScoreCutOff = 0.9)
    return K1_OutPutDF, K2_OutPutDF, K1_K2Mfolder, K2_K2Mfolder

