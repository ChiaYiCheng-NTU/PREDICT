import sys
import os
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
import subprocess
from Packages import Dealing_with_Folders_K2M
from Packages import Get_Inputs
from Packages import Compare_KmerToMotif

def main(Kmer, Motif, TopKmers = 0.1, KeepTopMotifs = 0.1, ScoreCutOff = 0.9):
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#
        ## 1. Copy "Kmer2Motif" folder(in "Demanded_Data" folder, old folder) to "Result" folder's new folder ##
    print("Dealing with folders...")
    new_folder = Dealing_with_Folders_K2M.main(Kmer, Motif)
    print("Dealing with folders Done!")
    print("=======================================================")
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#

    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#
        ## 2. Get "KmerList" and "MotifList" input
    print("Get_Inputs...")
    KmerList, MotifList = Get_Inputs.main(new_folder=new_folder, TopKmers=TopKmers)
    print("Get_Inputs Done!")
    print("=======================================================")
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#

    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#
        ## 3. Comparing inputs(Kmer to Motif), to determine Motifs from the Kmers
    print("Compare_Kmer_To_Motif...")
    OutPutDF = Compare_KmerToMotif.main(new_folder=new_folder, KmerList=KmerList, MotifList=MotifList, ScoreCutOff=ScoreCutOff, KeepTopMotifs=KeepTopMotifs)
    print("Compare_Kmer_To_Motif Done!")
    print("--+--+--+--+--+--+--+--+--E..N..D--+--+--+--+--+--+--+--+--")
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#
    return new_folder, OutPutDF

if __name__ == "__main__":
    main()
