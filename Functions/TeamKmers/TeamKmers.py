import sys
import os
sys.path.append(f"{os.path.dirname(os.path.realpath(__file__))}/Packages")
import Dealing_with_Folders
import Get_Inputs_GXK
import CalculatePvalue
import Generate_KmersMotif_input
import SentToKmer2Motif
import CalculateMotifCoScore 

def main(Motif=True, TPCSThreshold=0.75, TPTNCSThreshold=0.3):
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#
        # 1. Copy "TeamKmers" folder(in "Demanded_Data" folder, old folder) to "Result" folder's new folder ##
    print("Dealing_with_folders...")
    new_folder = Dealing_with_Folders.main()
    print("Dealing_with_folders Done!")
    print("=======================================================")
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#

    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#
        ## 2. Get "KmerXKmerTable" input
    print("Get_Inputs...")
    TPGeneXKmerDF, TNGeneXKmerDF = Get_Inputs_GXK.main(new_folder)
    print("Get_Inputs Done!")
    print("=======================================================")
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#

    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#
        ## 3. Calculate the Cosine similarity
    print("Evaluate the co-occurrence...")
    Pvalue_DF, CS_DF = CalculatePvalue.main(new_folder, TPGeneXKmerDF, TNGeneXKmerDF,
                                                            TPCSThreshold=TPCSThreshold, TPTNCSThreshold=TPTNCSThreshold)
    print("Evaluate the co-occurrence Done!")
    print("=======================================================")
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#
    
    if Motif == True:
        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#
            ## 4. Generate Kmers2Motif input
        print("Generate Kmers2Motif input...")
        Generate_KmersMotif_input.main(new_folder, CS_DF)
        print("Generate Kmers2Motif input Done!")
        print("=======================================================")
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#

        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#
            ## 5. Sent to Kmer2Motif
        print("Sent to Kmer2Motif input...")
        K1_OutPutDF, K2_OutPutDF, K1_K2Mfolder, K2_K2Mfolder = SentToKmer2Motif.main(new_folder)
        print("Sent to Kmer2Motif Done!")
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#

        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#
            ## 6. Calculate Motifs' co-occurrence score (max = 1)
        print("Calculate Motifs' co-occurrence score...")
        result = CalculateMotifCoScore.main(new_folder, CS_DF, K1_OutPutDF, K2_OutPutDF)
        print("Calculate Motifs' co-occurrence score Done!")
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#
        
if __name__ == "__main__":
    main()