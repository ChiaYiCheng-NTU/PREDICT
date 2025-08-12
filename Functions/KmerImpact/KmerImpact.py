import sys
import os
sys.path.append(f'{os.path.dirname(os.path.realpath(__file__))}/Packages')
from datetime import datetime
import Dealing_with_Folders
import Gff_to_Coord
import Coord_to_Fasta
import Split_data
import TPTN_to_Fasta
import Find_Enriches
import Balance_Data
import Enriches_to_DF
import DF_to_ML
import Integrate_ML_Results
import Generate_Other_Outputs
import warnings
warnings.filterwarnings("ignore")


def main(gff, genome, tp, tn, features = "gene", pThreshold = 0.01, alg = "RandomForest", up_stream = 1000, down_stream = 500):
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#
        ## 1. Copy "KmerImpact" folder(in "Demanded_Data" folder, old folder) to "Result" folder's new folder ##
    print("Dealing_with_folders...")
    new_folder = Dealing_with_Folders.main(gff, genome, tp, tn)
    print("Dealing_with_folders Done!")
    print("=======================================================")
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#

    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#
        ## 2. From Gff file extract sequnce coord(start nucleotide number & end nucleotide number) ##
    print("Gff_to_Coord...")
    Gff_to_Coord.main(new_folder, features=features, up_stream=up_stream, down_stream=down_stream)
    print ("Gff_to_Coord Done!")
    print("=======================================================")
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#

    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#
        ## 3. Use Coord to get exact sequence
    print("Coord_to_Fasta...")
    Coord_to_Fasta.main(new_folder)
    print("Coord_to_Fasta Done!")
    print("=======================================================")
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#

    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#
        ## 4. Split train/test set (70/30)
    print("Split_data...")
    copies = Split_data.main(new_folder, n_splits=5, seed=42)
    print("Split_data Done!")
    print("=======================================================")
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#

    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#
        ## 5. Generate train/test_set's fasta file
    print("TPTN_to_Fasta...")
    TPTN_to_Fasta.main(new_folder, For_BedData=False)
    print("TPTN_to_Fasta Done!")
    print("=======================================================")
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#

    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#
        ## 6. Find enriches from train_fasta
    print("Find_Enriches...")
    Find_Enriches.main(new_folder, pThreshold=pThreshold)
    print("Find_Enriches Done!")
    print("=======================================================")
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#

    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#
        ## 7. Balance TP/TN Data
    print("Balance_Data...")
    Balance_Data.main(new_folder)
    print("Balance_Data Done!")
    print("=======================================================")
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#

    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#
        ## 8. Get sequences of Balanced data
    print("BalancedData_to_Fasta...")
    TPTN_to_Fasta.main(new_folder, For_BedData=True)
    print("BalancedData_to_Fasta Done!")
    print("=======================================================")
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#

    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#
        ## 9. Make DataFrame
    print("Enriches_to_DF...")
    Enriches_to_DF.main(new_folder, copies)
    print("Enriches_to_DF Done!")
    print("=======================================================")
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#

    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#
        ## 10. Fit Machine Learning
    print("DF_to_ML...")
    BestCopy = DF_to_ML.main(new_folder, copies, alg=alg)
    print("DF_to_ML Done!")
    print("=======================================================")
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#

    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#
        ## 11. Integrate ML results
    print("Integrate_ML_results...")
    Integrate_ML_Results.main(new_folder, alg=alg)
    print("Integrate_ML_results Done!")
    print("=======================================================")
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#

    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#
        ## 12. Generate other useful outputs
    print("Generate other useful outputs...")
    Generate_Other_Outputs.main(new_folder, BestCopy)
    print("Generate other useful outputs Done!")
    print("--+--+--+--+--+--+--+--+--E..N..D--+--+--+--+--+--+--+--+--")
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#
    return new_folder

if __name__ == "__main__":
    StartTime = datetime.now()

    print(f"Folder: {main()}")

    EndTime = datetime.now()
    TimeDiff = int((EndTime - StartTime).total_seconds())
    print(f"Time cost: {TimeDiff//3600}hr {(TimeDiff % 3600)//60}min {TimeDiff%60}sec")
