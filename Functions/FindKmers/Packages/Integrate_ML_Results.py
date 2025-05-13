import os
import numpy as np

def get_all_copies_score():
    copies = 0
    ALL_copies_MeanScore = {
        "LogisticRegression":{"AUC":0, "F1":0, "MCC":0},
        "RandomForest":{"AUC":0, "F1":0, "MCC":0},
        "GradientBoosting":{"AUC":0, "F1":0, "MCC":0},
        "SVM":{"AUC":0, "F1":0, "MCC":0}
                            }

    for file in sorted(os.listdir(Train_Result_folder)):
        if file.endswith("score.txt"):
            copies += 1
            with open(f"{Train_Result_folder}/{file}") as ml_result:
                algs_result_list = "".join(ml_result.readlines()).split("---------------------------------\n")
                for one_alg_result in algs_result_list[:-1]:
                    model_name = one_alg_result.split("\n")[0].split(": ")[1]
                    three_scores_list = [float(score.split(": ")[1]) for score in one_alg_result.split("\n")[3:6]]
                    ALL_copies_MeanScore[model_name]["AUC"] += three_scores_list[0]
                    ALL_copies_MeanScore[model_name]["F1"] += three_scores_list[1]
                    ALL_copies_MeanScore[model_name]["MCC"] += three_scores_list[2]
    for keyi in ALL_copies_MeanScore.keys():
        for keyii in ALL_copies_MeanScore[keyi].keys():
            ALL_copies_MeanScore[keyi][keyii] /= copies
    return ALL_copies_MeanScore

def write_integrated_score(ALL_copies_MeanScore, alg):
    output = f"""Model: {alg}
Copy Number: ALL_COPIES_MEAN
AUC: {ALL_copies_MeanScore[alg]["AUC"]}
F1: {ALL_copies_MeanScore[alg]["F1"]}
MCC: {ALL_copies_MeanScore[alg]["MCC"]}
"""
    with open(f"{Train_Result_folder}/ALLCopies_mean_score.txt", "w") as f:
        f.write(output)
    
def get_all_copies_FeatureImportance(alg):
    result_dic = {}
    FeatureImportance_files = sorted([file for file in sorted(os.listdir(Train_Result_folder)) if file.startswith(alg)])
    MLScores_files = sorted([file for file in sorted(os.listdir(Train_Result_folder)) if file.startswith("copy_")])
    
    for Importance_file, Score_file in zip(FeatureImportance_files, MLScores_files):
        with open(f"{Train_Result_folder}/{Score_file}") as f:
            MCC = float([line for line in f.readlines() if line.startswith("MCC")][0].split(" ")[1])

        with open(f"{Train_Result_folder}/{Importance_file}") as f:
            lines = [line for line in f.readlines()[1:] if line!=""]
            for line in lines:
                Kmer, score, percentile = (line.split("\t")[1], line.split("\t")[2], line.split("\t")[3])
                score = float(score)
                percentile = float(percentile)
                if result_dic.get(Kmer, "") == "":
                    result_dic[Kmer] = [[percentile*MCC], 0, [score], 0, 1, MCC] ## [Weighted_percentile, percentile_sd, scores_list, scores_sd, count, MCC]
                else:
                    result_dic[Kmer][0].append(percentile*MCC) ## Weighted_percentile
                    result_dic[Kmer][2].append(score) ## ML_score_List
                    result_dic[Kmer][4] += 1 ## Count
                    result_dic[Kmer][5] += MCC ## MCC
    
    for Kmer in result_dic:
        Count = result_dic[Kmer][4]
        TotalWeight = result_dic[Kmer][5]
        percentile_sd = np.std(np.array(result_dic[Kmer][0]))
        result_dic[Kmer][1] = percentile_sd
        result_dic[Kmer][0] = sum(result_dic[Kmer][0]) / TotalWeight
        score_sd = np.std(np.array(result_dic[Kmer][2]))
        result_dic[Kmer][3] = score_sd
        result_dic[Kmer][2] = sum(result_dic[Kmer][2]) / Count
    # result_dic's shape: {"Kmer_1": [Weighted_percentile_1, Weighted_Percentile_sd1, Avg_Score_1, Scores_sd1, Count_1], "Kmer_2": [Weighted_percentile_2, Weighted_Percentile_sd2, Avg_Score_2, Scores_sd2, Count_2]...}
    sorted_result = sorted(result_dic.items(), key=lambda x: x[1][0])
    return sorted_result

def write_integrated_FeatureImportance(sorted_result):
    def Generate_KmerID(Kmer):
        mapping_dic = {"A": "1", "T": "2", "C": "3", "G": "4"}
        if Kmer.startswith("nt_"):
            orientation = "NonTemplate"
            Kmer = Kmer.replace("nt_", "")
            Kmer_ID = f"nt_K{str(len(Kmer)).zfill(2)}_"
        elif Kmer.startswith("t_"):
            orientation = "Template"
            Kmer = Kmer.replace("t_", "")
            Kmer_ID = f"t_K{str(len(Kmer)).zfill(2)}_"
        for mer in Kmer:
            Kmer_ID += mapping_dic[mer]
        
        return Kmer_ID, orientation

    with open(f"{Train_Result_folder}/AllCopies_FeatureImportance.txt", "w") as file:
        file.write("Final_Rank\tWeighted_percentile\tWeighted_Percentile_sd\tKmer\tKmer_ID\tOrientation\tAverage_Score\tScore_sd\tCounts\n")
        for index, (Kmer, values) in enumerate(sorted_result):
            Kmer_ID, orientation = Generate_KmerID(Kmer)
            file.write(f"{index+1}\t{round(values[0], 2)}\t{round(values[1], 2)}\t{Kmer}\t{Kmer_ID}\t{orientation}\t{values[2]}\t{values[3]}\t{values[4]}\n")

def main(new_folder, alg="RandomForest"):
    global Train_Result_folder
    Train_Result_folder = f"{new_folder}/Kmer/ML_Train_Result"

    ALL_copies_MeanScore = get_all_copies_score()
    write_integrated_score(ALL_copies_MeanScore, alg=alg)

    sorted_result = get_all_copies_FeatureImportance(alg=alg)
    write_integrated_FeatureImportance(sorted_result)