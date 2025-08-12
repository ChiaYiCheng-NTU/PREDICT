import os
import random

def Balance_train_folder(new_folder):
    Positives = sorted([f"{new_folder}/Kmer/Train/" + pos for pos in sorted(os.listdir(f"{new_folder}/Kmer/Train/")) if pos.endswith(".tp")])
    Negatives = sorted([f"{new_folder}/Kmer/Train/" + neg for neg in sorted(os.listdir(f"{new_folder}/Kmer/Train/")) if neg.endswith(".tn")])
    
    for index, value in enumerate(Positives):
        pos = value
        neg = Negatives[index]
        print(f"TP path: {pos}")
        print(f"TN path: {neg}")
        with open(pos, "r") as tp:
            tp_lines = tp.readlines()
        with open(neg, "r") as tn:
            tn_lines = tn.readlines()

        TP_balanced = ""
        TN_balanced = ""
        min_samples = min(len(tp_lines), len(tn_lines))
        if len(tp_lines) > min_samples:
            print(f"{pos.split('/')[-1]} is downsizing to {min_samples} samples")
            random.seed(42+index)
            TP_balanced = random.sample(tp_lines, min_samples)
            TN_balanced = tn_lines
        elif len(tn_lines) > min_samples:
            print(f"{neg.split('/')[-1]} is downsizing to {min_samples} samples")
            random.seed(42+index)
            TN_balanced = random.sample(tn_lines, min_samples)
            TP_balanced = tp_lines
        else:
            print("Data is already balanced, nothing needs to do.")
            TN_balanced = tn_lines
            TP_balanced = tp_lines
        
        with open(f"{pos.replace('.tp', '.Btp')}", "w") as Btp:
            Btp.writelines(TP_balanced)
        with open(f"{neg.replace('.tn', '.Btn')}", "w") as Btn:
            Btn.writelines(TN_balanced)
    
def main(new_folder):
    Balance_train_folder(new_folder=new_folder)
