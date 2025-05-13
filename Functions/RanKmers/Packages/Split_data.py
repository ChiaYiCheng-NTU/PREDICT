import os
from sklearn.model_selection import KFold

def main(new_folder, n_splits=5, seed=42):
    print(f"n_splits: {n_splits}")
    print(f"start_seed: {seed}")

    TP_and_TN_folder = f"{new_folder}/TP_and_TN"
    
    for copy_num in range(1, n_splits + 1):
        for file in sorted(os.listdir(TP_and_TN_folder)):
            if file.endswith(".tp"):
                with open(f"{TP_and_TN_folder}/{file}", "r") as old_tp:
                    old_tp_lines = old_tp.readlines()
                    old_tp_lines_n = len(old_tp_lines)

                    kf = KFold(n_splits=n_splits, shuffle=True, random_state=seed+copy_num)
                    for train_index, test_index in kf.split(old_tp_lines):
                        TP_train = [old_tp_lines[i] for i in train_index]
                        TP_test = [old_tp_lines[i] for i in test_index]

                        tp_out_middle_name = file.split(".")[0]

            if file.endswith(".tn"):
                with open(f"{TP_and_TN_folder}/{file}", "r") as old_tn:
                    old_tn_lines = old_tn.readlines()
                    old_tn_lines_n = len(old_tn_lines)

                    kf = KFold(n_splits=n_splits, shuffle=True, random_state=seed+copy_num)
                    for train_index, test_index in kf.split(old_tn_lines):
                        TN_train = [old_tn_lines[i] for i in train_index]
                        TN_test = [old_tn_lines[i] for i in test_index]

                        tn_out_middle_name = file.split(".")[0]
        
        os.makedirs(f"{new_folder}/Kmer/Test", exist_ok=True)
        with open(f"{new_folder}/Kmer/Test/TP_{tp_out_middle_name}_{copy_num:02d}.tp", "w") as Test_TP:
            Test_TP.writelines(TP_test)
        with open(f"{new_folder}/Kmer/Test/TN_{tn_out_middle_name}_{copy_num:02d}.tn", "w") as Test_TN:
            Test_TN.writelines(TN_test)

        os.makedirs(f"{new_folder}/Kmer/Train", exist_ok=True)
        with open(f"{new_folder}/Kmer/Train/TP_{tp_out_middle_name}_{copy_num:02d}.tp", "w") as Train_TP:
            Train_TP.writelines(TP_train)
        with open(f"{new_folder}/Kmer/Train/TN_{tn_out_middle_name}_{copy_num:02d}.tn", "w") as Train_TN:
            Train_TN.writelines(TN_train)

    return n_splits
