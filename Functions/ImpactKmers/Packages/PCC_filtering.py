import pandas as pd
import os



def main(new_folder, rThreshold):
    Train_folder = os.path.join(new_folder, "Kmer/Train")
    Test_folder = os.path.join(new_folder, "Kmer/Test")
    Test_X_dfs_list = sorted([os.path.join(Test_folder, file) for file in os.listdir(Test_folder) if file.startswith("Test_X")])
    Train_X_dfs_list = sorted([os.path.join(Train_folder, file) for file in os.listdir(Train_folder) if file.startswith("Train_X")])
    
    for index, X_df in enumerate(Train_X_dfs_list):
        print(f"processing fold{index+1}...")
        df = pd.read_csv(X_df, sep="\t")
        corr_df = df.corr().abs()
        
        features = df.columns
        features_count = len(features)

        filtered_features = set()
        for i in range(features_count):
            feature_1 = features[i]
            if feature_1 not in filtered_features:
                for j in range(i + 1, features_count-1):
                    feature_2 = features[j]
                    if feature_2 not in filtered_features:
                        pearson_cor = corr_df.iloc[i, j]
                        if pearson_cor > rThreshold:
                            filtered_features.add(feature_2)
        
        feature_keep = [feature for feature in features if feature not in filtered_features]
        print(f"filtered {len(filtered_features)} features")
        filterd_df = df[feature_keep]
        filterd_df.to_csv(X_df, sep="\t")

        test_df = pd.read_csv(Test_X_dfs_list[index], sep="\t")
        filtered_test_df = test_df[feature_keep]
        filtered_test_df.to_csv(Test_X_dfs_list[index], sep="\t")
