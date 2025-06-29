from sklearn.model_selection import GridSearchCV, PredefinedSplit
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.svm import SVC
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import make_scorer, roc_auc_score, f1_score, accuracy_score, matthews_corrcoef
import numpy as np
import pandas as pd
import os
import json
import warnings
from joblib import dump
warnings.filterwarnings("ignore")


def read_Train_test_files(copy_num):
    for file in sorted(os.listdir(Train_TPTN_folder)):
        if file.endswith(f"X_{copy_num:02d}_df.csv"):
            with open(f"{Train_TPTN_folder}/{file}") as X_Train_file:
                X_Train_list = [item.strip().split("\t") for item in X_Train_file.readlines()]
            col_names = X_Train_list[0]
            values = [row for row in X_Train_list[1:]]
            X_Train_df = pd.DataFrame(values, columns=['SeqName'] + col_names).set_index('SeqName')
            X_Train_df = X_Train_df.astype(int)
        elif file.endswith(f"y_{copy_num:02d}_df.csv"):
            with open(f"{Train_TPTN_folder}/{file}") as y_Train_file:
                y_Train_list = [item.strip().split("\t") for item in y_Train_file.readlines()]
            col_names = y_Train_list[0]
            values = [row for row in y_Train_list[1:]]
            y_Train_df = pd.DataFrame(values, columns=['SeqName'] + col_names).set_index('SeqName')
            y_Train_df = y_Train_df.astype(int).to_numpy().ravel()
    for file in sorted(os.listdir(Test_TPTN_folder)):
        if file.endswith(f"X_{copy_num:02d}_df.csv"):
            with open(f"{Test_TPTN_folder}/{file}") as X_Test_file:
                X_Test_list = [item.strip().split("\t") for item in X_Test_file.readlines()]
            col_names = X_Test_list[0]
            values = [row for row in X_Test_list[1:]]
            X_Test_df = pd.DataFrame(values, columns=['SeqName'] + col_names).set_index('SeqName')
            X_Test_df = X_Test_df.astype(int)
        elif file.endswith(f"y_{copy_num:02d}_df.csv"):
            with open(f"{Test_TPTN_folder}/{file}") as y_Test_file:
                y_Test_list = [item.strip().split("\t") for item in y_Test_file.readlines()]
            col_names = y_Test_list[0]
            values = [row for row in y_Test_list[1:]]
            y_Test_df = pd.DataFrame(values, columns=['SeqName'] + col_names).set_index('SeqName')
            y_Test_df = y_Test_df.astype(int).to_numpy().ravel()
    return X_Train_df, y_Train_df, X_Test_df, y_Test_df

def get_GridSearch_Hyperparameters():
    with open(f"{os.path.dirname(os.path.realpath(__file__))}/GridSearch_Hyperparameters.json") as Hyperpars_file:
        Hyperpars = json.load(Hyperpars_file)

    return Hyperpars

def GridSearch(X_Train, y_Train, models, key):
    
    scoring = {
    'AUC': make_scorer(roc_auc_score),
    'F1': make_scorer(f1_score),
    'MCC': make_scorer(matthews_corrcoef)
    }
    cv = StratifiedKFold(n_splits=5)
    model = models[key]
    model_name = key
    model_Hyperpars = Hyperpars[model_name]
    model_Hyperpars = {key: value for key, value in model_Hyperpars.items()}
    pipe = Pipeline([
    ('scaler', MinMaxScaler()),  # Scaler for LogisticRegression
    ('model', model)  # Placeholder for the model
    ])
    grid_search = GridSearchCV(pipe, model_Hyperpars, cv=cv, scoring=scoring, refit="F1", n_jobs=-1)
    grid_search.fit(X_Train, y_Train)
    return grid_search

def get_GridSearch_Result(grid_search, copy_num, results_dict, X_Test, y_Test):
    
    best_model = grid_search.best_estimator_
    model_name = best_model.named_steps['model']
    best_params = grid_search.best_params_
    best_scores = grid_search.cv_results_


    y_test_pred = best_model.predict(X_Test) ## X_Test has no need to scale manually. Pipeline() will do it for you.
    f1_test = f1_score(y_Test, y_test_pred)
    y_test_pred_proba = best_model.predict_proba(X_Test)[:, 1] if hasattr(best_model, "predict_proba") else None
    aucroc_test = roc_auc_score(y_Test, y_test_pred_proba) if y_test_pred_proba is not None else 'N/A'
    mcc_test = matthews_corrcoef(y_Test, y_test_pred)

    output = f"Model: {results_dict['model_name']}\n"
    output += f"Copy Number: {copy_num:02d}\n"
    output += "Best Parameters: {}\n".format(best_params)
    output += f"AUC: {aucroc_test}\n"
    output += f"F1: {f1_test}\n"
    output += f"MCC: {mcc_test}\n"
    output += "---------------------------------\n"


    if hasattr(model_name, 'feature_importances_'):
        importances = model_name.feature_importances_
    elif hasattr(model_name, 'coef_'):
        importances = np.abs(model_name.coef_[0])
    else:
        raise ValueError("no feature_importances founded!")

    if importances is not None:
        feature_names = X_Test.columns if hasattr(X_Test, 'columns') else np.arange(X_Test.shape[1])
        importance_df = pd.DataFrame({'Feature': feature_names, 'Importance': importances})
        importance_df = importance_df.sort_values(by='Importance', ascending=False)
        importance_df.insert(0, "Rank", range(1, len(importance_df["Feature"])+1))
        importance_df["Percentile"] = importance_df["Importance"].rank(pct=True, ascending=False) * 100

        feature_importance_file = f"{Train_Result_folder}/{results_dict['model_name']}_feature_importances_copy_{copy_num:02d}.txt"
        importance_df.to_csv(feature_importance_file, sep='\t', index=False)

    results_dict['output'] += output
    return best_model, aucroc_test, mcc_test

def main(new_folder, copies, alg="RandomForest", random_state=42):
    global Train_TPTN_folder
    global Test_TPTN_folder
    global Train_Result_folder
    Train_TPTN_folder = f"{new_folder}/Kmer/Train"
    Test_TPTN_folder = f"{new_folder}/Kmer/Test"
    Train_Result_folder = f"{new_folder}/Kmer/ML_Train_Result"

    global Hyperpars
    Hyperpars = get_GridSearch_Hyperparameters()

    Models_Dict = {}
    for copy_num in range(1, copies+1):
        print(f"Fitting Set{copy_num:02d}")
        X_Train_df, y_Train_df, X_Test_df, y_Test_df = read_Train_test_files(copy_num)


        models = {
        'LogisticRegression': LogisticRegression(max_iter=5000),
        'RandomForest': RandomForestClassifier(random_state=random_state),
        'GradientBoosting': GradientBoostingClassifier(random_state=random_state),
        'SVM': SVC(random_state=random_state)
        }
        results_dict = {'output': ''}
        for key in models:
            if alg == key:
                print(f"Fitting {key}... Random State: {random_state}")
                grid_search = GridSearch(X_Train_df, y_Train_df, models, key)
                results_dict['model_name'] = key
                best_model, aucroc_test, mcc_test = get_GridSearch_Result(grid_search, copy_num, results_dict, X_Test_df, y_Test_df)
                Models_Dict.update({str(copy_num): [best_model, aucroc_test, mcc_test]})
                print(f"{key} Done!")


        result_filename = f"{Train_Result_folder}/copy_{copy_num:02d}_all_models_score.txt"
        with open(result_filename, 'w') as f:
            f.write(results_dict['output'])
    
    aucroc_values = [value[1] for value in Models_Dict.values()]
    mcc_values = [value[2] for value in Models_Dict.values()]
    aucroc_rank = np.argsort(-np.array(aucroc_values))
    mcc_rank = np.argsort(-np.array(mcc_values))
    total_rank = aucroc_rank + mcc_rank
    BestCopy = str(np.argmin(total_rank)+1)
    BestCopy_Model = Models_Dict[BestCopy][0]
    dump(BestCopy_Model, f'{Train_Result_folder}/BestModel_Copy{BestCopy}.joblib')

    return BestCopy

