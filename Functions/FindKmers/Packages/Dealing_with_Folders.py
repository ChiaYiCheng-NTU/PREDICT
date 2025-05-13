import shutil
import os
from datetime import datetime

def main():
    now = datetime.now()
    formatted_date = now.strftime("%Y_%m_%d")

    i = 1
    while True:
        if i < 10:
            pre_destination_dir = f"./Results/{formatted_date}_FindKmers0{i}"
            if not os.path.exists(pre_destination_dir):
                destination_dir = pre_destination_dir
                break
        elif i >= 10:
            pre_destination_dir = f"./Results/{formatted_date}_FindKmers{i}"
            if not os.path.exists(pre_destination_dir):
                destination_dir = pre_destination_dir
                break
        i += 1

    source_dir = './InputData/FindKmers'
    shutil.copytree(source_dir, destination_dir)
    print(f"Destination_dir({destination_dir}) copied!")
    shutil.rmtree(source_dir)
    print(f"Source_dir({source_dir}) removed!")

    To_make_dirs = ["./InputData/FindKmers/Gff_and_Genome", "./InputData/FindKmers/TP_and_TN", f"{pre_destination_dir}/Kmer/ML_Train_Result", f"{pre_destination_dir}/UsefulOutputs"]
    for To_make_dir in To_make_dirs:
        os.makedirs(To_make_dir)
    print(f"Source_dir({source_dir}) recreated!")

    return destination_dir