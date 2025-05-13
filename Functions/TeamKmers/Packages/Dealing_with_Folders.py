import shutil
import os
from datetime import datetime

def main():
    now = datetime.now()
    formatted_date = now.strftime("%Y_%m_%d")

    
    i = 1
    while True:
        pre_destination_dir = f"./Results/{formatted_date}_TeamKmers{i:02d}"
        if not os.path.exists(pre_destination_dir):
            destination_dir = pre_destination_dir
            break
        i += 1
    source_dir = './InputData/TeamKmers'
    shutil.copytree(source_dir, destination_dir)
    print(f"Destination_dir({destination_dir}) copied!")
    shutil.rmtree(source_dir)
    print(f"Source_dir({source_dir}) removed!")

    To_make_dirs = ["./InputData/TeamKmers/GeneXKmerTable",
                    "./InputData/TeamKmers/MotifList",
                    f"{pre_destination_dir}/Outputs/Kmer2Motif_inputs/K1",
                    f"{pre_destination_dir}/Outputs/Kmer2Motif_inputs/K2"]
    for To_make_dir in To_make_dirs:
        os.makedirs(To_make_dir)
    print(f"Source_dir({source_dir}) recreated!")
    return destination_dir
