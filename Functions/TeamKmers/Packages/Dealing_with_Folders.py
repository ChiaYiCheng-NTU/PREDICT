import shutil
import os
from datetime import datetime

def main():
    now = datetime.now()
    formatted_date = now.strftime("%Y_%m_%d")

    
    i = 1
    while True:
        pre_destination_dir = f"{os.path.dirname(os.path.realpath(__file__))}/../../../Results/{formatted_date}_TeamKmers{i:02d}"
        if not os.path.exists(pre_destination_dir):
            destination_dir = pre_destination_dir
            break
        i += 1
    source_dir = f'{os.path.dirname(os.path.realpath(__file__))}/../../../InputData/TeamKmers'
    shutil.copytree(source_dir, destination_dir)
    print(f"Destination_dir({destination_dir}) copied!")
    shutil.rmtree(source_dir)
    print(f"Source_dir({source_dir}) removed!")

    To_make_dirs = [f"{source_dir}/GeneXKmerTable",
                    f"{source_dir}/MotifList",
                    f"{pre_destination_dir}/Outputs/Kmer2Motif_inputs/K1",
                    f"{pre_destination_dir}/Outputs/Kmer2Motif_inputs/K2"]
    for To_make_dir in To_make_dirs:
        os.makedirs(To_make_dir)
    print(f"Source_dir({source_dir}) recreated!")
    return destination_dir
