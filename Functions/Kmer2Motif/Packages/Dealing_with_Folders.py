import shutil
import os
from datetime import datetime

def main():
    now = datetime.now()
    formatted_date = now.strftime("%Y_%m_%d")

    
    i = 1
    while True:
        pre_destination_dir = f"./Results/{formatted_date}_Kmer2Motif{i:02d}"
        if not os.path.exists(pre_destination_dir):
            destination_dir = pre_destination_dir
            break
        i += 1
    source_dir = './InputData/Kmer2Motif'
    shutil.copytree(source_dir, destination_dir)
    print(f"Destination_dir({destination_dir}) copied!")
    shutil.rmtree(source_dir)
    print(f"Source_dir({source_dir}) removed!")

    To_make_dirs = ["./InputData/Kmer2Motif/KmerList", "./InputData/Kmer2Motif/MotifList", f"{destination_dir}/temp_Rscripts"]
    for To_make_dir in To_make_dirs:
        os.makedirs(To_make_dir)
    print(f"Source_dir({source_dir}) recreated!")
    return destination_dir
