import shutil
import os
from datetime import datetime

def FileExtention_changer(ParentFolder, FE, ori_path):
    basename = os.path.basename(ori_path)
    return os.path.join(ParentFolder, f"{basename}{FE}")

def ClearFolder(source_dir, ChildFolder):
    for filename in os.listdir(os.path.join(source_dir, ChildFolder)):
        file_path = os.path.join(source_dir, ChildFolder, filename)
        if os.path.isfile(file_path) or os.path.islink(file_path):
            os.unlink(file_path)
        elif os.path.isdir(file_path):
            shutil.rmtree(file_path)

def main(tp, tn, Motif):
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
    ClearFolder(source_dir, "GeneXKmerTable")
    ClearFolder(source_dir, "MotifList")

    shutil.copyfile(tp, os.path.join(source_dir, FileExtention_changer("GeneXKmerTable", "", tp)))
    shutil.copyfile(tn, os.path.join(source_dir, FileExtention_changer("GeneXKmerTable", "", tn)))
    if Motif != "":
        shutil.copyfile(Motif, os.path.join(source_dir, FileExtention_changer("MotifList", ".motif", Motif)))

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
