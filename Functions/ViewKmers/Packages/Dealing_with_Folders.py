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

def main(gff, genome, gene, Kmer, Motif):
    now = datetime.now()
    formatted_date = now.strftime("%Y_%m_%d")

    i = 1
    while True:
        pre_destination_dir = f"{os.path.dirname(os.path.realpath(__file__))}/../../../Results/{formatted_date}_ViewKmers{i:02d}"
        if not os.path.exists(pre_destination_dir):
            destination_dir = pre_destination_dir
            break
        i += 1

    source_dir = f'{os.path.dirname(os.path.realpath(__file__))}/../../../InputData/ViewKmers'
    ClearFolder(source_dir, "GeneList")
    ClearFolder(source_dir, "Gff_and_Genome")
    ClearFolder(source_dir, "KmerList")
    ClearFolder(source_dir, "MotifList")
    
    shutil.copyfile(gff, os.path.join(source_dir, FileExtention_changer("Gff_and_Genome", ".gff", gff)))
    shutil.copyfile(genome, os.path.join(source_dir, FileExtention_changer("Gff_and_Genome", ".fasta", genome)))
    shutil.copyfile(gene, os.path.join(source_dir, FileExtention_changer("GeneList", ".gene", gene)))
    shutil.copyfile(Kmer, os.path.join(source_dir, FileExtention_changer("KmerList", ".kmer", Kmer)))
    shutil.copyfile(Motif, os.path.join(source_dir, FileExtention_changer("MotifList", ".motif", Motif)))

    shutil.copytree(source_dir, destination_dir)
    print(f"Destination_dir({destination_dir}) copied!")
    shutil.rmtree(source_dir)
    print(f"Source_dir({source_dir}) removed!")

    To_make_dirs = [f"{source_dir}/Gff_and_Genome",
                    f"{source_dir}/GeneList",
                    f"{source_dir}/KmerList",
                    f"{source_dir}/MotifList"]
    for To_make_dir in To_make_dirs:
        os.makedirs(To_make_dir)
    print(f"Source_dir({source_dir}) recreated!")

    return destination_dir