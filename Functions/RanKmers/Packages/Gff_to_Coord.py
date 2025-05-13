##### OUTLINE #####
#============================================================#
    ## 1. import part ##
#============================================================#
    ## 2. def part ##
        ## i)   gff extract coord  ##
        ## ii)  get gff_file path  ##
        ## iii) main               ##
#============================================================#

##### OUTLINE END#####


#============================================================#
#============================================================#
## 1. import part ##
import argparse
import os

#============================================================#
## 2. def part ##
def gff_prom_to_coord_5utr(gff, features, up_stream, down_stream):
    print(f"features: {features}")
    inp = open(gff)
    oup = open(gff+"_prom-5utr.coord","w")
    inl = inp.readline()
    tmp_list =[]
    count = 0
    while inl != "":
        one_row = inl.strip().split("\t")              # a row in gff
        if len(one_row) > 6:
            Chromosome = one_row[0] # one_row[0] == Chromosome number e.g.Chr1, Chr2...
            #if ("Chr" in one_row[0] or "chr" in one_row[0])else "Chr"+one_row[0] 
            Feature = one_row[2]                       # one_row[2] == Feature type (gene, mRNA...)
            if count == 0 and (Feature in features):  # "feature" is a inputed argument
                count += 1
                tmp_list=[]

                if one_row[6] == "+":                # forward orientation
                    Left_Coord = one_row[3]              # left coord
                    Right_Coord = one_row[4]             # right coord
                    prom1= int(Left_Coord)-up_stream    # get 1000 bp upstream of TSS
                    prom2= int(Left_Coord)+down_stream     # get 200bp downstream of TSS

                elif one_row[6] == "-":              # reverse orientation
                    Left_Coord = one_row[4]
                    Right_Coord = one_row[3]
                    prom1 = int(Left_Coord)+up_stream
                    prom2 = int(Left_Coord)-down_stream

                if prom1 > 0:
                    tmp_list=[prom2, prom1]
                else:
                    tmp_list=[prom2, 0]
                Sequce_name = ""                     # sequence name
                sequence_IDs_cell = one_row[-1].split(";")     # get gene name. "cell" means a cross of a column and a row.
                if "Name" in one_row[-1]:            # has Name tag, # refers to last item in list
                    for name in sequence_IDs_cell:
                        if "Name" in name:
                            Sequce_name = name.split("=")[1]
                            break
                        elif one_row[-1] != "":            # no name tag but not empty, use 1st
                            Sequce_name = sequence_IDs_cell[0].split("=")[1]
                        else:
                            print ("No desc:",one_row)
                if Sequce_name != "":
                    L2= tmp_list[0]
                    prom= tmp_list[1]
                    oup.write("%s\t%s,%s\t%s\n" % (Chromosome,prom,L2,Sequce_name))
                count = 0
        
        else:
            pass
        inl = inp.readline()
    return

def get_gff_file_path(new_folder):
    cwd = os.getcwd()
    print(f"Gff_to_Coord working direction: {cwd}")
    Gff_and_Genome_folder = f"{new_folder}/Gff_and_Genome"
    Gff_and_Genome_files = sorted(os.listdir(Gff_and_Genome_folder))
    gff_file = ""
    for file in Gff_and_Genome_files:
        if file.endswith(".gff") or file.endswith(".gff3"):
            gff_file = os.path.join(cwd, Gff_and_Genome_folder, file)
            break
        else:
            continue

    if gff_file == "":
        raise NameError("Can't get gff file name. Make sure it's end with '.gff' or '.gff3'")
    else:
        print(f"gff_filename: {gff_file}")
    return gff_file

def main(new_folder, features, up_stream, down_stream):
    print("Now executing Gff_to_Coord.py")
    gff_prom_to_coord_5utr(gff=get_gff_file_path(new_folder=new_folder), features=features, up_stream=up_stream, down_stream=down_stream)