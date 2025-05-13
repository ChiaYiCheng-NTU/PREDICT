import os
import pandas

def get_GeneXKmerDF(new_folder):
    GeneXKmerDF_Folder = os.path.join(new_folder, "GeneXKmerTable")
    TPGeneXKmerDF_File = sorted(os.listdir(GeneXKmerDF_Folder))[0]
    TNGeneXKmerDF_File = sorted(os.listdir(GeneXKmerDF_Folder))[1]
    TPGeneXKmerDF_Path = os.path.join(GeneXKmerDF_Folder, TPGeneXKmerDF_File)
    TNGeneXKmerDF_Path = os.path.join(GeneXKmerDF_Folder, TNGeneXKmerDF_File)
    with open(TPGeneXKmerDF_Path) as f:
        TPlines = [line.strip().split("\t") for line in f.readlines()]
        TPindex = [line[0] for line in TPlines[1:]]
    with open(TNGeneXKmerDF_Path) as f:
        TNlines = [line.strip().split("\t") for line in f.readlines()]
        TNindex = [line[0] for line in TNlines[1:]]
    TPcolumns = TPlines[0]
    TNcolumns = TNlines[0]
    TPdata = [line[1:] for line in TPlines[1:]]
    TNdata = [line[1:] for line in TNlines[1:]]
    TPGeneXKmerDF = pandas.DataFrame(data=TPdata, columns=TPcolumns, index=TPindex).astype(int)
    TNGeneXKmerDF = pandas.DataFrame(data=TNdata, columns=TNcolumns, index=TNindex).astype(int)
    return TPGeneXKmerDF, TNGeneXKmerDF


def main(new_folder):
    TPGeneXKmerDF, TNGeneXKmerDF = get_GeneXKmerDF(new_folder)
    return TPGeneXKmerDF, TNGeneXKmerDF

