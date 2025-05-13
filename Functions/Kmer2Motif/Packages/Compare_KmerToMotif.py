import subprocess
import os
import pandas as pd
from multiprocessing import Pool


def Make_CompareMotifs_command(new_folder, CreateMotif_command, Kmer):
    motifs_list = [motif_ID.split(" <- ")[0] for motif_ID in CreateMotif_command.split("\n")]
    vector = "c(" + ", ".join(motifs_list) + ")"
    command_str = f"compare_motifs_result <- compare_motifs({vector}, compare.to = 1, max.p = 0.01, max.e = Inf, nthreads = 0, tryRC = FALSE)\n"
    command_str += "Output_Str <- ''\n"
    command_str += "for (rowNum in 1:nrow(compare_motifs_result)) {\n    OutMotif <- compare_motifs_result[rowNum, 3]\n    OutScore <- compare_motifs_result[rowNum, 5]\n    Output_Str <- paste(Output_Str, paste(OutMotif, OutScore, sep='>'), sep=',')}\n"
    command_str += "print(Output_Str)"
    return command_str

def Make_CreateMotif_command(KmersOrMotifs, MotifInfoList):
    commands_list = []
    for index, element in enumerate(KmersOrMotifs):
        if index == 0:
            command = f"Motif_{index+1:03d} <- create_motif(input='{element}', name='{element}')"
        else:
            command = f"Motif_{index+1:03d} <- create_motif(input='{element}', name='{element}<{MotifInfoList[index-1]}')"
        commands_list.append(command)
    command_str = "\n".join(commands_list)
    return command_str

def Make_tempRscript(new_folder, CreateMotif_command, CompareMotifs_command, Ori_Kmer):
    tempRscript_path = os.path.join(new_folder, f"temp_Rscripts/{Ori_Kmer}_temp_Rscript.R")
    with open(tempRscript_path, "w") as f:
        f.write("library(universalmotif)\n")
        f.write(f"{CreateMotif_command}\n")
        f.write(CompareMotifs_command)
    return tempRscript_path

class MotifList_preprocess:
    def __init__(self, MotifList, Kmer):
        if self.MotifListIsFasta(MotifList) == True:
            ## MotifList Format: [[Motif001, MotifFamily001], [Motif002, MotifFamily002], ...]
            self.MotifList = [Combination[0] for Combination in MotifList]
            self.MotifInfo = [Combination[1] for Combination in MotifList]
        else:
            self.MotifList = MotifList
            self.MotifInfo = [None] * len(MotifList)
        self.Kmer = Kmer
        self.Klength = len(self.Kmer)
        self.Remove_SideNs()
        self.Remove_TooManyNsMotif()
        self.Remove_SideNs()
        self.Insert_Kmer()

    def MotifListIsFasta(self, MotifList):
        if type(MotifList[0]) == list:
            return True
        else:
            return False

    def Remove_SideNs(self):
        self.MotifList = [motif.strip("N") for motif in self.MotifList]

    def Remove_TooManyNsMotif(self):
        # self.lessN_MotifList = [motif for motif in self.MotifList if motif.count("N"*(self.Klength-1)) == 0]
        # self.MuchN_MotifList = [motif for motif in self.MotifList if motif.count("N"*(self.Klength-1)) != 0]
        # self.Mod_MuchN_MotifList = self.Keep_LongestFragment(self.MuchN_MotifList)
        # self.MotifList = self.lessN_MotifList + self.Mod_MuchN_MotifList

        lessN_MotifList = [(idx, motif) for idx, motif in enumerate(self.MotifList) if motif.count("N" * (self.Klength - 1)) == 0]
        muchN_MotifList = [(idx, motif) for idx, motif in enumerate(self.MotifList) if motif.count("N" * (self.Klength - 1)) != 0]
        mod_MuchN_MotifList = self.Keep_LongestFragment([motif for _, motif in muchN_MotifList])
        merged_MotifList = lessN_MotifList + [(idx, motif) for (idx, _), motif in zip(muchN_MotifList, mod_MuchN_MotifList)]
        merged_MotifList.sort(key=lambda x: x[0])
        self.MotifList = [motif for _, motif in merged_MotifList]
        
    def Keep_LongestFragment(self, MuchN_MotifList):
        Mod_MuchN_MotifList = [max(motif.split("NNN"), key=len) for motif in MuchN_MotifList]
        return Mod_MuchN_MotifList

    def Insert_Kmer(self):
        self.MotifList.insert(0, self.Kmer)

def GenerateOutputFormat(CompareMotifs_output, Ori_Kmer, KmerRank, ScoreCutOff, KeepTopMotifs):
    CompareMotifs_Result_List = [line.strip("\n").strip('"').split(">") 
                                for line in CompareMotifs_output.stdout.split(",")[1:]]
                                ## Format: [[Motif001, Score001], [Motif002, Score002], ...]

    ToMakeDF_List = []
    for Motif in CompareMotifs_Result_List:
        MotifSeq = Motif[0].split("<")[0]
        MotifInfo = Motif[0].split("<")[1]
        MotifScore = Motif[1]
        Motif = [Ori_Kmer, KmerRank+1, MotifSeq, MotifInfo, MotifScore]
        ToMakeDF_List.append(Motif)

    df = pd.DataFrame(ToMakeDF_List, columns=["KmerSeq", "KmerRank", "MotifSeq", "MotifInfo", "Score"])
    df['Score'] = df['Score'].astype(float)
    df = df[df['Score'] > ScoreCutOff]
    num_rows_to_keep = int(len(df) * KeepTopMotifs)
    df = df.head(num_rows_to_keep)
    return df

def Compare_Process(args):
    Ori_Kmer, KmerRank, MotifList, new_folder, ScoreCutOff, KeepTopMotifs = args
    
    Kmer = Ori_Kmer.strip("nt_").strip("t_")
    Motif_Infomations = MotifList_preprocess(MotifList, Kmer)
    temp_MotifList = Motif_Infomations.MotifList
    temp_MotifInfoList = Motif_Infomations.MotifInfo

    CreateMotif_command = Make_CreateMotif_command(temp_MotifList, temp_MotifInfoList)
    CompareMotifs_command = Make_CompareMotifs_command(new_folder, CreateMotif_command, Kmer)
    tempRscript_path = Make_tempRscript(new_folder, CreateMotif_command, CompareMotifs_command, Ori_Kmer)
    R_command = f"conda run -n PREDICTm Rscript {tempRscript_path}"
    
    CompareMotifs_output = subprocess.run(R_command, shell=True, capture_output=True, text=True)
    
    UniversalMotif_Err = CompareMotifs_output.stderr
    if UniversalMotif_Err != "":
        raise ValueError(f"***If you got this Error, it's very likely that you inputed some Motifs with unrecognized alphabets***\n{UniversalMotif_Err}")

    temp_OutPutDF = GenerateOutputFormat(
        CompareMotifs_output, Ori_Kmer, KmerRank, ScoreCutOff, KeepTopMotifs
    )
    
    return temp_OutPutDF



def main(new_folder, KmerList, MotifList, ScoreCutOff, KeepTopMotifs):
    
    tasks = [
        (Ori_Kmer, KmerRank, MotifList, new_folder, ScoreCutOff, KeepTopMotifs)
        for KmerRank, Ori_Kmer in enumerate(KmerList)
    ]
    
    with Pool() as pool:
        results = pool.map(Compare_Process, tasks)
    OutPutDF = pd.concat(results, ignore_index=True)
    # print(OutPutDF)
    OutPutDF.to_csv(f"{new_folder}/Kmer2Motif.tsv", sep="\t")
    return OutPutDF









