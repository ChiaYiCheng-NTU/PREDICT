import os
import ahocorasick

def RC(kmer):
    complement = str.maketrans("ACGT", "TGCA")  # DNA complement mapping
    return kmer.translate(complement)[::-1]

def build_automaton(kmers):
    """
    Build an Aho-Corasick automaton for Kmers and their reverse complements.
    Handles potential collision Kmer.
    """
    auto = ahocorasick.Automaton()
    for kmer in kmers:
        auto.add_word(kmer, kmer)
    auto.make_automaton()
    return auto

def KmerCount_inSeq(temp_tp_dic, seq_name, Auto, Kmer_list):
    seq = temp_tp_dic[seq_name][0]
    rc_seq = RC(seq)
    KmerCount_inSeq_Dict = {kmer:0 for kmer in Kmer_list}
    KmerCount_inRCSeq_Dict = {kmer:0 for kmer in Kmer_list}
    for _, kmer in Auto.iter(seq):
        KmerCount_inSeq_Dict[kmer] += 1
    for _, kmer in Auto.iter(rc_seq):
        KmerCount_inRCSeq_Dict[kmer] += 1

    KmerCount_inSeq_list = list(KmerCount_inSeq_Dict.values()) + (list(KmerCount_inRCSeq_Dict.values()))

    return KmerCount_inSeq_list

def get_tp(Train_TPTN_folder, copy_num):
    tp_dic = {}
    Kmer_list = get_Kmer_list(Train_TPTN_folder, copy_num)
    TP_Auto = build_automaton(Kmer_list)
    ### get tp sequence_name & sequnce
    for file in sorted(os.listdir(Train_TPTN_folder)):
        if file.endswith(f"{copy_num:02d}.Btp.fa"):
            with open(f"{Train_TPTN_folder}/{file}") as tpfa:
                lines = tpfa.readlines()
                for index, line in enumerate(lines):
                    if ">" in line:
                        temp_tp_dic = {}
                        seq_name = line.split(">")[1].split(" ")[0]   ## pure sequence name
                        temp_tp_dic[seq_name] = [lines[index+1]]   ## temp_tp_dic's key==SeqName, value==[seq,].
                        KmerCount_list = KmerCount_inSeq(temp_tp_dic, seq_name, TP_Auto, Kmer_list)
                        temp_tp_dic[seq_name].append(KmerCount_list)    ## temp_tp_dic's key==SeqName, value==[seq, [KmerCount_inBothOrientationSeq]].
                        temp_tp_dic[seq_name].append(1)    ## temp_TestTn_dic's key==SeqName, value==[seq, [KmerCount_inBothOrientationSeq], Class==1].
                        tp_dic.update(temp_tp_dic)
    return tp_dic

def get_tn(Train_TPTN_folder, copy_num):
    tn_dic = {}
    Kmer_list = get_Kmer_list(Train_TPTN_folder, copy_num)
    TN_Auto = build_automaton(Kmer_list)
    ### get tn sequence_name & sequnce
    for file in sorted(os.listdir(Train_TPTN_folder)):
        if file.endswith(f"{copy_num:02d}.Btn.fa"):
            with open(f"{Train_TPTN_folder}/{file}") as tnfa:
                lines = tnfa.readlines()
                for index, line in enumerate(lines):
                    if ">" in line:
                        temp_tn_dic = {}
                        seq_name = line.split(">")[1].split(" ")[0]   ## pure sequence name
                        temp_tn_dic[seq_name] = [lines[index+1]]   ## temp_tn_dic's key==SeqName, value==[seq,].
                        KmerCount_list = KmerCount_inSeq(temp_tn_dic, seq_name, TN_Auto, Kmer_list)
                        temp_tn_dic[seq_name].append(KmerCount_list)    ## temp_tn_dic's key==SeqName, value==[seq, [KmerCount_inBothOrientationSeq]].
                        temp_tn_dic[seq_name].append(0)    ## temp_TestTn_dic's key==SeqName, value==[seq, [KmerCount_inBothOrientationSeq], Class==0].
                        tn_dic.update(temp_tn_dic)
    return tn_dic

def get_Kmer_list(Train_TPTN_folder, copy_num):
    ## extract Kmers from Enriched_Kmer file
    Kmer_list = []
    for file in sorted(os.listdir(Train_TPTN_folder)):
        if file.endswith(f"{copy_num:02d}.pcre_WRSresults.txt"):
            Kmer_file = os.path.join(Train_TPTN_folder, file)
            break
    with open(Kmer_file) as Kmer_file:
        lines = Kmer_file.readlines()
        lines = lines[1:]
        for line in lines:
            Kmer = line.split("\t")[0]
            Kmer_list.append(Kmer)
    return Kmer_list

def get_test_tp(Test_TPTN_folder, Train_TPTN_folder, copy_num):
    TestTp_dic = {}
    Kmer_list = get_Kmer_list(Train_TPTN_folder, copy_num)
    Test_TP_Auto = build_automaton(Kmer_list)
    for file in sorted(os.listdir(Test_TPTN_folder)):
        if file.endswith(f"{copy_num:02d}.tp.fa"):
            with open(f"{Test_TPTN_folder}/{file}") as TestTpfa:
                lines = TestTpfa.readlines()
                for index, line in enumerate(lines):
                    if ">" in line:
                        temp_TestTp_dic = {}
                        seq_name = line.split(">")[1].split(" ")[0]   ## pure sequence name
                        temp_TestTp_dic[seq_name] = [lines[index+1]]   ## temp_TestTp_dic's key==SeqName, value==[seq,].
                        KmerCount_list = KmerCount_inSeq(temp_TestTp_dic, seq_name, Test_TP_Auto, Kmer_list)
                        temp_TestTp_dic[seq_name].append(KmerCount_list)    ## temp_TestTp_dic's key==SeqName, value==[seq, [KmerCount_inBothOrientationSeq]].
                        temp_TestTp_dic[seq_name].append(1)    ## temp_TestTn_dic's key==SeqName, value==[seq, [KmerCount_inBothOrientationSeq], Class==1].
                        TestTp_dic.update(temp_TestTp_dic)
    return TestTp_dic

def get_test_tn(Test_TPTN_folder, Train_TPTN_folder, copy_num):
    TestTn_dic = {}
    Kmer_list = get_Kmer_list(Train_TPTN_folder, copy_num)
    Test_TN_Auto = build_automaton(Kmer_list)
    for file in sorted(os.listdir(Test_TPTN_folder)):
        if file.endswith(f"{copy_num:02d}.tn.fa"):
            with open(f"{Test_TPTN_folder}/{file}") as TestTnfa:
                lines = TestTnfa.readlines()
                for index, line in enumerate(lines):
                    if ">" in line:
                        temp_TestTn_dic = {}
                        seq_name = line.split(">")[1].split(" ")[0]   ## pure sequence name
                        temp_TestTn_dic[seq_name] = [lines[index+1]]   ## temp_TestTn_dic's key==SeqName, value==[seq,].
                        KmerCount_list = KmerCount_inSeq(temp_TestTn_dic, seq_name, Test_TN_Auto, Kmer_list)
                        temp_TestTn_dic[seq_name].append(KmerCount_list)    ## temp_TestTn_dic's key==SeqName, value==[seq, [KmerCount_inBothOrientationSeq]].
                        temp_TestTn_dic[seq_name].append(0)    ## temp_TestTn_dic's key==SeqName, value==[seq, [KmerCount_inBothOrientationSeq], Class==0].
                        TestTn_dic.update(temp_TestTn_dic)
    return TestTn_dic

def Make_DFs(Train_TPTN_folder, Test_TPTN_folder, copy_num):
    import pandas as pd
    Train_TP_dic = get_tp(Train_TPTN_folder, copy_num)
    Train_TN_dic = get_tn(Train_TPTN_folder, copy_num)
    Test_TP_dic = get_test_tp(Test_TPTN_folder, Train_TPTN_folder, copy_num)
    Test_TN_dic = get_test_tn(Test_TPTN_folder, Train_TPTN_folder, copy_num)
    Kmer_list = ["nt_" + Kmer for Kmer in get_Kmer_list(Train_TPTN_folder, copy_num)] + ["t_" + Kmer for Kmer in get_Kmer_list(Train_TPTN_folder, copy_num)]

    Train_X_df = pd.DataFrame(columns=Kmer_list)
    Train_y_df = pd.DataFrame(columns=["Class"])
    
    for key in Train_TP_dic.keys():
        Train_X_df.loc[key] = Train_TP_dic[key][1]  ## pure "KmerCount_inBothOrientationSeq" counts list
        Train_y_df.loc[key] = Train_TP_dic[key][2]  ## class "1" of SourceData
    for key in Train_TN_dic.keys():
        Train_X_df.loc[key] = Train_TN_dic[key][1]  ## pure "KmerCount_inBothOrientationSeq" counts list
        Train_y_df.loc[key] = Train_TN_dic[key][2]  ## class "0" of SourceData

    Test_X_df = pd.DataFrame(columns=Kmer_list)
    Test_y_df = pd.DataFrame(columns=["Class"])
    for key in Test_TP_dic.keys():
        Test_X_df.loc[key] = Test_TP_dic[key][1]  ## pure "KmerCount_inBothOrientationSeq" counts list
        Test_y_df.loc[key] = Test_TP_dic[key][2]  ## class "1" of SourceData
    for key in Test_TN_dic.keys():
        Test_X_df.loc[key] = Test_TN_dic[key][1]  ## pure "KmerCount_inBothOrientationSeq" counts list
        Test_y_df.loc[key] = Test_TN_dic[key][2]  ## class "0" of SourceData
        

    return Train_X_df, Train_y_df, Test_X_df, Test_y_df

def main(new_folder, copies):
    Train_TPTN_folder = f"{new_folder}/Kmer/Train"
    Test_TPTN_folder = f"{new_folder}/Kmer/Test"
    for copy_num in range(1, copies+1):
        print(f"Making fold{copy_num:02d} DF...")
        Train_X_df, Train_y_df, Test_X_df, Test_y_df = Make_DFs(Train_TPTN_folder, Test_TPTN_folder, copy_num=copy_num)
        
        Train_X_df.to_csv(f"{Train_TPTN_folder}/Train_X_{copy_num:02d}_df.csv", sep="\t")
        Train_y_df.to_csv(f"{Train_TPTN_folder}/Train_y_{copy_num:02d}_df.csv", sep="\t")
        Test_X_df.to_csv(f"{Test_TPTN_folder}/Test_X_{copy_num:02d}_df.csv", sep="\t")
        Test_y_df.to_csv(f"{Test_TPTN_folder}/Test_y_{copy_num:02d}_df.csv", sep="\t")