import os
import ahocorasick
import pandas as pd
import itertools

def RC(kmer: str) -> str:
    """Return reverse complement of DNA sequence."""
    complement = str.maketrans("ATCG", "TAGC")
    return kmer.translate(complement)[::-1]

def build_automaton(kmers):
    auto = ahocorasick.Automaton()
    for kmer in kmers:
        if "n" not in kmer:
            auto.add_word(kmer, kmer)
        else:
            n_count = kmer.count("n")
            for repl in itertools.product("ATCG", repeat=n_count):
                expanded = list(kmer)
                repl_iter = iter(repl)
                expanded = "".join(next(repl_iter) if c == "n" else c for c in expanded)
                auto.add_word(expanded, kmer)
    auto.make_automaton()
    return auto

def KmerCount_inSeq(seq: str, auto, kmer_list):
    """Count Kmers in seq and its reverse complement."""
    rc_seq = RC(seq)
    count_seq = {k: 0 for k in kmer_list}
    count_rc = {k: 0 for k in kmer_list}

    for _, kmer in auto.iter(seq):
        count_seq[kmer] += 1
    for _, kmer in auto.iter(rc_seq):
        count_rc[kmer] += 1

    return list(count_seq.values()) + list(count_rc.values())

def get_Kmer_list(folder, copy_num, pThreshold):
    """Extract Kmers from Enriched_Kmer file."""
    for file in sorted(os.listdir(folder)):
        if file.endswith(f"{copy_num:02d}.ac_WRSresults.txt"):
            with open(os.path.join(folder, file)) as f:
                lines = f.readlines()[1:]
                return [line.split("\t")[0] for line in lines if (float(line.split("\t")[-3]) < pThreshold and 
                                                                  (line.split("\t")[-1].strip()=="TP" or line.split("\t")[-1].strip()=="TN"))]
                # return [line.split("\t")[0] for line in lines]
    return []

def parse_fasta_and_count(folder, suffix, auto, kmer_list, label, copy_num):
    """Parse fasta file and count Kmers."""
    result = {}
    for file in sorted(os.listdir(folder)):
        if file.endswith(f"{copy_num:02d}{suffix}"):
            with open(os.path.join(folder, file)) as f:
                lines = f.readlines()
                for i, line in enumerate(lines):
                    if line.startswith(">"):
                        seq_name = line[1:].split()[0]
                        seq = lines[i + 1].strip()
                        kmer_counts = KmerCount_inSeq(seq, auto, kmer_list)
                        result[seq_name] = [seq, kmer_counts, label]
    return result

def get_sequences(folder, copy_num, suffix, label, pThreshold, kmer_source_folder=None):
    """Generalized sequence loader."""
    if kmer_source_folder is None:
        kmer_source_folder = folder
    kmers = get_Kmer_list(kmer_source_folder, copy_num, pThreshold)
    auto = build_automaton(kmers)
    return parse_fasta_and_count(folder, suffix, auto, kmers, label, copy_num)

def Make_DFs(train_folder, test_folder, copy_num, pThreshold):
    kmers = get_Kmer_list(train_folder, copy_num, pThreshold)
    kmer_cols = ["nt_" + k for k in kmers] + ["t_" + k for k in kmers]

    Train_TP = get_sequences(train_folder, copy_num, ".Btp.fa", 1, pThreshold)
    Train_TN = get_sequences(train_folder, copy_num, ".Btn.fa", 0, pThreshold)
    Test_TP = get_sequences(test_folder, copy_num, ".tp.fa", 1, pThreshold, kmer_source_folder=train_folder)
    Test_TN = get_sequences(test_folder, copy_num, ".tn.fa", 0, pThreshold, kmer_source_folder=train_folder)

    Train_X = pd.DataFrame.from_dict({k:v[1] for k,v in {**Train_TP, **Train_TN}.items()},
                                     orient="index", columns=kmer_cols)
    Train_y = pd.DataFrame.from_dict({k:v[2] for k,v in {**Train_TP, **Train_TN}.items()},
                                     orient="index", columns=["Class"])

    Test_X = pd.DataFrame.from_dict({k:v[1] for k,v in {**Test_TP, **Test_TN}.items()},
                                    orient="index", columns=kmer_cols)
    Test_y = pd.DataFrame.from_dict({k:v[2] for k,v in {**Test_TP, **Test_TN}.items()},
                                    orient="index", columns=["Class"])

    return Train_X, Train_y, Test_X, Test_y

def main(new_folder, copies, pThreshold):
    train_folder = f"{new_folder}/Kmer/Train"
    test_folder = f"{new_folder}/Kmer/Test"
    for copy_num in range(1, copies + 1):
        print(f"Making fold{copy_num:02d} DF...")
        Train_X, Train_y, Test_X, Test_y = Make_DFs(train_folder, test_folder, copy_num, pThreshold)
        Train_X.to_csv(f"{train_folder}/Train_X_{copy_num:02d}_df.csv", sep="\t")
        Train_y.to_csv(f"{train_folder}/Train_y_{copy_num:02d}_df.csv", sep="\t")
        Test_X.to_csv(f"{test_folder}/Test_X_{copy_num:02d}_df.csv", sep="\t")
        Test_y.to_csv(f"{test_folder}/Test_y_{copy_num:02d}_df.csv", sep="\t")
