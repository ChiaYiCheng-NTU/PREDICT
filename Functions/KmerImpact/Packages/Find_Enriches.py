import os
import ahocorasick
import multiprocessing
import itertools
from itertools import product
import pandas as pd
import numpy as np
from scipy.stats import norm
from scipy.stats import ranksums
from statsmodels.stats.multitest import multipletests

# -----------------------------
# basic tools
# -----------------------------
def RC(seq):
    return seq.translate(str.maketrans("ATCG", "TAGC"))[::-1]

def atcg_sort_key(kmer):
    order = {"A": 0, "T": 1, "C": 2, "G": 3}
    return (len(kmer), [order.get(b, 99) for b in kmer])

# -----------------------------
# Kmer & Automaton
# -----------------------------
def Generate_Gapped_Kmers(k1=3, k2=3, gap=2):
    heads = Generate_Kmers(k1)
    tails = Generate_Kmers(k2)
    return [h + "n"*gap + t for h in heads for t in tails]

def Generate_Kmers(k):
    return [''.join(p) for p in product("ATCG", repeat=k)]

def Build_Automaton(kmers):
    auto = ahocorasick.Automaton()
    for idx, kmer in enumerate(kmers):
        if "n" not in kmer:
            auto.add_word(kmer, idx)
        else:
            n_count = kmer.count("n")
            for repl in itertools.product("ATCG", repeat=n_count):
                repl_iter = iter(repl)
                expanded = "".join(next(repl_iter) if c == "n" else c for c in kmer)
                auto.add_word(expanded, idx)
    auto.make_automaton()
    return auto

def FindExistedKmers(sequences, k):
    kmers = {sequence[i:i+k] for sequence in sequences for i in range(len(sequence)-k+1)}
    return sorted(kmers, key=atcg_sort_key)

# -----------------------------
# Sequence IO
# -----------------------------
def Get_Seqences(folder, fold_num):
    tp = sorted(f for f in os.listdir(f"{folder}/Kmer/Train") if f.endswith(".tp.fa"))[fold_num]
    tn = sorted(f for f in os.listdir(f"{folder}/Kmer/Train") if f.endswith(".tn.fa"))[fold_num]
    
    def read_fasta(path):
        with open(path) as f:
            return [l.strip().upper() for l in f if not l.startswith(">")]
    
    return read_fasta(f"{folder}/Kmer/Train/{tp}"), read_fasta(f"{folder}/Kmer/Train/{tn}")

# -----------------------------
# Kmer statistics (Optimized)
# -----------------------------
def Auto_CountKmers(auto, pos_seqs, neg_seqs, n_kmers):
    def count(seqs):
        # initiate numpy array: [Seq_Count x Kmer_Count]
        nt = np.zeros((len(seqs), n_kmers), dtype=np.int32)
        t  = np.zeros((len(seqs), n_kmers), dtype=np.int32)
        
        for i, seq in enumerate(seqs):
            for _, idx in auto.iter(seq):
                nt[i, idx] += 1
            for _, idx in auto.iter(RC(seq)):
                t[i, idx] += 1
        
        twostrand = nt + t
        return nt, t, twostrand

    pos_res = count(pos_seqs) # (nt_pos, t_pos, ts_pos)
    neg_res = count(neg_seqs) # (nt_neg, t_neg, ts_neg)
    
    return pos_res, neg_res

def WilcoxonRankSum(kmers, pos_counts, neg_counts):
    nt_pos, t_pos, ts_pos = pos_counts
    nt_neg, t_neg, ts_neg = neg_counts
    
    n_pos, n_kmers = nt_pos.shape
    n_neg = nt_neg.shape[0]
    total_samples = n_pos + n_neg
    total_pairs = float(n_pos * n_neg)
    
    out = []

    def fast_mannwhitney(count_pos, count_neg):
        max_val = int(max(count_pos.max(), count_neg.max())) + 5
        fx = np.bincount(count_pos, minlength=max_val)
        fy = np.bincount(count_neg, minlength=max_val)
        
        cum_fy_counts = np.cumsum(fy)
        total_fy = cum_fy_counts[-1]
        
        strict_smaller_fy = np.zeros_like(fy)
        strict_smaller_fy[1:] = cum_fy_counts[:-1]
        strict_larger_fy = total_fy - cum_fy_counts
        
        n_greater = np.sum(fx * strict_smaller_fy)
        n_lesser  = np.sum(fx * strict_larger_fy)
        n_ties = total_pairs - n_greater - n_lesser
        
        U = n_greater + 0.5 * n_ties
        
        mu_U = total_pairs / 2.0
        
        f_combined = fx + fy
        t_indices = f_combined > 1
        t = f_combined[t_indices]
        if len(t) > 0:
            tie_term = np.sum(t**3 - t) / (total_samples * (total_samples - 1))
        else:
            tie_term = 0
            
        sigma_U = np.sqrt((n_pos * n_neg / 12.0) * (total_samples + 1 - tie_term))
        
        if sigma_U == 0:
            return 0.0, 1.0
        
        z = (U - mu_U) / sigma_U
        p = 2 * norm.sf(abs(z))
        return z, p

    for i, k in enumerate(kmers):
        z_nt, p_nt = fast_mannwhitney(nt_pos[:, i], nt_neg[:, i])
        z_t,  p_t  = fast_mannwhitney(t_pos[:, i],  t_neg[:, i])
        z_ts, p_ts = fast_mannwhitney(ts_pos[:, i], ts_neg[:, i])
        
        tp_sum = ts_pos[:, i].sum()
        tn_sum = ts_neg[:, i].sum()
        
        if tp_sum == tn_sum:
            enriched = "None"
            logfc = 0.0
        elif tp_sum > tn_sum:
            enriched = "TP"
            logfc = np.log2((tp_sum + 1e-6) / (tn_sum + 1e-6))
        else:
            enriched = "TN"
            logfc = np.log2((tp_sum + 1e-6) / (tn_sum + 1e-6))

        out.append({
            "Kmer": k,
            "nt_TPCount": nt_pos[:, i].sum(),
            "t_TPCount":  t_pos[:, i].sum(),
            "2strand_TPCount": tp_sum,
            "nt_TNCount": nt_neg[:, i].sum(),
            "t_TNCount":  t_neg[:, i].sum(),
            "2strand_TNCount": tn_sum,
            "nt_Pvalue": p_nt,
            "t_Pvalue":  p_t,
            "2strand_Pvalue": p_ts,
            "2strand_Z": z_ts,
            "nt_FDR": p_nt, 
            "t_FDR":  p_t,
            "2strand_FDR": p_ts,
            "2strand_Log2FC": logfc,
            "enriched_in": enriched,
        })
        
    return pd.DataFrame(out)

def FDR(df):
    df = df.copy()
    for klen in df["Kmer"].str.len().unique():
        mask = df["Kmer"].str.len() == klen
        sub = df.loc[mask]

        for col in ["nt_FDR", "t_FDR", "2strand_FDR"]:
            _, adj, _, _ = multipletests(sub[col], method="fdr_bh")
            df.loc[mask, col] = adj
    return df

# -----------------------------
# Multiprocessing wrapper
# -----------------------------
def process_large_kmers(kmers, pos, neg):
    # 1. Build Automaton (Map string -> index)
    auto = Build_Automaton(kmers)
    n_kmers = len(kmers)
    
    # 2. Count (Return Numpy Arrays)
    pos_counts, neg_counts = Auto_CountKmers(auto, pos, neg, n_kmers)
    
    if (len(pos_counts) == 0) or (len(neg_counts) == 0):
        raise ValueError("No sequences to count kmers. Please check whether the IDs in the .tp/.tn files match the IDs in the .gff3 file.")

    # 3. Significance Test (Return DataFrame)
    return WilcoxonRankSum(kmers, pos_counts, neg_counts)

# -----------------------------
# Main
# -----------------------------
def main(folder, pThreshold, k1=None, k2=None, gap=None):
    
    for fold in range(5):
        pos, neg = Get_Seqences(folder, fold)

        if all([k1, k2, gap]):
            print(f"Finding Enrich gapped {k1}{'n'*gap}{k2} from fold{fold+1}...")
            kmers = Generate_Gapped_Kmers(k1=k1, k2=k2, gap=gap)
            chunks = [kmers[i::multiprocessing.cpu_count()] for i in range(multiprocessing.cpu_count())]
            with multiprocessing.Pool() as pool:
                results = pool.starmap(process_large_kmers, [(c, pos, neg) for c in chunks])
            df = pd.concat(results)

        else:
            Fivemers = Generate_Kmers(5)
            print(f"Finding Enrich 5mers from fold{fold+1}...")
            chunks = [Fivemers[i::multiprocessing.cpu_count()] for i in range(multiprocessing.cpu_count())]
            with multiprocessing.Pool() as pool:
                results = pool.starmap(process_large_kmers, [(c, pos, neg) for c in chunks])
            df = pd.concat(results)

            for plus in range(6):
                K = 6+plus
                print(f"Finding Enrich {K}mers from fold{fold+1}...")

                kmers = Generate_Kmers(K) if K <= 8 else FindExistedKmers(pos, K)
                print(f"Found {len(kmers)} {K}mers")
                
                if not Kmers:
                    continue

                chunks = [kmers[i::num_processes] for i in range(num_processes)]
                with multiprocessing.Pool() as pool:
                    results = pool.starmap(process_large_kmers, [(c, pos, neg) for c in chunks if c])
                if results:
                    df = pd.concat([df] + results)

        df = FDR(df)
        df = df.sort_values(by="Kmer", key=lambda col: col.map(atcg_sort_key)).reset_index(drop=True)
        
        out = f"{folder}/Kmer/Train/TPTN_{fold+1:02d}.ac_WRSresults.txt"
        df.to_csv(out, sep="\t", index=False)
        print(f"Saved: {out}")