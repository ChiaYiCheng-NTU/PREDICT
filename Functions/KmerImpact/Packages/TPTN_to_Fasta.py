import os

def rmlb(astr):
    if astr[-2:] == "\r\n":
        astr = astr[:-2]
    elif astr[-1] == "\n":
        astr = astr[:-1]
    return astr

def getseq(fasta,name,tokens="",verbose=0):
    
    if tokens == 0 or tokens == "":
        tokens = [0]
    else:
        tmp = tokens.split(",")
        tokens = []
        for i in tmp:
            tokens.append(int(i))
    
    inp = open(name)
    inl = inp.readline()
    ndict = {}
    countR = 0
    rlist  = []
    while inl != "":    
        L = inl.strip().split("\t")
        for j in tokens:
            if L[j] == "":
                pass
            elif L[j] not in ndict:
                ndict[L[j]] = 0
            else:
                countR += 1
                rlist.append(L[j])
        inl = inp.readline()
    countN = len(ndict.keys())
    oup = open(name+".fa","w")
    countA = 0
    countF = 0
    for i in fasta.split(","):
        if verbose:
            print ("",i)
        inp = open(i)
        inl = inp.readline()
        flagW  = 0
        while inl != "":
            if ">" in inl:
                countA += 1
                if countF == countN:
                    break
                N = rmlb(inl)[1:]
                N = N.split(" ")[0]
                if N in ndict and ndict[N] == 0:
                    ndict[N] = 1
                    if verbose:
                        if countF % 1e2 == 0:
                            print (" %i x100" % (countF/1e2))
                    countF += 1
                    flagW = 1
                else:
                    flagW = 0
            if flagW:
                oup.write(inl.upper())

            inl = inp.readline()
            
    oup.close()

    missing = []
    if countF < len(ndict.keys()):
        for i in ndict:
            if ndict[i] == 0:
                missing.append(i)

    print (" Redun:%i, Unique:%i, Found:%i, Missed:%i" % \
        (countR,countN,countF,len(missing)))

    if len(missing) > 0:
        print (f"You have {len(missing)} gene(s) missing in {len(ndict.keys())} input ({len(missing)/len(ndict.keys())*100:.2f}%)")
    
    if countF == 0:
        print("******************* ERROR *******************")
        print("******************* ERROR *******************")
        print("No sequences were found in the input files.")
        print("It is likely that the identifiers in the .tp/.tn files do not match any of the identifiers in the .gff3 file.")
        print("***************** ERROR END *****************")
        
    return missing

def main(new_folder, For_BedData):
    if For_BedData == False:
        fasta_folder = f"{new_folder}/Gff_and_Genome"
        fasta_name = [fasta for fasta in sorted(os.listdir(f"{new_folder}/Gff_and_Genome")) if fasta.endswith(".coord.fa")][0]
        fasta_path = os.path.join(fasta_folder, fasta_name)
        print(f"Fasta_source_path: {fasta_path}")
        for folder_path in [f"{new_folder}/Kmer/Test", f"{new_folder}/Kmer/Train"]:
                all_files = sorted(os.listdir(folder_path))
                tptn_files = [file for file in all_files if file.endswith(".tp") or file.endswith(".tn")]
                for file_name in tptn_files:
                    file_path = os.path.join(folder_path, file_name)
                    getseq(fasta=fasta_path, name=file_path)
    
    if For_BedData == True:
        fasta_folder = f"{new_folder}/Gff_and_Genome"
        fasta_name = [fasta for fasta in sorted(os.listdir(f"{new_folder}/Gff_and_Genome")) if fasta.endswith(".coord.fa")][0]
        fasta_path = os.path.join(fasta_folder, fasta_name)
        print(f"Fasta_source_path: {fasta_path}")
        for folder_path in [f"{new_folder}/Kmer/Train"]:
                all_files = sorted(os.listdir(folder_path))
                tptn_files = [file for file in all_files if file.endswith(".Btp") or file.endswith(".Btn")]
                for file_name in tptn_files:
                    file_path = os.path.join(folder_path, file_name)
                    getseq(fasta=fasta_path, name=file_path)