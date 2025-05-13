import os
import Translation

def rmlb(astr):
    if astr[-2:] == "\r\n":
        astr = astr[:-2]
    elif astr[-1] == "\n":
        astr = astr[:-1]
    return astr

def fasta_to_dict(fasta,dflag=0,verbose=0,newline=0):

    # each idx should only occur once.
    inp = open(fasta,"r")
    inl = inp.readline()
    fdict = {} # idx as key, seq as value
    c = 0
    N = 0 
    while inl != "":
        inl = rmlb(inl)
        next = 0
        idx  = ""
        desc = ""
        if inl == "":
            pass
        elif inl[0] == ">":
            if verbose and c%1e3 == 0:
                print (" %i k" % (c/1e3))
            c += 1
            # rid of anything after space if asked
            if dflag and inl.find(" ") != -1:
                desc = idx[idx.find(" ")+1:]
                idx = idx[:idx.find(" ")]
            else:
                idx = inl[1:]
            # count lines and store seq into a list
            seq_list = []
            inl = inp.readline()
            while inl[0] != ">":
                inl = inl.strip()
                if newline:
                    inl = inl + "\n"
                seq_list.append(inl)
                inl = inp.readline()
                if inl == "":
                    break
            seq = "".join(seq_list)
            
            if idx in fdict.keys():
            #if fdict.has_key(idx):
                if verbose:
                    print ("Redundant_id:",idx,)
                if dflag:
                    if len(fdict[idx][1]) < len(seq):
                        fdict[idx] = [desc,seq]
                        if verbose:
                            print ("longer")
                    else:
                        if verbose:
                            print ("shorter")
                else:
                    if len(fdict[idx]) < len(seq):
                        fdict[idx] = seq
                        if verbose:
                            print ("longer")
                    else:
                        if verbose:
                            print ("shorter")
            else:
                N += 1
                if dflag:
                    fdict[idx] = [desc,seq]
                else:
                    fdict[idx] = seq
            next = 1
        
        # so no extra line is read, because of the innder while
        if not next:
            inl = inp.readline()
    inp.close()
    if verbose:
        print ("Total %i sequences, %i with non-redun names" % (c,N))
    
    
    return fdict

def get_stretch4(fasta,coords,seqid=0):
    
    print ("Sequence to dict...")
    seq = fasta_to_dict(fasta,0)
    # fasta_to_dict got rid of "\n" already

    new_dict={}
    for i in seq.keys():
        newseq= i.strip().split(" ")[0]
        value= seq[i]
        new_dict[newseq]=value

    c = 0 # count total
    m = 0 # count not in fasta
    if coords.find(",") == -1:
        print ("Read coordinates...")
        inp = open(coords)
        oup = open(coords+".fa","w")
        oup2= open(coords+".missing","w")
        inl = inp.readline()
        while inl != "":    # Go through each coord
            # if c % 1000 == 0:
            #     print (" %i k" % (c/1000))
            c  += 1
            L  = inl.strip().split("\t")
            seqName = L[0]	# Sequnece name
            if seqName in new_dict:
                if len(L) >= 2:
                    # L = [name, L, R], some may have 4th col which is IDs to be given.
                    
                    # if L[1].find(",") == -1:
                    #     # Deal with reverse ori
                    #     ori= 1; cL = int(L[1]); cR = int(L[2])
                    #     if cL > cR:
                    #         ori = -1
                    #     # Get sequence
                    #     if ori == -1:
                    #         S = new_dict[seqName][cR-1:cL]
                    #     else:
                    #         S = new_dict[seqName][cL-1:cR]

                    #     if S == "":
                    #         print ("ERR COORD: %s,[%i,%i]" % (seqName,cL,cR))
                    #     else:
                    #         if ori == -1:
                    #             S = trans.rc(S)
                            
                    #         # If there is 4th column, use them as sequence IDs.
                    #         if len(L) == 4 and seqid:
                    #             oup.write(">%s\n%s\n" % (L[3],S))
                    #         else:
                    #             oup.write(">%s|%i-%i\n%s\n" %	(seqName,cL,cR,S))
                    # name <\t> "L1,R1,L2,R2..." <\t> whatever
                    if L[1].find(",") != -1:
                        genename= L[2]
                        coordList = L[1].split(",")
                        S = ""
                        # Set orientation, only consider the first pair
                        if int(coordList[0]) < int(coordList[1]):
                            ori = 1
                        else:
                            ori = -1
                        for j in range(0,len(coordList),2):
                            cL = int(coordList[j]); cR = int(coordList[j+1])
                            if ori == -1:
                                cL = int(coordList[j+1])
                                cR = int(coordList[j])
                                S = new_dict[seqName][cL-1:cR] + S
                                S = RC(S)
                            else:
                                S += new_dict[seqName][cL-1:cR]
                        oup.write(">%s %s|%s\n%s\n" % (genename, seqName,"-".join(coordList),S))
                    else:
                        print ("Unknown cooord format:",L)
                        print ("Quit!")
                        sys.exit(0)

            else:
                m += 1
                oup2.write(inl)
            inl = inp.readline()
        print ("Total coords:",c)
        print ("Not in seq  :",m)
        oup2.close()

    # coordinates are passed
    else:
        coords = coords.split(",")
        print ("Coords:",coords)
        oup = open("%s_%s.fa" % (fasta,"-".join(coords)),"w")
        C = []
        for i in coords:
            C.append(int(i))
        for i in seq:
            s = ""
            for j in range(0,len(C),2):
                s += seq[i][C[j]-1:C[j+1]]
            oup.write(">%s\n%s\n" % (i,s))
                
    oup.close()
    print ("Done!")	

def get_genome_file_path(new_folder):
    cwd = os.getcwd()
    print(f"Coord_to_Fasta working direction: {cwd}")
    Gff_and_Genome_folder = f"{new_folder}/Gff_and_Genome"
    Gff_and_Genome_files = sorted(os.listdir(Gff_and_Genome_folder))
    genome_file = ""
    for file in Gff_and_Genome_files:
        if not file.endswith(".gff") and not file.endswith(".gff3") and not file.endswith(".coord") and not file.endswith(".coord.fa") and not file.endswith(".missing"):
            genome_file = os.path.join(cwd, Gff_and_Genome_folder, file)
            break
        else:
            continue
    
    if genome_file == "":
        raise NameError(f"Can't get genome file name. Make sure the folder({Gff_and_Genome_folder}) have only 1 .gff file & 1 genome file(file which not end with '.gff' will be considered as a genome file.")
    else:
        print(f"genome_file: {genome_file}")
    return genome_file

def get_coords_file_path(new_folder):
    cwd = os.getcwd()
    Gff_and_Genome_folder = f"{new_folder}/Gff_and_Genome"
    Gff_and_Genome_files = sorted(os.listdir(Gff_and_Genome_folder))
    coords_file = ""
    for file in Gff_and_Genome_files:
        if file.endswith(".coord"):
            coords_file = os.path.join(cwd, Gff_and_Genome_folder, file)
            break
        else:
            continue
    
    if coords_file == "":
        raise NameError(f"Can't get coords file name. The '.coord' file is automatically generated by 'Gff_to_Coord.py'")
    else:
        print(f"coords_file: {coords_file}")
    return coords_file

def RC(kmer):
    complement = str.maketrans("ACGTacgt", "TGCAtgca")  # DNA complement mapping
    return kmer.translate(complement)[::-1]

def main(new_folder):
    trans = Translation.translate()
    fasta = get_genome_file_path(new_folder)
    coords = get_coords_file_path(new_folder)
    get_stretch4(fasta=fasta, coords=coords)