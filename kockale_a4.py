#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 26 09:07:44 2021

@author: logankocka
"""
import pandas as pd
import sys
# sys.path.append("/Users/krunal/Desktop/CSE 466/Proj04")
import SmithWaterman_class as sw


class MyGene:

    def __init__(self, File1Name):
            self.FastaFileName = File1Name
            self.SeqDict = {} # A dictionary for gene ID and sequence  
            self.SppDict = {} # A dictionary for gene ID and names of gene and species
            self.bs_list = ['TATTGTTTATT',
                            'TATTGTTTATT',
                            'TATTGTTTACT',
                            'AATTGTTTATT',
                            'AATTGTTTATT',
                            'CATTGTTTATT',
                            'TATTGTTTATT',
                            'GATTGTTTACT',
                            'AATTGTTTATT',
                            'AATTGTTTAGT'] # a list of default binding site
            self.Consensus = ''
            self.IUPAC = ''
            self.PFM = pd.DataFrame()
            self.PPM = pd.DataFrame()
            self.__extract_file()
    
    
    def __extract_file(self):
       fh = open(self.FastaFileName, 'r')
       name, desc, seq=('', '', '')
       for line in fh:
           line = line.strip()

           if line.startswith('>'):
               if seq !='':
                   name, desc, seq = ('', '', '')
               name = line[1:line.index(' ')]
               desc = line[line.index(' ')+1:]
               self.SppDict[name] = desc #SppDict is key=name, value=desc
           else:
               seq=seq+line
               self.SeqDict[name] = seq #SeqDict is key=name, value=seq
               
    # >GD0001 (name) ABC-0 Species01 (desc)
    # TGTCCAACGGGCCGAGGTTGTCTCTTTCGAGATCTTGTCGCGGGGGGGGGCTGCCTGTGGCGGTGGGTGG
    # AGTGCGGGTCACAGATGCCTCGCACCTATTGTTTATTCCTCCCCCGACAATGTGGCCCGTATGGAGGGTCTCG (seq)
    
    
    def set_motif(self, bindSites):
        self.bs_list = bindSites
        print("[1]. The protein binding sequences are: ")
        for i, seq in enumerate(self.bs_list):
            print(i+1, seq)

    
    def set_Consensus(self, bindSites):
        n = len(bindSites[0])
        prof = { 'T':[0]*n,'G':[0]*n ,'C':[0]*n,'A':[0]*n } 
        
        for seq in bindSites:
            for i, a in enumerate(seq):
                prof[a][i] += 1
        cons = ""
        for i in range(n):
            max_count = 0
            max_nt = 'x'
            for nt in "ACGT":
                if prof[nt][i] > max_count:
                    max_count = prof[nt][i]
                    max_nt = nt
            cons += max_nt
        self.Consensus = cons
        return self.Consensus
       
   
    def set_IUPAC(self):
        # IDK the rules for this so im returning the highest probability nucleotide
        # if there are different rule this will be wrong
        
        A = []; T = []; G = []; C = []; # initialize lists
        A_count = G_count = C_count = T_count = 0

        temp_list = []
        new_dict = {}
        for seq in self.bs_list: # this loops thru the sequences
            for char in seq: # this loops thru the chars in a sequence
                temp_list.append(char)
            new_dict[seq] = temp_list
            temp_list = []
       
        values_only = new_dict.values() # get list of values for dataframe
        df = pd.DataFrame(values_only)
       
        for col in df.columns:
            # A_count = df[col].value_counts()["A"]
            A_count = len(df[df[col] == 'A'])
            A.append(A_count)
            T_count = len(df[df[col] == 'T'])
            T.append(T_count)
            G_count = len(df[df[col] == 'G'])
            G.append(G_count)
            C_count = len(df[df[col] == 'C'])
            C.append(C_count)
            # reset back to 0
            A_count = G_count = C_count = T_count = 0 
       
        list_new = [A, T, G, C]
        df_new = pd.DataFrame(list_new)
        df_new.index = ["A", "T", "G", "C"]
    
        n = sum(df_new[:][0])
        PPM = df_new / int(n)
      
        front = PPM.idxmax()
        back = PPM.iloc[::-1].idxmax()
        
        front_maxes = front.values.tolist()
        back_maxes = back.values.tolist()
        maxes = [0] * len(front_maxes)
        
        for i in range(len(front_maxes)):
            # print("front", front_maxes[i])
            # print("back", back_maxes[i])
            if front_maxes[i] == back_maxes[i]:
                maxes[i] = front_maxes[i]
            else:
                if front_maxes[i] == 'A' and back_maxes[i] == 'G': maxes[i] = 'R'
                elif front_maxes[i] == 'T' and back_maxes[i] == 'C': maxes[i] = 'Y'
                elif front_maxes[i] == 'G' and back_maxes[i] == 'C': maxes[i] = 'S'
                elif front_maxes[i] == 'A' and back_maxes[i] == 'T': maxes[i] = 'W'
                elif front_maxes[i] == 'T' and back_maxes[i] == 'G': maxes[i] = 'K'
                elif front_maxes[i] == 'A' and back_maxes[i] == 'C': maxes[i] = 'M'
                else: maxes[i] = 'N'
        
        self.IUPAC = "".join(maxes)
        print(self.IUPAC)
       
       
    def set_PFM(self): # find the position frequency matrix
        #loop thru each position of each item in the list
        A = []; T = []; G = []; C = []; # initialize lists
        A_count = G_count = C_count = T_count = 0

        temp_list = []
        new_dict = {}
        for seq in self.bs_list: # this loops thru the sequences
            for char in seq: # this loops thru the chars in a sequence
                temp_list.append(char)
            new_dict[seq] = temp_list
            temp_list = []
        
        values_only = new_dict.values() # get list of values for dataframe
        df = pd.DataFrame(values_only)
        
        for col in df.columns:
            # A_count = df[col].value_counts()["A"]
            A_count = len(df[df[col] == 'A'])
            A.append(A_count)
            T_count = len(df[df[col] == 'T'])
            T.append(T_count)
            G_count = len(df[df[col] == 'G'])
            G.append(G_count)
            C_count = len(df[df[col] == 'C'])
            C.append(C_count)
            # reset back to 0
            A_count = G_count = C_count = T_count = 0 
        
        list_new = [A, T, G, C]
        df_new = pd.DataFrame(list_new)
        df_new.index = ["A", "T", "G", "C"]
        # df_new.columns = [''] * len(df_new.columns)
        
        self.PFM = df_new
        print(self.PFM)
    
    
    def set_PPM(self):
        df = self.PFM 
        n = sum(df[:][0])
        self.PPM = df / n
        print(self.PPM)
    
    
    def find_bs(self): 
        # loop thru the given dictionaries to extract one seq, gene desc, and species name
        for name in self.SeqDict.keys():
            seq = self.SeqDict[name]
            desc = self.SppDict[name]
            gene_desc = desc.split(" ")[0]
            species = desc.split(" ")[1]
            start, end = self.scan_seq(seq)
            print("\nSequence Num:\t{0}".format(species[-2:]))
            print("GeneID:\t\t{0}".format(name))
            print("Gene Name:\t{0}".format(gene_desc))
            print("Species Name:\t{0}".format(species))
            
            if start >= 2:
               ext_start = start - 2
            else:
                ext_start = start
            if end != len(seq):
                ext_end = end + 2
            else:
                ext_end = end
            ext_seq = seq[ext_start:ext_end]
            
            aligned_seq1, aligned_seq2, identity_perc, gap_perc = self.align_seq(ext_seq, self.Consensus)
            print("Identity:\t{0:.1%}".format(identity_perc))
            print("Gaps:\t\t{0:.1%}\n".format(gap_perc))
            
            # print(type(start))
            # print(type(name))
            # print(type(aligned_seq1))
            self.get_alignment(start, (start+len(aligned_seq1)), name, aligned_seq1, aligned_seq2)
            seq_print = seq[:start] + seq[start:end].lower() + seq[:end]
            printWithRuler(name, seq, gene_desc, "Y", "")
        
    def scan_seq(self, seq):
        # returns positions of seq input that are most likely to be the binding site
      
        d = {}
        seqLen = len(seq)
        bsLen = len(self.bs_list[0])
        for i in range(0, seqLen - bsLen):
            seq_fragment = seq[i:i+bsLen]
            frag_prob = 0
            for j in range(0, bsLen):
                nucleo = seq_fragment[j]
                ppm = self.PPM.iloc[:,j]
                nuc_prob = ppm[nucleo]
                frag_prob += nuc_prob
            d[i] = frag_prob
        maxS = max(d, key=d.get)
       
        return maxS, maxS+bsLen
        
          
    def align_seq(self, ext_seq, motif): 
        # this uses the smith waterman algorithm
        aligned_seq1, aligned_seq2, identity_perc, gap_perc = sw.Smith_Waterman(ext_seq, motif).give_final_result()
        
        return aligned_seq1, aligned_seq2, identity_perc, gap_perc
    
    
    def get_alignment(self, position1, position2, gene_id, gene_seq, bs_seq):
        align = ""
        for i in range(0,len(bs_seq),1):
            if gene_seq[i] == bs_seq[i]:
                align += "|"
            elif gene_seq[i] == "-" or bs_seq[i] == "-":
                align += " "
            else:
                align += ":"
        print("{0}\t\t{1}\t{2}\t{3}".format(gene_id, position1, gene_seq, bs_seq))
        print("\t\t\t\t{0}".format(align))
        print("bs_site\t\t1\t{0}\t{1}\n".format(bs_seq, (position2-position1)))
        return
        
    
def printWithRuler(seqName, sequence, seqDescription, spacerChoice, indexes):
    if indexes == "":
        NON_FASTA_LINE_NUM = 10
        if spacerChoice == "Y":
            spacer = " "
        if spacerChoice == "N":
            spacer = ""
        
        # print first line and numbers with spaces
        print('>' + seqName + ' ' + seqDescription + '\n')
        print(15*' ', end="")
        nums = list(range(1,11))
        for i in nums:
            print(nums[i-1], end="")
            if i < 9:
                print(9*' ' + spacer, end="")
            if i == 9:
                print(8*' ' + spacer, end="")
        print('\n' + " Line", end=" ")
        for i in nums:
            print("1234567890" + spacer, end="")
            if i == 10: print()
        
        # print each line with spacer every ten characters
        if spacerChoice == "N":
            # format the list
            NON_FASTA_LINE_NUM = 100 # split every 100 characters
            l_N = [sequence[i:i+NON_FASTA_LINE_NUM] for i in range(1, len(sequence), NON_FASTA_LINE_NUM)]
            #print them
            for count,val in enumerate(l_N):
                print('{:>5} {:<100}'.format(count+1,val))
        
        if spacerChoice == "Y":
            # format the list
            l_Y = [sequence[i:i+10] for i in range(1, len(sequence), 10)]
            l_Y2 = []; start = 0
            for i in l_Y:
                str = " ".join(l_Y[start:start+10])
                l_Y2.append(str)
                start += 10
            l_Y3 = list(filter(None, l_Y2))
            # print them
            for count,val in enumerate(l_Y3):
                print('{:>5} {:<110}'.format(count+1,val))

    else:
        NON_FASTA_LINE_NUM = 10
        if spacerChoice == "Y":
            spacer = " "
        if spacerChoice == "N":
            spacer = ""
        
        # print first line and numbers with spaces
        print('>' + seqName + ' ' + seqDescription + '\n')
        print(15*' ', end="")
        nums = list(range(1,11))
        for i in nums:
            print(nums[i-1], end="")
            if i < 9:
                print(9*' ' + spacer, end="")
            if i == 9:
                print(8*' ' + spacer, end="")
        print('\n' + " Line", end=" ")
        for i in nums:
            print("1234567890" + spacer, end="")
            if i == 10: print()
        
        # find and replace the motifs with lower case letters
        seq_list = list(sequence)
        for i in indexes:
            # for j in range(indexes[i][0],indexes[i][-1]):
            
            seq_list[i[0]:i[1]] = [j.lower() for j in seq_list[i[0]:i[1]]]
            
        print("TEST", ''.join(seq_list))
        sequence = ''.join(seq_list)
        # print each line with spacer every ten characters
        if spacerChoice == "N":
            # format the list
            NON_FASTA_LINE_NUM = 100 # split every 100 characters
            l_N = [sequence[i:i+NON_FASTA_LINE_NUM] for i in range(1, len(sequence), NON_FASTA_LINE_NUM)]
            #print them
            for count,val in enumerate(l_N):
                print('{:>5} {:<100}'.format(count+1,val))
        
        if spacerChoice == "Y":
            # format the list
            l_Y = [sequence[i:i+10] for i in range(1, len(sequence), 10)]
            l_Y2 = []; start = 0
            for i in l_Y:
                str = " ".join(l_Y[start:start+10])
                l_Y2.append(str)
                start += 10
            l_Y3 = list(filter(None, l_Y2))
            # print them
            for count,val in enumerate(l_Y3):
                print('{:>5} {:<110}'.format(count+1,val))

    return


       
def main():
    # print the main welcome message
    stars = "********************************************************************************************"
    msg = "Use this program to detect to detect potential protein binding sites for a set of sequences."
    print(stars, '\n', msg, '\n', stars, '\n', sep="")
    # get inputs 
    args = sys.argv
    program = args[0]
    user = program.split("_")[0]
    file = args[1]
    # define the thing
    mygene = MyGene(file)
    numSeqs =  len(mygene.SeqDict.keys()) # get number of sequences
    # print the program info
    print("Program: ", program)
    print("Developer: ", user)
    print("Input File: ", file, '\t', "(There are a total of ", numSeqs, " sequences)")
    
    # set binding site motifs
    bindSite_file = input("Please enter the name of the file that contains binding-site sequences (e.g., bs.txt) ")
    bindSites = [] # get the bind sites list
    if bindSite_file:
        if bindSite_file == "bs.txt":
            fh = open(bindSite_file, 'r')
            for line in fh:
                bindSites.append(line.strip())
        else: bindSites = mygene.bs_list
    else: bindSites = mygene.bs_list
    # call set_motif
    mygene.set_motif(bindSites) # this prints the protein binding sequences to main
    consensus = mygene.set_Consensus(bindSites)
    print("\nThe consensus sequence is", consensus)
    # print("why is this happening")
    print("The IUPAC sequence is ", end="")
    mygene.set_IUPAC()
    # print matrices
    print("\n[2]. The Position Frequency Matrix for the protein binding sequences is:")
    mygene.set_PFM(); print()
    print("\n[3]. The Position Probability Matrix for the protein binding sequences is:")
    mygene.set_PPM(); print();
    
    print("\n[4]. Results of sequence scan")
    mygene.find_bs()
    
    

if __name__ == '__main__':
  main()
  
  