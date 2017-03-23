import pandas as pd
import numpy as np
import re

#reverse compliment of dna string
def revComp(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in dna[::-1]])

"""This class turns all of the lines into lists, separated by tabs. the np.vstack is used
to create a matrix. this will be used later for selecting the column desired.dimesnion is
len(header)x Number of variants"""

class brcaParse:
    def __init__(self, input, output1, output2):
        f = open(input, "r")
        #f_out = open(output1, "w")
        rawLabels = f.readline()
        labels = rawLabels.split("\t")

        #imports the sequence for BRCA1/2 into the program for use.
        self.BRCA1hg38Seq = open("brca1_hg38.txt", "r").read()
        self.BRCA1hg38Start = 43000000

        self.BRCA2hg38Seq = open("brca2_hg38.txt", "r").read()
        self.BRCA2hg38Start = 32300000

        self.BRCA1 = {"hg38": {"start": 43000000, "sequence": open("brca1_hg38.txt", "r").read()},
                      "hg19": {"start": 41100000, "sequence": open("brca1_hg19.txt", "r").read()}}
        self.BRCA2 = {"hg38": {"start": 32300000, "sequence": open("brca2_hg38.txt", "r").read()},
                      "hg19": {"start": 32800000, "sequence": open("brca2_hg19.txt", "r").read()}}

        # finds column in matrix associated with header label
        a = labels.index("Alt")
        b = labels.index("Ref")
        c = labels.index("Pos")
        d = labels.index("Clinical_significance_ENIGMA")
        e = labels.index("Gene_Symbol")
        g = labels.index("id")
        #print(labels)
        buildMat = []
        buildMat.append(labels)
        #f_out.write(rawLabels)

        # for loop used to append a complete data set associated with an id number. Later put into matrix/array.
        #A is all of the input data as a matrix. matOut is the selected columns as lists.
        for lines in f:
            l = f.readline()
            if len(l.split("\t")) == len(labels):
                buildMat.append(l.split('\t'))
            #f_out.write(l)
        L = np.vstack(buildMat)
        self.A = L
        #f_out.close()

        # use column definition to get list of a, b, c, and d. Now this is modification, reference, position and significance.
        self.Alt = brcaParse.column(self.A, a)  # positon of mutation column
        self.Ref = brcaParse.column(self.A, b)  # reference sequence column
        self.Pos = brcaParse.column(self.A, c)  # alteration of sequence column
        self.Sig = brcaParse.column(self.A, d)  # Clinical significance column
        self.Gene = brcaParse.column(self.A, e)  # BRCA1/2 column for direction of sequence
        self.id = brcaParse.column(self.A, g) #id number column
        #print("{0}\n{1}\n{2}\n{3}\n{4}".format(self.Alt, self.Ref, self.Pos, self.Sig, self.Gene))
        self.matOut = np.vstack((self.Alt, self.Ref, self.Pos, self.Sig, self.Gene))

        #f_out = open(output2, "w")
        #f_out.write(str(self.matOut))
        #f_out.close()

    #gets the desired data and puts to object
    def getDat(self):
        return (self.A, self.matOut)

    #Creates a column from the matrix of data.
    def column(matrix, i):
        return [row[i] for row in matrix]


    #Parses through the string to make all the new variant strings. If the variant is an indel
    #rather than a SNP, the definition will make more iterations to account for the new modifications.
    def maxEntForm(self,output, lengthSeq):
        f = open(output, 'w')
        lenSplice = lengthSeq

        for i in range(0,len(self.Gene)):
            if self.Gene[i] == "BRCA1":
                f.write("> " + self.id[i] +"\n" )
                loc = (int(self.Pos[i]) - int(self.BRCA1hg38Start))
                tempSeq = self.BRCA1hg38Seq[:loc-1] + self.Alt[i] + self.BRCA1hg38Seq[loc+len(self.Ref[i])-1:]
                for j in range(0,lenSplice- 1 + (len(self.Ref[i]))):
                    n = loc-lenSplice+j
                    o = loc+j
                    #print(tempSeq[n:o])
                    f.write(self.BRCA1hg38Seq[n:o] + "\n")
                    f.write(tempSeq[n:o] + "\n")

            if self.Gene[i] == "BRCA2":#fix for reverse compliment and such
                f.write("> " + self.id[i] + "\n")
                loc = (int(self.Pos[i]) - self.BRCA2hg38Start)
                tempSeq = self.BRCA2hg38Seq[:loc - 1] + self.Alt[i] + self.BRCA2hg38Seq[loc + len(self.Ref[i]) - 1:]
                loc2 = len(self.BRCA2hg38Seq) - loc
                tempSeqComp = revComp(tempSeq)
                BRCA2Comp = revComp(self.BRCA2hg38Seq)
                for j in range(0, lenSplice -1 + (len(self.Ref[i]))):
                    n = loc2 - lenSplice + j
                    o = loc2 + j
                    f.write(BRCA2Comp[n:o] + "\n")
                    f.write(tempSeqComp[n:o] + "\n")

    #takes the raw scores in str format from MaxEntScan and puts them through a few parsing loops.
    #the data is eventually transformed into something in float/int format for use in other programs.
    def getMaxEntMax(self, input):
        f = open(input, 'r')
        rawScoresStr = []
        self.intScore = []

        for line in f:
            m = re.findall("[-+]?\d+[\.]?[-+]?\d*", line)
            rawScoresStr.append(m)

        for i in range(0,len(rawScoresStr)-1):
            try:
                if type(int(rawScoresStr[i][0]))==int:
                    self.intScore.append(int(rawScoresStr[i][0]))
            except ValueError:
                if type(float(rawScoresStr[i][0]))==float:
                    self.intScore.append(float(rawScoresStr[i][0]))
        print(self.intScore)

def main():
    a = brcaParse("variants.tsv", "1", "2")
    A, matOut = a.getDat()
    maxEnt = 'maxEntScanList'
    a.maxEntForm(maxEnt,9)#get donor splice sites 9mers
    #a.maxEntForm(maxEnt, 23) #get acceptor splice sites 23mers UNCOMMENT FOR 23 mer!
    #a.getMaxEntMax('maxEntOut') #get list of max ent scores


# Make loops. if BRCA1 or 2, reference BRCA1 or BRCA1, then replace string with alt(eration).

if __name__ == "__main__":
    main()

