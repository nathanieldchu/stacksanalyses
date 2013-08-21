#!/usr/bin/env python
"""
Created on Sat Jul 13 16:51:17 2013

@author: nathanieldchu
Takes a directory full of STACKS haplotype export files and makes them into phylip files with IUPAC ambiguity code
This assumes that they are all appropriately formatted STACKS haplotype files in tsv format

When run as a bash script, this will require one argument, the path of the directory in which ALL files will be converted to IUPAC Phylip files.
"""
def trimsnplabels(x):
    """(list) -> (list)
    Takes a split list of from the SNPs column and return a list of integers
    >>> trimsnplabels(['40,T>C', '42,C>T', '43,G>C'])
    [40, 42, 43]
    """
    import re
    out = []
    if x == ['']:
        pass
    else:
        for i in x:
            out.append(int(re.sub(r"\D", "", i)))
    return out
    
def LengthIUPAC(radtagsfilename):
    """(file, str) -> int
    Return length of resulting IUPAC sequence for all tags with 3 or less SNPs and 2 pops represented
    """
    import csv
    with open(radtagsfilename, "U") as f:
        radtags = csv.reader(f, delimiter="\t")
        radheader = radtags.next()
        snpindex = radheader.index("Num SNPs")
        consensus = radheader.index("Consensus Sequence")
        count = 0
        for i in radtags:
            if int(i[snpindex]) < 4:
                count += len(i[consensus])
    return count
    
def IUPACfileOnePop(radtagsfilename, pop):
    """(file, str) -> (str, file)
    Convert all rad tags from one population to a fasta format file for all tags with more than 1 population represented
    
    """
    iupac = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'AC': 'M', 'AG': 'R', 'AT': 'W', 'CG': 'S', 'CT': 'Y', 'GT': 'K', 'ACG': 'V', 'ACT': 'H', 'AGT': 'D', 'CGT': 'B', 'ACGT': 'N', 'N': 'N'}
    import csv
    with open(radtagsfilename, "U") as f:
        radtags = csv.reader(f, delimiter="\t")
        radheader = radtags.next()
        snpindex = radheader.index("Num SNPs")
        consensus = radheader.index("Consensus Sequence")
        radindex = radheader.index(pop)
        snplabelsindex = radheader.index("SNPs")
        snps = []
        conseq = []
        radpop = []
        snplabels = []
        for i in radtags:
            if int(i[snpindex]) < 4:
                snps.append(int(i[snpindex]))
                conseq.append(list(i[consensus]))
                radpop.append(i[radindex].split("/"))
                snplabels.append(trimsnplabels(i[snplabelsindex].split(";")))
    seqlist = []
    for j in range(len(radpop)):
        if radpop[j] == ['']:
            seqlist.append("-" * len(conseq[j]))
        elif snps[j] == 0:
            seqlist.append(''.join(conseq[j]))
        else:
            iupacallele = []
            for k in range(snps[j]):
                iupacnuc = iupac[''.join(sorted(list(set(x[k] for x in radpop[j]))))]
                iupacallele.extend(iupacnuc)
            for l in range(snps[j]):
                conseq[j][snplabels[j][l]] = iupacallele[l]
            seqlist.append(''.join(conseq[j]))
    return ''.join(seqlist)

def IUPACFileMultiPopPHYLIP(radtagsfilename, poplist, output):
    with open(output, "w+") as h:
        h.write(str(len(poplist)) + " " + str(LengthIUPAC(radtagsfilename)) + "\n")
        for i, j in enumerate(poplist):
            sample = j[(j.find("sample_") + 7):]
            label = sample + (" " * (10 - len(sample)))
            h.write(label + IUPACfileOnePop(radtagsfilename, j))
            if i < len(poplist) - 1:
                h.write("\n")

if __name__ == "__main__":
    import sys, os
    poplist = ["sample_CAGT.WHD", "sample_CAGT.QHC", "sample_ACGT.ANP", "sample_GGTT.DMC", "sample_ACGT.CBC", "sample_AATT.BVR"]
    os.chdir(sys.argv[1])
    for filename in os.listdir(sys.argv[1]):
        IUPACFileMultiPopPHYLIP(filename, poplist, filename + "Phylipfile.phy") 
