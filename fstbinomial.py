#!/usr/bin/env python
"""
Created on Wed Mar 20 16:40:08 2013

@author: nathanieldavidchu
Take a STACKS haplotype export and allele depth exports calculates pairwise FST values for each pair of populations in poplist using a binomial approach as in Hohenlohe 2010.

Usage:
./fstbinomial.py -i <inputhaplotypefile> -h <inputalleledepthsfile> -o <outputfile> -p <poplist>

where <poplist> is a string of the populations to be compared
"pop1 pop2 pop3 ..."
"""
def binomial(n, k):
   """Compute n factorial by a direct multiplicative method."""
   if k > n - k: 
	k = n - k
   accum = 1
   for i in range(1, k + 1):
      accum *= (n - (k - i))
      accum /= i
   return accum
   
def nucleodiv(x):
    """Takes a list of tuples containing the allele and read count and returns pi"""
    numer = 0
    totalreads = sum([y[1] for y in x])
    deno = binomial(totalreads, 2)
    for i in x:
        numer += binomial(i[1], 2)
    return 1 - (numer/deno)
    
def mergetuplelist(x, y):
    #Return a merged tuple list retaining unique values and summing read counts
    merged = x + y
    alleles = list(set([z[0] for z in x] + [w[0] for w in y]))
    counts = []
    for i in alleles:
        ind = [z for z in range(len(merged)) if merged[z][0] == i]
        total = 0
        for j in ind:
            total += merged[j][1]
        counts.append(total)
    return zip(alleles, counts)

def FstRADHohen(radtagsfilename, haplotagsfilename, pop1, pop2):
    #Returns SNP identities for 2 given pooled poputlations
    import csv
    count = 0
    markersused = 0
	#Isolate all tag identities
    with open(radtagsfilename, "U") as f:
        radtags = csv.reader(f, delimiter="\t")
        radheader = radtags.next()
        radindex1 = radheader.index(pop1)
        radindex2 = radheader.index(pop2)
        snpindex = radheader.index("Num SNPs")
        radpopA = []
        radpopB = []
        snps = []
        for i in radtags:
            if i[radindex1] <> "" and i[radindex2] <> "":
                radpopA.append(i[radindex1].split("/"))
                radpopB.append(i[radindex2].split("/"))
                snps.append(int(i[snpindex]))
	#Isolate allele frequencies
    with open(haplotagsfilename, "r+") as g:
        haplotags = csv.reader(g, delimiter = "\t")
        haploheader = haplotags.next()
        haploindex1 = haploheader.index(pop1)
        haploindex2 = haploheader.index(pop2)
        happopA = []
        happopB = []
        for j in haplotags:
            if j[haploindex1] <> "" and j[haploindex2] <> "":
                happopA.append(map(int, j[haploindex1].split("/")))
                happopB.append(map(int, j[haploindex2].split("/")))
	#Calculate FST for all markers that are not missing data for either population
    for k in range(len(radpopA)):
        if radpopA[k] == radpopB[k] and len(radpopA[k]) == 1 and len(radpopB[k]) == 1:
            markersused += 1
        elif radpopA[k] != radpopB[k] and len(radpopA[k]) == 1 and len(radpopB[k]) == 1:
            count += 1
            markersused += 1
		#Filter out stacks with too many polymorphic sites (cutoff at 3)
        elif snps[k] < 4:
			#Zip together each individual SNP position with its allele frequency
            combA = sorted(zip(radpopA[k], happopA[k]))
            combB = sorted(zip(radpopB[k], happopB[k]))
            radpopAset = set(radpopA[k])
            radpopBset = set(radpopB[k])
			#Check that each zipped allele list has no redundancies
            if len(radpopAset) == len(radpopA[k]):
                pass
            else:
                haplodepths = []
                for l in radpopAset:
                    indicies = [m for m, x in enumerate(radpopA[k]) if x == l]
                    readcount = 0
                    for n in indicies:
                        readcount += happopA[k][n]
                    haplodepths.append(readcount)
                combA = zip(radpopAset, haplodepths)
            if len(radpopBset) == len(radpopB[k]):
                pass
            else:
                haplodepths = []
                for l in radpopBset:
                    indicies = [m for m, x in enumerate(radpopB[k]) if x == l]
                    readcount = 0
                    for n in indicies:
                        readcount += happopB[k][n]
                    haplodepths.append(readcount)
                combB = zip(radpopBset, haplodepths)
			#Sum total read counts and calculate ratios for each allele from each population.
            merged = mergetuplelist(combA, combB)
            sumallelA = len(combA)
            sumallelB = len(combB)
            numerator = ((binomial(sumallelA, 2) * nucleodiv(combA)) + (binomial(sumallelB, 2) * nucleodiv(combB)))
            denomenator = nucleodiv(merged) * (binomial(sumallelA, 2) + binomial(sumallelB, 2))
			#Add up Fst values per loci and keep track of how many have been analysed (markersused)
            if denomenator == 0:
                markersused += 1
            elif numerator/denomenator > 1:
                pass
            else:
                count += 1 - numerator/denomenator
                markersused += 1
	#Calculate the average Fst value over all loci
    return count/markersused
    print "radtags used =", markersused
    print "Fst =", count/markersused

def FstHohenTxt(radtagsfilename, haplotagsfilename, poplist, output):
	"""Calculate Fst with radtagsfilename, haplotagsfilename, and poplist and export restuls to file 'output'"""
    with open(output, "w+") as h:
        for i in range(len(poplist)):
            for j in range(i + 1, len(poplist)):
                h.write(poplist[i] + " " + poplist[j] + " " + str(FstRADHohen(radtagsfilename, haplotagsfilename, poplist[i], poplist[j])))
                h.write("\n")

if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(description="Calculate FST from molecular diversity")
	parser.add_argument('-i', '--inputhaplo')
	parser.add_argument('-h', '--inputdepths')
	parser.add_argument('-o', '--output')
	parser.add_argument('-p', '--poplist')
	args = parser.parse_args()
	FstHohenTxt(args.inputhaplo, args.inputdepths, args.poplist, args.output)     
