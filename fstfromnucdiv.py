#!/usr/bin/env python
"""
Created on Wed Apr 17 13:59:07 2013

@author: nathanieldchu
Take a STACKS haplotype export and allele depth exports calculates pairwise FST values for each pair of populations in poplist from molecular diversity

Usage:
./fstfromnucdiv.py -i <inputhaplotypefile> -h <inputalleledepthsfile> -o <outputfile> -p <poplist>

where <poplist> is a string of the populations to be compared
"pop1 pop2 pop3 ..."
"""
def NucDiv(x):
	"""Takes a list of tuples and calculates NucDiv based on pairwise differences and frequencies
	"""
	import itertools
	Nucdiv = 0
	sumreads = sum([z[1] for z in x])
	for i in itertools.combinations(x, 2):
		pi = sum(1 for m, n in zip(i[0][0], i[1][0]) if m != n) / 42
		Nucdiv += pi * (i[0][1]/sumreads) * (i[1][1]/sumreads)
	return float(Nucdiv)

def mergetuplelist(x, y):
	"""Return a merged tuple list retaining unique values and summing read counts"""
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

def FstFromNucDiv(radtagsfilename, haplotagsfilename, pop1, pop2):
	"""Calculate Fst based on NucDiv, read counts, and pairwise differences"""
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
			if i[radindex1] != "" and i[radindex2] != "":
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
			if j[haploindex1] != "" and j[haploindex2] != "":
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
			sumreadcountA = sum([x[1] for x in combA])
			sumreadcountB = sum([x[1] for x in combB])
			sumreadcount = sum([x[1] for x in merged])
			numerator = (sumreadcountA*NucDiv(combA) + sumreadcountB*NucDiv(combB))
			denomenator = sumreadcount*NucDiv(merged)
			#Add up Fst values per loci and keep track of how many have been analysed (markersused)
			if denomenator == 0:
				markersused += 1
			else:
				count += 1 - numerator/denomenator
				markersused += 1
	#Calculate the average Fst value over all loci
	return float(count/markersused)
	print "radtags used =", markersused
	print "Fst =", count/markersused
	
def FSTfromNucDivAllpop(radtagsfilename, haplotagsfilename, poplist, output):
	"""Calculate Fst with radtagsfilename, haplotagsfilename, and poplist and export restuls to file 'output'"""
	with open(output, "w+") as h:
		for i in range(len(poplist.split(" "))):
			for j in range(i + 1, len(poplist)):
				h.write(poplist[i] + " " + poplist[j] + " " + str(FstFromNucDiv(radtagsfilename, haplotagsfilename, poplist[i], poplist[j])))
				h.write("\n")

if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(description="Calculate FST from molecular diversity")
	parser.add_argument('-i', '--inputhaplo')
	parser.add_argument('-h', '--inputdepths')
	parser.add_argument('-o', '--output')
	parser.add_argument('-p', '--poplist')
	args = parser.parse_args()
	FSTfromNucDivAllpop(args.inputhaplo, args.inputdepths, args.poplist, args.output)
