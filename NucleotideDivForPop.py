#!/usr/bin/env python
"""
Created on Wed Apr 17 09:10:45 2013

@author: nathanieldavidchu
Calculates nucleotide diversity for all populations.

Usage:
./NucleotideDivForPop.py -i <inputhaplotypefile> -h <inputhaplodepthsfile> -p <poplist> -o <output>
"""

def NucleotideDiv(radtagsfilename, haplotagsfilename, pop1):
	""" (file, file, str) -> (float)
	Return nucleotide diversity for a given population, using only RAD tags that have fewer than 4 SNPs
	
	"""
	import csv
	import itertools
	count = 0
	markersused = 0
	with open(radtagsfilename, "U") as f:
		radtags = csv.reader(f, delimiter="\t")
		radheader = radtags.next()
		radindex1 = radheader.index(pop1)
		snpindex = radheader.index("Num SNPs")
		radpopA = []
		snps = []
		for i in radtags:
			if int(i[snpindex]) < 4 and i[radindex1] != "":
				radpopA.append(i[radindex1].split("/"))
				snps.append(int(i[snpindex]))
	with open(haplotagsfilename, "r+") as g:
		haplotags = csv.reader(g, delimiter = "\t")
		haploheader = haplotags.next()
		haploindex1 = haploheader.index(pop1)
		snpindex = radheader.index("Num SNPs")
		happopA = []
		for j in haplotags:
			if int(j[snpindex]) < 4 and j[haploindex1] != "":
				happopA.append(map(int, j[haploindex1].split("/")))
	for k in range(len(radpopA)):
		if len(radpopA[k]) == 1:
			markersused += 1
		else:
			combA = sorted(zip(radpopA[k], happopA[k]))
			radpopAset = set(radpopA[k])
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
			sumreadsA = sum([x[1] for x in combA])
			Nucdiv = 0
			for o in itertools.combinations(combA, 2):
				pi = sum(1 for x, y in zip(o[0][0], o[1][0]) if x != y) / 42
				Nucdiv += pi * (o[0][1]/sumreadsA) * (o[1][1]/sumreadsA)
			count += Nucdiv
			markersused += 1
	return count/markersused

	print "radtags used = ", markersused
	print "Nucleodiversity =", count/markersused
	
def NucleotideDivAllPop(radtagsfilename, haplotagsfilename, poplist, output):
	with open(output, "w+") as h:
		for i in poplist:
			h.write(i + " " + str(NucleotideDiv(radtagsfilename, haplotagsfilename, i)))
			h.write("\n")

if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(description="Calculate nucleotide diversity for all populations")
	parser.add_argument('-i', '--inputhaplo')
	parser.add_argument('-h', '--inputdepths')
	parser.add_argument('-p', '--poplist')
	parser.add_argument('-o', '--output')
	args = parser.parse_args()
	NucleotideDivAllPop(args.inputhaplo, args.inputdepths, args.poplist, args.output)
