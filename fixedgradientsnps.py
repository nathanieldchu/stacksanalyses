#!/usr/bin/env python
"""
Created on Tue Aug 13 09:36:49 2013

@author: nathanieldchu
Takes a STACKS haplotypes export file and extracts stacks containing polymorphisms that follow a "gradient", meaning that at that locus, the identity of the polymorphism changes from one identity to another, without reverting to a previous identity. At each end of the gradient, the polymorphism identity must be "fixed", so that neither shares an allele with the other. This script is meant for DNA sequence data.

Usage:
./fixedgradientsnps.py -i <inputfile> -o <inputfile>

"""
def setwithorder(sequence):
	"""
	Make a set while retaining order.
	"""
   	seen = set()
	seen_add = seen.add
	return [x for x in sequence if x not in seen and not seen_add(x)]

def reorderkeepingorder(tuplist, numsnps):
	"""(tuple list) -> (tuple list)
	Takes a tuples list of alleles from populations after being split from STACKS and returns a list of lists, where each list is the sorted alleles for each SNP position, as the original order that they appeared
	reorder([['TGA', 'ACT'], ['TGA', 'ACG'], ['TGG', 'TCT']], 3)
		[[['T', 'A'], ['T', 'A'], ['T']], [['G', 'C'], ['G', 'C'], ['G', 'C']], [['A', 'T'], ['A', 'G'], ['G', 'T']]]
	"""
	alleles1 = []
	for i in range(numsnps):
		alleles2 = []
		for j in tuplist:
			alleles2.append([x[i] for x in j])
		alleles1.append([setwithorder(y) for y in alleles2])
	return alleles1
	
def ExtractLatitudinalPolymorphisms(radtagsfile, output):
	"""(file) -> (file)
	Take a radtags file and extract all stacks that vary latitudinally and are fixed at least at the terminal ends of the represented populations
	"""
	import csv
	with open(radtagsfile, "U") as f:
		radtags = csv.reader(f, delimiter="\t")
		radheader = radtags.next()
		snpindex = radheader.index("Num SNPs")
		catalogidindex = radheader.index("Catalog ID")
		consensusindex = radheader.index("Consensus Sequence")
		fixedcatalogIDs = []
		consensusseqs = []
		numreps = []
		for i in radtags:
			numrep = sum(1 for x in range(12, 18) if i[x] != "")
			if 0 < int(i[snpindex]) < 4 and numrep > 2:
				#all radtags that have less than 4 SNPs and more than 2 populations represented. Change the SNP cutoff if needed.
				allelespre = [i[y].split("/") for y in [16, 15, 13, 17, 14, 12] if i[y] != ""]
				#In line 47 the order of numbers represents the population "cline" that you are interested in. In this example, looking at latitudinal clines, 16 was the northernmost populations, and 12 was the southernmost. It could work in reverse as well. Change this to suite your needs.
				allelesord = reorderkeepingorder(allelespre, int(i[snpindex]))
				latvary = False
				fixedends = False
				for site in allelesord:
					if 'N' not in str(site):
						valdic = {}
						valuelist = []
						#Create a dictionary where each allele identity is assigned a value based on when that allele appeared
						for character in str(site):
							if character.isalpha() and character not in valdic:
								if valdic == {}:
									valdic[character] = 1
								else:
									valdic[character] = valdic[max(valdic)] + 1
						#Use the average of allele values to determine whether the populations are incrementally going from one identity to another
						for pop in site:
							avgvalue = sum(valdic[nucleotide] for nucleotide in pop) / float(len(pop))
							valuelist.append(avgvalue)
						if all(valuelist[x + 1] >= valuelist[x] for x in range(len(valuelist) - 1)):
							latvary = True
						if set(site[0]).intersection(set(site[-1])) == set([]):
							fixedends = True
				if latvary == True and fixedends == True:
					fixedcatalogIDs.append(i[catalogidindex])
					consensusseqs.append(i[consensusindex])
					numreps.append(numrep)
	with open(output, "w+") as g:
		for k in range(len(fixedcatalogIDs)):
			g.write(fixedcatalogIDs[k] + "\t" + consensusseqs[k] + "\t" + str(numreps[k]) + "\n")

if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(description="Extract fixed SNPs along a population gradient")
	parser.add_argument('-i', '--input')
	parser.add_argument('-o', '--output')
	args = parser.parse_args()
	ExtractLatitudinalPolymorphisms(args.input, args.output)
