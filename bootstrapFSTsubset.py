#!/usr/bin/env python
"""
Created on Sun Sep  1 19:28:28 2013

@author: nathanieldavidchu
Calculate Fst using a binomial procedure as in Hohenlohe 2010. Then calculate the Fst value of a subset of markers. Using the full list of Fst values for each stack, then create a normal curve of mean Fst values when a subset of equal length to the subset of interest is bootstrap sampled from the total data set. Finally, the mean Fst of the subset in compared against this normal curve and the fraction of values that are above that subset value is returned. This script will use any number of cores to increase processing speed.

Usage:
./bootstrapFST.py -i <inputhaplotypefile> -h <inputalleledepthsfile> -o <outputfile> -d <datasubsetsize> -b <numberofbootstrapreplications> -c <listofcatalogidsforsubsetofinterest> -m <numberofcorestouse>

where -c is in a list format
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

def nucdivtotal(tuplelist, numsnps):
	"""Calculate pooled nucleotide diversity (binomial) for all populations at a given marker
	Take a list of tuples and export an int
	
	>>> nucdivtotal([('C/T', '6/4'), ('C/T', '3/8'), ('T', '3'), ('T', '4'), ('C/T', '4/3'), ('T', '10')], 1)
	0.010004810004810006
	>>> nucdivtotal([('CA/TC', '6/4'), ('CA/TA', '3/8')], 2)
	0.019954648526077097
	"""
	if numsnps == 1:
		alleledic1 = {}
		for pop in tuplelist:
			zipped = zip(pop[0].split('/'), pop[1].split('/'))
			for allele in zipped:
				if allele[0] in alleledic1:
					alleledic1[allele[0]] += int(allele[1])
				else:
					alleledic1[allele[0]] = int(allele[1])
		sumalleles = sum([alleleidepth for allelei, alleleidepth in alleledic1.items()])
		nucdiv = 1 - (sum([binomial(alleleidepth, 2) for allelei, alleleidepth in alleledic1.items()]) / binomial(sumalleles, 2))
		return nucdiv / 42
	elif numsnps == 2:
		alleledic1 = {}
		alleledic2 = {}
		for pop in tuplelist:
			zipped = zip(pop[0].split('/'), pop[1].split('/'))
			for allele in zipped:
				if allele[0][0] in alleledic1:
					alleledic1[allele[0][0]] += int(allele[1])
				else:
					alleledic1[allele[0][0]] = int(allele[1])
				if allele[0][1] in alleledic2:
					alleledic2[allele[0][1]] += int(allele[1])
				else:
					alleledic2[allele[0][1]] = int(allele[1])
		sumalleles = sum([alleleidepth for allelei, alleleidepth in alleledic1.items()])
		nucdiv1 = 1 - (sum([binomial(alleleidepth, 2) for allelei, alleleidepth in alleledic1.items()]) / binomial(sumalleles, 2))
		nucdiv2 = 1 - (sum([binomial(alleleidepth, 2) for allelei, alleleidepth in alleledic2.items()]) / binomial(sumalleles, 2))
		return (nucdiv1 + nucdiv2) / 42
	elif numsnps == 3:
		alleledic1 = {}
		alleledic2 = {}
		alleledic3 = {}
		for pop in tuplelist:
			zipped = zip(pop[0].split('/'), pop[1].split('/'))
			for allele in zipped:
				if allele[0][0] in alleledic1:
					alleledic1[allele[0][0]] += int(allele[1])
				else:
					alleledic1[allele[0][0]] = int(allele[1])
				if allele[0][1] in alleledic2:
					alleledic2[allele[0][1]] += int(allele[1])
				else:
					alleledic2[allele[0][1]] = int(allele[1])
				if allele[0][2] in alleledic3:
					alleledic3[allele[0][2]] += int(allele[1])
				else:
					alleledic3[allele[0][2]] = int(allele[1])
		sumalleles = sum([alleleidepth for allelei, alleleidepth in alleledic1.items()])
		nucdiv1 = 1 - (sum([binomial(alleleidepth, 2) for allelei, alleleidepth in alleledic1.items()]) / binomial(sumalleles, 2))
		nucdiv2 = 1 - (sum([binomial(alleleidepth, 2) for allelei, alleleidepth in alleledic2.items()]) / binomial(sumalleles, 2))
		nucdiv3 = 1 - (sum([binomial(alleleidepth, 2) for allelei, alleleidepth in alleledic3.items()]) / binomial(sumalleles, 2))
		return (nucdiv1 + nucdiv2 + nucdiv3) / 42

def nucdivind(tuplelist, numsnps):
	"""Calculate pooled nucleotide diversity (binomial) for each population individually
	Take a list of tuples and return a list of int tuples, each a nucleotide diversity for a population and the total of alleles/reads for that population
	
	>>> nucdivind([('C/T', '6/4'), ('C/T', '3/8'), ('T', '3'), ('T', '4'), ('C/T', '4/3'), ('T', '10')], 1)
	[(0.012698412698412698, 10), (0.01038961038961039, 11), (0.0, 3), (0.0, 4), (0.013605442176870748, 7), (0.0, 10)]
	>>> nucdivind([('CA/TC', '6/4'), ('CA/TA', '3/8')], 2)
	[(0.025396825396825397, 10), (0.01038961038961039, 11)]
	"""
	if numsnps == 1:
		nucdivlist = []
		for pop in tuplelist:
			alleledic1 = {}
			zipped = zip(pop[0].split('/'), pop[1].split('/'))
			for allele in zipped:
				if allele[0] in alleledic1:
					alleledic1[allele[0]] += int(allele[1])
				else:
					alleledic1[allele[0]] = int(allele[1])
			sumalleles = sum([alleleidepth for allelei, alleleidepth in alleledic1.items()])
			nucdiv = 1 - (sum([binomial(alleleidepth, 2) for allelei, alleleidepth in alleledic1.items()]) / binomial(sumalleles, 2))
			nucdivlist.append(((nucdiv / 42), sumalleles))
		return nucdivlist
	elif numsnps == 2:
		nucdivlist = []
		for pop in tuplelist:
			alleledic1 = {}
			alleledic2 = {}
			zipped = zip(pop[0].split('/'), pop[1].split('/'))
			for allele in zipped:
				if allele[0][0] in alleledic1:
					alleledic1[allele[0][0]] += int(allele[1])
				else:
					alleledic1[allele[0][0]] = int(allele[1])
				if allele[0][1] in alleledic2:
					alleledic2[allele[0][1]] += int(allele[1])
				else:
					alleledic2[allele[0][1]] = int(allele[1])
			sumalleles = sum([alleleidepth for allelei, alleleidepth in alleledic1.items()])
			nucdiv1 = 1 - (sum([binomial(alleleidepth, 2) for allelei, alleleidepth in alleledic1.items()]) / binomial(sumalleles, 2))
			nucdiv2 = 1 - (sum([binomial(alleleidepth, 2) for allelei, alleleidepth in alleledic2.items()]) / binomial(sumalleles, 2))
			nucdivlist.append((((nucdiv1 + nucdiv2) / 42), sumalleles))
		return nucdivlist
	elif numsnps == 3:
		nucdivlist = []
		for pop in tuplelist:
			alleledic1 = {}
			alleledic2 = {}
			alleledic3 = {}
			zipped = zip(pop[0].split('/'), pop[1].split('/'))
			for allele in zipped:
				if allele[0][0] in alleledic1:
					alleledic1[allele[0][0]] += int(allele[1])
				else:
					alleledic1[allele[0][0]] = int(allele[1])
				if allele[0][1] in alleledic2:
					alleledic2[allele[0][1]] += int(allele[1])
				else:
					alleledic2[allele[0][1]] = int(allele[1])
				if allele[0][2] in alleledic3:
					alleledic3[allele[0][2]] += int(allele[1])
				else:
					alleledic3[allele[0][2]] = int(allele[1])
			sumalleles = sum([alleleidepth for allelei, alleleidepth in alleledic1.items()])
			nucdiv1 = 1 - (sum([binomial(alleleidepth, 2) for allelei, alleleidepth in alleledic1.items()]) / binomial(sumalleles, 2))
			nucdiv2 = 1 - (sum([binomial(alleleidepth, 2) for allelei, alleleidepth in alleledic2.items()]) / binomial(sumalleles, 2))
			nucdiv3 = 1 - (sum([binomial(alleleidepth, 2) for allelei, alleleidepth in alleledic3.items()]) / binomial(sumalleles, 2))
			nucdivlist.append((((nucdiv1 + nucdiv2 + nucdiv3) / 42), sumalleles))
		return nucdivlist

def fstbycatalog(radtagsfilename, radtagdepthsfilename, cataloglistfile):
	""" (file, file, int, int) -> (float)
	Calculate binomial FST for a subset of the stacks file
	
	"""
	#make a list of catalogids for the subset
	catalogids = []
	with open(cataloglistfile, "r+") as f:
		for line in f:
			catalogids.append(line.split('\t')[0])
	import itertools
	#Calculate FST for a subset of stacks that match the list of catalogids just made
	fstaccumulator = 0
	numstacks = 0
	for radtag, raddepth in itertools.izip(open(radtagsfilename), open(radtagdepthsfilename)):
		radtagsplit = radtag.rstrip('\n').split('\t')
		raddepthsplit = raddepth.rstrip('\n').split('\t')
		numrep = sum([1 for x in [12, 13, 14, 15, 16, 17] if radtagsplit[x] != ""])
		if radtagsplit[0] == 'Catalog ID':
			pass
		elif radtagsplit[0] not in catalogids:
			#filter out stacks not in catalogids
			pass
		elif int(radtagsplit[7]) == 0 and numrep > 1:
			numstacks += 1
		elif 0 < int(radtagsplit[7]) < 4 and numrep > 1:
			#calculate FST for stack
			fstnum = 0
			fstdenom = 0
			zipped = zip([radtagsplit[x] for x in [12, 13, 14, 15, 16, 17] if radtagsplit[x] != ""], [raddepthsplit[x] for x in [12, 13, 14, 15, 16, 17] if raddepthsplit[x] != ""])
			totalnucdiv = nucdivtotal(zipped, int(radtagsplit[7]))
			for pop in nucdivind(zipped, int(radtagsplit[7])):
				fstnum += binomial(pop[1], 2) * pop[0]
				fstdenom += totalnucdiv * binomial(pop[1], 2)
			if fstdenom == 0:
				numstacks += 1
			else:
				fstaccumulator += 1 - (fstnum / fstdenom)
				numstacks += 1
	return fstaccumulator / numstacks
	
def listfst(radtagsfilename, radtagdepthsfilename):
	""" (file, file, int, int) -> (list)
	List Fst values for each catalogid
	"""
	import itertools
	#make a list of fst per stack for the full dataset
	fstlist = []			
	for radtag, raddepth in itertools.izip(open(radtagsfilename), open(radtagdepthsfilename)):
		radtagsplit = radtag.rstrip('\n').split('\t')
		raddepthsplit = raddepth.rstrip('\n').split('\t')
		numrep = sum([1 for x in [12, 13, 14, 15, 16, 17] if radtagsplit[x] != ""])
		if radtagsplit[0] == 'Catalog ID':
			pass
		elif int(radtagsplit[7]) == 0 and numrep > 1:
			fstlist.append(0)
		elif 0 < int(radtagsplit[7]) < 4 and numrep > 1:
			#calculate FST for stack
			fstnum = 0
			fstdenom = 0
			zipped = zip([radtagsplit[x] for x in [12, 13, 14, 15, 16, 17] if radtagsplit[x] != ""], [raddepthsplit[x] for x in [12, 13, 14, 15, 16, 17] if raddepthsplit[x] != ""])
			totalnucdiv = nucdivtotal(zipped, int(radtagsplit[7]))
			for pop in nucdivind(zipped, int(radtagsplit[7])):
				fstnum += binomial(pop[1], 2) * pop[0]
				fstdenom += totalnucdiv * binomial(pop[1], 2)
			if fstdenom == 0:
				fstlist.append(0)
			else:
				fst = 1 - (fstnum / fstdenom)
				fstlist.append(fst)
	return fstlist

def bootstrapfstworker(queue, fstlist, datasetsize, numberofreps, value):
	"""bootstrapfst worker function. Sums a number of random bootstrap samples from entries in fstlist. Number of samples is equal to numberofreps
	"""
	import random
	numberover = 0
	for k in range(numberofreps):
		fstvalue = float(sum(random.choice(fstlist) for l in range(datasetsize)) / datasetsize)
		if fstvalue >= value:
			numberover += 1
	queue.put(numberover)
	
def semiequalparts(x, num):
	"""take an int and divide it into num roughly equal parts as a list. This splits the replications for parallelization.
	"""
	equalparts = [(x // num) for i in range(num)]
	for i in range(x % num):
		equalparts[i] += 1
	return equalparts
	
def dump_queue(queue):
    """
    Empties all pending items in a queue and returns them in a list.
    """
    result = []
    for i in iter(queue.get, 'STOP'):
        result.append(i)
    return result
	
def bootstrapmultiprocess(radtagsfilename, radtagdepthsfilename, output, datasetsize, numberofbootstrapreps, cataloglistfile, cores):
	#Create bootstrap normal distribution of nucleotide diversity by randomly choosing tuples from the list and averaging them
	#take as many as the value of datasetsize
	#Will work in concert across multiple cores (denoted by the input variable 'cores')
	value = fstbycatalog(radtagsfilename, radtagdepthsfilename, cataloglistfile)
	import multiprocessing
	fstlist = listfst(radtagsfilename, radtagdepthsfilename)
	numberoverqueue = multiprocessing.Queue()
	bootstrapsplit = semiequalparts(numberofbootstrapreps, cores)
	for k in bootstrapsplit:
		bootp = multiprocessing.Process(target=bootstrapfstworker, args = (numberoverqueue, fstlist, datasetsize, k, value))
		bootp.start()
	bootp.join()
	numberoverqueue.put('STOP')
	numberover = sum(dump_queue(numberoverqueue))
	with open(output, "w+") as h:
		h.write(str(numberover / numberofbootstrapreps))
	return numberover / numberofbootstrapreps
	
if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(description="Bootstrap Fst Resampling")
	parser.add_argument('-i', '--inputhaplo')
	parser.add_argument('-h', '--inputdepths')
	parser.add_argument('-o', '--output')
	parser.add_argument('-d', '--datasubsetsize')
	parser.add_argument('-b', '--numberofbootstrapreplications')
	parser.add_argument('-c', '--listofcatalogidsforsubsetofinterest')
	parser.add_argument('-m', '--numberofcorestouse')
	args = parser.parse_args()
	bootstrapmultiprocess(args.inputhaplo, args.inputdepths, args.output, args.datasubsetsize, args.numberofbootstrapreplications, args.listofcatalogidsforsubsetofinterest, args.numberofcorestouse)
