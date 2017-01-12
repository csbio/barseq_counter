#!/usr/bin/env python
# -*- coding: utf-8 -*-
import scipy, os

assert os.getenv('AGREP') is not None, 'AGREP environment variable not defined. Please ensure you entered "export AGREP=/path/to/agrep" before running this command.'

#AGREP = '/project/csbio/barseq_counter/agrep-3.41/agrep'
AGREP = os.getenv('AGREP')

def nmismatches(string1, string2):
	count = 0
	for i in range(len(string1)):
		if string1[i] != string2[i]:
			count += 1
	return count

def writeList(lst, filename):
	f = open(filename, 'w')
	for e in lst:
		f.write('%s\n' % (str( e ))) 
	f.close()

#def getData(barseqfile, multiplextags):
	#mtagbarcode2count = {}
	#f= open(barseqfile)
	#for line in f: 
		#string = line.strip()
		#if nmismatches( string[8:26], 'GATGTCCACGAGGTCTCT' ) > 2:
			#continue
		#mtag = string[:8]
		#if mtag not in multiplextags:
			#continue
		#barcode = string[26:46]
		##observedbarcodes.append(barcode)
		#key = mtag + '\t' + barcode
		#try: 
			#mtagbarcode2count[key] += 1
		#except KeyError:
			#mtagbarcode2count[key] = 1
	#return mtagbarcode2count

def getData(barseqfile, multiplextags):
	mtagbarcode2count = {}
	f= open(barseqfile)
	count = 0
	for line in f:
		count += 1
		#print count, line.strip()
		[mtag, barcode] = line.strip().split()
		if mtag not in multiplextags:
			continue
		key = mtag + '\t' + barcode
		try: 
			mtagbarcode2count[key] += 1
		except KeyError:
			mtagbarcode2count[key] = 1
	return mtagbarcode2count

# 1 mismatch for mtags is illegal operation because there are conflicts even in cases there are perfect matches. Examples:
#TCTCATAT ['TCTCATAT', 'TCTTCATA']
#GTATCGCA ['GTATGCGC', 'GTATCGCA']
#CACCGTCC ['CACCGTCC', 'CTACCGTC']
#AACTTAAT ['ACTTCAAT', 'AACTTAAT']
#GTGTTGAA ['TGTTGGAA', 'GTGTTGAA']
#TGGTGGGT ['TGGGTGGG', 'TGGTGGGT']
# AGAIN NOT USING THE FOLLOWING FUNCTION - ILLEGAL TO USE 
def merge1MismatchMtag(mtagbarcode2count, multiplextags):
	mtagsindata2count = {}
	for mtagbarcode, count in mtagbarcode2count.iteritems():
		mtag = mtagbarcode.split()[0]
		try: 
			mtagsindata2count[mtag] += count
		except KeyError:
			mtagsindata2count[mtag] = count
	mtagsindata = set( mtagsindata2count.keys() ) - set(multiplextags)
	writeList(mtagsindata, 'MtagsInData.txt')
	mtag2original = {}
	for mtag in multiplextags:
		cmd = '%s -w -1 \'%s\' %s' % (AGREP, mtag, 'MtagsInData.txt' )
		lines = os.popen(cmd).readlines()
		for line in lines[:-1]:
			modified = line.strip()
			try:
				mtag2original[modified] += [mtag]
			except KeyError:
				mtag2original[modified] = [mtag]
	writeList(multiplextags, 'MtagsExpected.txt')
	modifiedmtags = sorted(list( set(mtagsindata) - set(mtag2original) ), key = lambda mtag:mtagsindata2count[mtag], reverse = True )
	print '#tags in data, #tags mapped, #not mapped (stage 1)= %d, %d, %d' % (  len(mtagsindata), len(mtag2original), len(modifiedmtags) )
	print '# tags proper match %d' % (len([modified for modified, originals in mtag2original.iteritems() if len(originals) == 1]))
	#for modified in modifiedmtags:
		##print mtagsindata2count[modified], modified
		#cmd = '/project/csbio/raamesh/projects/laptopprojects/chemicalGenomics/barcodeCounting/Pipeline/tools/agrep-3.41/agrep -1 \'%s\' %s' % (modified, 'MtagsExpected.txt' )
		#lines = os.popen(cmd).readlines()
		#if len(lines) == 2:
			#mtag2original[modified] = [lines[0].strip()]
		##if len(lines) > 1:
			##continue
	mtagstobepopped = []
	for modified, originals in mtag2original.iteritems():
		if len(originals) != 1:
			mtagstobepopped.append(modified)
			#print modified, originals
	for modified in mtagstobepopped:
		mtag2original.pop(modified)
	
	for mtag in multiplextags:
		mtag2original[mtag] = [mtag]
	print "Mapping dictionary has %d mappings" % len(mtag2original)
	newmtagbarcode2count = {}
	count1, count2, count3 = 0, 0, 0
	for mtagbarcode, count in mtagbarcode2count.iteritems():
		mtag, barcode = mtagbarcode.split()
		try:
			newmtag = mtag2original[mtag][0]
			if newmtag != mtag: count1+=1
			else: count2+=1
		except KeyError:
			count3+=1
			continue
		try: 
			newmtagbarcode2count[newmtag+'\t'+barcode] += count
		except KeyError:
			newmtagbarcode2count[newmtag+'\t'+barcode] = count
	#print count1, count2, count3
	#print sum(mtagbarcode2count.values()), sum(newmtagbarcode2count.values())
	return newmtagbarcode2count

def getbarcode2original(barcodes):
	barcode2original = {}
	for barcode in barcodes:
		barcode2original[barcode] = [barcode]
		for i in range(len(barcode)):
			for char in set( ['A', 'T', 'G', 'C', 'N'] ) - set([barcode[i]]):
				newbarcode = barcode[:i] + char + barcode[i+1:]
				if newbarcode not in barcode2original:
					barcode2original[newbarcode] = [barcode]
				else:
					barcode2original[newbarcode] = barcode2original[newbarcode] + [barcode]
				for j in range(i):
					for char in set( ['A', 'T', 'G', 'C', 'N'] ) - set([barcode[j]]):
						newnewbarcode = newbarcode[:j] + char + newbarcode[j+1:]
						if newnewbarcode not in barcode2original:
							barcode2original[newnewbarcode] = [barcode]
						else:
							barcode2original[newnewbarcode] = barcode2original[newnewbarcode] + [barcode]
	return barcode2original

def getbarcode2original1(barcodes, mtagbarcode2count):
	allbarcodes_dict = getAllBarcodesInTheData(mtagbarcode2count)
	allbarcodes = list(set(allbarcodes_dict))
	f = open('barcodesInData.txt', 'w')
	for barcode in allbarcodes:
		f.write('%s\n' % barcode)
	f.close()
	
	f = open('barcodesExpected.txt', 'w')
	for barcode in barcodes:
		f.write('%s\n' % barcode)
	f.close()
	
	barcode2original = {}
	for i, original in enumerate(barcodes):
		print i
		cmd = '%s -w -2 \'%s\' %s' % (AGREP, original, 'barcodesInData.txt' )
		lines = os.popen(cmd).readlines()
		for line in lines[:-1]:
			barcode = line.strip()
			if barcode not in barcode2original:
				barcode2original[barcode] = [original]
			else:
				barcode2original[barcode] += [original]
	barcodestobepopped = []
	for barcode, originals in barcode2original.iteritems():
		if len(originals) != 1:
			barcodestobepopped.append(barcode)
	for barcode in barcodestobepopped:
		barcode2original.pop(barcode)
	print len(sorted([(barcode, allbarcodes_dict[barcode]) for barcode in set(allbarcodes) - set(barcode2original) if allbarcodes_dict[barcode] > 100], key = lambda x:x[1], reverse = True))
	for i, (barcode, count) in  enumerate( sorted([(barcode, allbarcodes_dict[barcode]) for barcode in set(allbarcodes) - set(barcode2original) if allbarcodes_dict[barcode] > 100], key = lambda x:x[1], reverse = True) ):
		print i
		for n in range(1, 4):
			cmd = '%s -w -%d \'%s\' %s' % (AGREP, n, barcode, 'barcodesExpected.txt' )
			lines = os.popen(cmd).readlines()
			if len(lines) == 2:
				#print barcode, count, lines[0].strip(), n
				barcode2original[barcode] = [lines[0].strip()]
			if len(lines) > 1:
				break
	return barcode2original


def getAllOriginalBarcodesInTheData(mtagbarcode2count, barcode2original):
	allbarcodes= [] # in this data
	for tagbarcode, count in mtagbarcode2count.iteritems():
		[tag, barcode] = tagbarcode.split()
		try:
			if len(set(barcode2original[barcode])) > 1:
				print "Warning: barcode conflict !!", barcode2original[barcode], barcode
			originalbarcode = barcode2original[barcode][0]
			allbarcodes.append(originalbarcode)
		except KeyError:
			continue
	allbarcodes = list(set(allbarcodes))
	return allbarcodes
	
def getAllBarcodesInTheData(mtagbarcode2count):
	allbarcodes_dict= {} # in this data
	for tagbarcode, count in mtagbarcode2count.iteritems():
		[tag, barcode] = tagbarcode.split()
		try:
			allbarcodes_dict[barcode] += count
		except KeyError:
			allbarcodes_dict[barcode] = count
	return allbarcodes_dict

def getMapping(lst):
	return dict([(item, i) for i, item in enumerate(lst)])

def getMatrix(allbarcodes, alltags, mtagbarcode2count, barcode2original):
	tagmapping = getMapping(alltags)
	barcodemapping = getMapping(allbarcodes)
	matrix = scipy.zeros((len(allbarcodes), len(alltags)))
	for tagbarcode, count in mtagbarcode2count.iteritems():
		[tag, barcode] = tagbarcode.split()
		try:
			if len(set(barcode2original[barcode])) > 1:
				print "Warning: barcode conflict !!", barcode2original[barcode], barcode
			originalbarcode = barcode2original[barcode][0]
			matrix[ barcodemapping[originalbarcode], tagmapping[tag] ] += int(count)
		except KeyError:
			continue
	return matrix

def makeMatrix(mtagbarcode2count, multiplextag2condition, barcode2original, barcode2gene, outputfile):
	alltags = list(multiplextag2condition)
	allbarcodes = getAllOriginalBarcodesInTheData(mtagbarcode2count, barcode2original)
	matrix = getMatrix(allbarcodes, alltags, mtagbarcode2count, barcode2original)
	
	print scipy.sum(matrix)
	f = open(outputfile, 'w')
	f.write('%s\t' % ('Genes \ Multiplex tags'))
	for i, tag in enumerate(alltags):
		if i == (len(alltags) - 1):
			f.write('%s_%s' % (multiplextag2condition[tag], tag))
		else:
			f.write('%s_%s\t' % (multiplextag2condition[tag], tag))  
	f.write('\n')
	for i, barcode in enumerate(allbarcodes):
		gene = '%s_%s' % (barcode2gene[barcode], barcode)
		#gene = barcode2gene[barcode]
		f.write('%s\t' % gene)
		for j in range(len(alltags)):
			if j == range(len(alltags))[-1]:
				f.write('%d' % matrix[i, j])
			else:
				f.write('%d\t' % matrix[i, j])
		f.write('\n')
	f.close()


def loadDict(filename):
	dictionary = {}
	f = open(filename)
	for line in f:
		words = line.strip().split('\t')
		if len(words) == 1:
			dictionary[words[0]] = ''
		else:
			dictionary[words[0]] = words[1]
	f.close()
	return dictionary

def main(barseqfile, multiplextagfile, barcodefile, outputfile):
	multiplextag2condition = loadDict(multiplextagfile)
	barcode2gene = loadDict(barcodefile)
	mtagbarcode2count = getData(barseqfile, list(multiplextag2condition))
	#mtagbarcode2count = merge1MismatchMtag(mtagbarcode2count, list(multiplextag2condition)) # illegal
	barcode2original = getbarcode2original1(list(barcode2gene), mtagbarcode2count)
	makeMatrix(mtagbarcode2count, multiplextag2condition, barcode2original, barcode2gene, outputfile)

import sys
if len(sys.argv) < 4:
	print "Usage: python processBarSeq_rd.py <barseqfile> <multiplex tag file> <barcode file> <outputfile>"
else:
	barseqfile = sys.argv[1]
	multiplextagfile = sys.argv[2]
	barcodefile = sys.argv[3]
	outputfile = sys.argv[4]
	main(barseqfile, multiplextagfile, barcodefile, outputfile)
