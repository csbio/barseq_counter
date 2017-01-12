#!/usr/bin/env python
# -*- coding: utf-8 -*-
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')

import os, sys

assert os.getenv('BARSEQ_DIR') is not None, 'BARSEQ_DIR environment variable not defined. Please ensure you entered "export BARSEQ_DIR=/path/to/barseq_counter" before running this command.'

sys.path.append(os.path.join(os.getenv('BARSEQ_DIR'), 'lib'))

import scipy, matplotlib.pyplot as plt

def numberOfReads(folder):
	files = os.listdir(folder)
	totalcount = 0
	for filename in files:
		if filename.endswith('txt'): # check this condition before everyrun
			#print filename
			count = os.popen('wc %s/%s' % (folder, filename)).readlines()
			totalcount += int( count[0].strip().split()[0] )
	return totalcount

def numberOfLines(filename):
	count = os.popen('wc %s' % (filename)).readlines()
	return int( count[0].strip().split()[0] )

def numberOfCounts(filename):
	f = open(filename)
	count = 0
	for line in f:
		for word in line.strip().split('\t'):
			try:
				count += int(word)
			except ValueError:
				pass
	f.close()
	return count

def getData(filename):
	f = open(filename)
	mtags = f.readline().strip().split('\t')[1:]
	genes = []
	for line in f:
		words = line.strip().split('\t')
		genes.append(words[0])
	f.close()
	
	matrix = scipy.zeros((len(genes), len(mtags)))
	f = open(filename)
	mtags = f.readline().strip().split('\t')[1:]
	rowcount = 0
	for line in f:
		words = line.strip().split('\t')
		#print words[1:]
		for j, val in enumerate(words[1:]):
			matrix[rowcount, j] = float(val)
		rowcount += 1
	f.close()
	return genes, mtags, matrix
	

def writeCounts(items, counts, filename):
	Is = sorted(range(len(items)), key = lambda x:counts[x], reverse = True)
	f = open(filename, 'w')
	for i in Is:
		f.write('%s\t%s\n' % (items[i], counts[i]))
	f.close()

def makeHistogram(items, counts, filename):
	plt.figure()
	counts = sorted(counts, reverse = True)
	plt.bar(range(len(counts)), counts)
	if filename.find('mtag') != -1:
		plt.xlabel('Multiplex Tags (sorted)')
	else: 
		plt.xlabel('Barcodes (sorted)')
	plt.ylabel('Frequency')
	plt.savefig(filename)

def distribution(dataset, prefix):
	genes, drugs, matrix = dataset
	#print matrix
	#print drugs
	mtagcounts = scipy.sum(matrix, axis = 0)
	barcodecounts = scipy.sum(matrix, axis = 1)
	writeCounts(drugs, mtagcounts, '%s_mtagcounts.txt' % prefix)
	writeCounts(genes, barcodecounts, '%s_barcodecounts.txt' % prefix)
	makeHistogram(drugs, mtagcounts, '%s_mtagcountshist.png' % prefix)
	makeHistogram(genes, barcodecounts, '%s_barcodecountshist.png' % prefix)

def writeSummary(count1, count2, count3, dataset, prefix):
	f = open('%s_summary.txt' % prefix, 'w')
	f.write('Total read counts = %d\n' % count1)
	f.write('Readable read counts = %d (%.2f %%)\n' % (count2, float(count2)*100/count1) )
	f.write('Counted read counts = %d (%.2f %%)\n' % (count3, float(count3)*100/count1) )
	genes, drugs, matrix = dataset
	f.write('Size of the data = %d genes x %d conditions\n' % (len(genes), len(drugs)))
	f.close()

def main(seqfolder, barseq_raw_file, barseq_processed_file, prefix):
	count1 = numberOfReads(seqfolder)
	count2 = numberOfLines(barseq_raw_file)
	count3 = numberOfCounts(barseq_processed_file)
	dataset = getData(barseq_processed_file)
	#print count1, count2, '(%.2f %%)' % float(count2)/count1,  count3, '(%.2f %%)' % float(count3)/count1
	writeSummary(count1, count2, count3, dataset, prefix)
	distribution(dataset, prefix)


import sys

if len(sys.argv) < 5:
	print "Usage: python generateReport.py <seq folder> <barseq file> <barseq processed file> <output prefix>"
else:
	seqfolder = sys.argv[1]
	barseq_raw_file = sys.argv[2]
	barseq_processed_file = sys.argv[3]
	prefix = sys.argv[4]
	main(seqfolder, barseq_raw_file, barseq_processed_file, prefix)
	
