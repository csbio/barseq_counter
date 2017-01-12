#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os

assert os.getenv('AGREP') is not None, 'AGREP environment variable not defined. Please ensure you entered "export AGREP=/path/to/agrep" before running this command.'

AGREP = os.getenv('AGREP')

def demultiplex(filename):
	f = open(filename)
	strings= []
	cmd = AGREP + ' -2 \'GATGTCCACGAGGTCTCT\' %s' % (filename)
	lines = os.popen(cmd).readlines()
	for line in lines[:-1]:
		string = line.strip()
		mtag = string[:10]
		barcode = string[28:48]
		key = mtag + '\t' + barcode
		strings.append(key)
	
	cmd = AGREP + ' -2 \'CGCTCCCGCCTTACTTCGCATTTAAA\' %s' % (filename) # for pombe
	lines = os.popen(cmd).readlines()
	for line in lines[:-1]:
		string = line.strip()
		mtag = string[:10]
		barcode = string[36:56]
		key = mtag + '\t' + barcode
		strings.append(key)
	
	cmd = AGREP + ' -2 \'AATCTTCGGTAGTCCAGCG\' %s' % (filename) # for ecoli
	lines = os.popen(cmd).readlines()
	for line in lines[:-1]:
		string = line.strip()
		mtag = string[:10]
		barcode = string[29:49]
		key = mtag + '\t' + barcode
		strings.append(key)
	
	return strings

def mergeFiles(filenames, outputfile):
	w = open(outputfile, 'w')
	for filename in filenames:
		w.write( "\n".join(demultiplex(filename)) + '\n')
	w.close()

def makeTxtFileForReport(filenames, outputfile):
	w = open(outputfile, 'w')
	for filename in filenames:
		f = open(filename)
		count = 0
		for line in f:
			count += 1
			if count % 4 == 2:
				w.write(line)
		f.close()
	w.close()

def main(barseqfolder, outputfile):
	filenames = ['%s/%s' % (barseqfolder, filename) for filename in os.listdir(barseqfolder) if filename.endswith('fastq')]
	mergeFiles(filenames, outputfile)
	makeTxtFileForReport(filenames, '%s/rawdataAndjustdata.txt' % barseqfolder)
	
import sys
if len(sys.argv) < 3:
	print "Usage: python preprocess_MIseq.py <folder containing barseq files (fastq extension)> <outputfile>"
	print "The barseq data folder should not contain any txt extension files. This program will generate a txt file which would be used by report generation program."
else:
	barseqfolder = sys.argv[1]
	outputfile = sys.argv[2]
	barseqfolder = barseqfolder[:-1] if barseqfolder[-1] == '/' else barseqfolder
	main(barseqfolder, outputfile)
