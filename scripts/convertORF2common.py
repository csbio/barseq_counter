#!/usr/bin/env python
# replaces all occurrences of ORFs to common names
import os, sys

assert os.getenv('BARSEQ_DIR') is not None, 'BARSEQ_DIR environment variable not defined. Please ensure you entered "export BARSEQ_DIR=/path/to/barseq_counter" before running this command.'

sys.path.append(os.path.join(os.getenv('BARSEQ_DIR'), 'lib'))

import orf_common as orthology, re


def getPattern(orfs):
	#pattern = re.compile('Y\w{2}\d{3}\w')
	#otherorfs = ['Y\w{2}\d{3}\w']
	#for orf in orfs:
		#if re.match(pattern, orf) == None:
			#otherorfs.append(orf)
	#return '|'.join(otherorfs)
	return 'Y\w{2}\d{3}\w(-\w){0,1}'


def convertORF2common(inputfilename, outputfilename):
	text = open(inputfilename).read()
	orf2common = orthology.getorf2common()
	def repl(matchobj):
		try: 
			return orf2common[matchobj.group(0)]
		except KeyError:
			return matchobj.group(0)
	#pattern = '|'.join(orf2common.keys()) # OverflowError: regular expression code size limit exceeded
	pattern = getPattern(orf2common.keys())
	text = re.sub(pattern, repl, text)
	# too slow
	#for orf, common in orf2common.iteritems():
		#print orf
		#text = re.sub(orf, common, text)
	#
	
	open(outputfilename, 'w').write(text)
	
import sys
if len(sys.argv) < 3:
	print "Usage: python convertORF2common.py <inputfile> <outputfile>"
	print "This script will replaces all instances of ORF with corresponding common names"
else:
	inputfilename = sys.argv[1]
	outputfilename = sys.argv[2]
	convertORF2common(inputfilename, outputfilename)
