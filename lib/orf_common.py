# -*- coding: utf-8 -*-

import re, os

assert os.getenv('BARSEQ_DIR') is not None, 'BARSEQ_DIR environment variable not defined. Please ensure you entered "export BARSEQ_DIR=/path/to/barseq_counter" before running this command.'

folder = os.path.join(os.getenv('BARSEQ_DIR'), 'data')


def getorf2common():
	f = open(os.path.join(folder, 'SGD_features.tab'))
	orf2common = {}
	for line in f:
		words = line.strip().split('\t')
		try:
			if words[3] != '' and words[4] != '':
				#print words[3], words[9], words[10]
				orf2common[words[3]] = words[4]
		except IndexError:
			pass
		except ValueError:
			pass
	f.close()
	return orf2common


def getcommon2orf():
	f = open(os.path.join(folder, 'SGD_features.tab'))
	common2orf = {}
	for line in f:
		words = line.strip().split('\t')
		try:
			if words[3] != '' and words[4] != '':
				#print words[3], words[9], words[10]
				common2orf[words[4]] = words[3]
		except IndexError:
			pass
		except ValueError:
			pass
	f.close()
	return common2orf
	
