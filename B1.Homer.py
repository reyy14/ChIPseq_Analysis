#!/usr/bin/python
import sys
from string import *
import os
import os.path
sys.path.append("/home/won/WON/")
import Util
from multiprocessing import Pool

FILES = ["/home/skw/Desktop/REy/Ben_Rey/SRR1130791.bowtie","/home/skw/Desktop/REy/Ben_Rey/H3K9ac.bam","/home/skw/Desktop/REy/Ben_Rey/H3K9me3.bam"]

def Job(files):
	odirec = files.replace("bowtie","homer")
	odirec = odirec.replace("bam","homer")
	cmd = 'makeTagDirectory %s %s -tbp 1 '%(odirec, files)
	print cmd
	#sys.exit(1)
	os.system(cmd)
	return

if __name__=='__main__':
	count = 0
	pool = Pool(processes = 3)
	pool.map(Job,FILES)
