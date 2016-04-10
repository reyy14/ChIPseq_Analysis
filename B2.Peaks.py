#!/usr/bin/python
import sys
from string import *
import os
import os.path
sys.path.append("/home/won/WON/")
import Util
from multiprocessing import Pool

HomerDirec = ["/home/skw/Desktop/REy/Ben_Rey/H3K9ac.homer","/home/skw/Desktop/REy/Ben_Rey/H3K9me3.homer","SRR1130791.homer"]

def Job(direc):
	cmd = 'findPeaks %s -style factor -o auto '%direc
	print cmd
	os.system(cmd)
	return

if __name__=='__main__':
	count = 0
	pool = Pool(processes = 3)
	pool.map(Job,HomerDirec)
