#!/usr/bin/env python

"""
This is ClusterBOSS (Cluster Based On Sequence Similarity), a Python tool designed to cluster sequencing data from in vitro selection experiments into families of sequence similarity.
Input count files are assumed to be 'galaxy-type', that is, to have 3 head lines: number of unique sequences, total number of molecules and an empty line.

Author: Celia Blanco 
Email: celiablanco@ucla.edu

V1.0 2021-04-01 : Created	

Dependencies:
	• python-Levenshtein

How to run it:

$ python ClusterBOSS.py input_file d_cutoff n_min a_min c_min rec(y/n) keep_not_clustered(y/n)

e.g.:
$ python ClusterBOSS.py input_file 3 1 1 10 n n 

"""

import Levenshtein
import sys
import os
from operator import itemgetter, attrgetter

print("""
  ____ _           _            ____   ___  ____ ____  
 / ___| |_   _ ___| |_ ___ _ __| __ ) / _ \/ ___/ ___| 
| |   | | | | / __| __/ _ \ '__|  _ \| | | \___ \___ \ 
| |___| | |_| \__ \ ||  __/ |  | |_) | |_| |___) |__) |
 \____|_|\__,_|___/\__\___|_|  |____/ \___/|____/____/    
""")


def main():
	######################################################
	################### PARAMETERS #######################
	######################################################

	f_name_in= sys.argv[1]		   	# name of the input file

	dist_cutoff = int(sys.argv[2]) 	# d_cutoff: cutoff distance used to cluster 
	min_seqs = int(sys.argv[3])    	# n_min: minimum number of sequences per peak
	min_abd = int(sys.argv[4])     	# a_min: minimum abundance of sequences included in peaks - either 1 (includes singletons) or >1 
	min_abd_center = int(sys.argv[5])     	# c_min: minimum abundance of centers
	rec = str(sys.argv[6])			# y/n to whether sequences should be used in more than one peak
	keep_not_clustered = str(sys.argv[7])			# y/n to whether an output file with unclustered sequences should be generated

	print("Parameters:")
	print("")
	print("Input file:", f_name_in)
	print("Cutoff distance:", dist_cutoff)
	print("Minimum number of sequences per peak:", min_seqs)
	print("Minimum abundance of sequences in peaks:", min_abd)
	print("Minimum abundance of center sequences:", min_abd_center)
	if rec =="y":
		print("Recycle: required")
	if rec =="n":
		print("Recycle: not required")
	if keep_not_clustered =="y":
		print("Generate file with unclustered sequences: required")
	if keep_not_clustered =="n":
		print("Generate file with unclustered sequences: not required")
	print("")
	print("------------------------------------------------------")
	print("")
	if not os.path.exists("e"+str(dist_cutoff)):
		os.makedirs("e"+str(dist_cutoff))

	######################################################
	################## CREATE OUTPUT FILES ###############
	######################################################

	if rec == 'y':
		rec_num = 2
		f_name_out= "e"+str(dist_cutoff)+ "/" + f_name_in.split(".")[0] + "_e"+str(dist_cutoff)+"_ms"+str(min_seqs)+"_ma"+str(min_abd) + "_mac"+str(min_abd_center)+ "_rec_peaks.txt"    # name of the output file
		out=open(f_name_out, 'w')
		if keep_not_clustered == 'y':
			f_name_out2= "e"+str(dist_cutoff)+ "/" + f_name_in.split(".")[0] + "_e"+str(dist_cutoff)+"_ms"+str(min_seqs)+"_ma"+str(min_abd) + "_mac"+str(min_abd_center) + "_rec_nopeaks.txt"    # name of the output file
			out2=open(f_name_out2, 'w')
	
	if rec == 'n':
		rec_num = 1
		f_name_out= "e"+str(dist_cutoff)+ "/" + f_name_in.split(".")[0] + "_e"+str(dist_cutoff)+"_ms"+str(min_seqs)+"_ma"+str(min_abd) + "_mac"+str(min_abd_center) + "_norec_peaks.txt"    # name of the output file
		out=open(f_name_out, 'w')
		if keep_not_clustered == 'y':
			f_name_out2= "e"+str(dist_cutoff)+ "/" + f_name_in.split(".")[0] + "_e"+str(dist_cutoff)+"_ms"+str(min_seqs)+"_ma"+str(min_abd) + "_mac"+str(min_abd_center) + "_norec_nopeaks.txt"    # name of the output file
			out2=open(f_name_out2, 'w')
	
	#######################################################################################################
	################################## READ INPUT (COUNTS) FILE ###########################################
	#######################################################################################################

	print("Reading input counts file ...")

	tot = 0
	max_len=0
	with open(f_name_in) as readin:
		sequences = []
		next(readin) #these 3 lines are to skip the header (the count files in galaxy format have 3 lines with total counts before the first line with a seq)
		next(readin) #these 3 lines
		next(readin) #these 3 lines
		for line in readin:
			linesp = line.split()
			if len(linesp) == 2:
				sequences.append({'seq': linesp[0], 'abd': int(linesp[1]), 'dist': 0, 'status': 1, 'len': len(linesp[0])})    #sequence, abundance, distance to center (for the future), status (for the future, will change to 2 once a seq has been used in a previous cluster)
				len_seq=len(linesp[0]) #these 3 lines are to calculate the max lenght in the sequences
				tot += int(linesp[1])
				if len_seq > max_len:  #these 3 lines
					max_len = len_seq  #these 3 lines
	readin.close()

	sorted_seq=[]
	sorted_seq=sorted(sequences, reverse=True, key=lambda x: x['abd'])

	################################################
	#################### CLUSTER ###################
	################################################

	print("Clustering sequences ...")

	peak_number = 1     # peak_number will be the number of peaks once the script finishes
	seq_number=1        # seq_number will be the number of seqs per peak (it will change for every peak and that's why we reset it to 1 when we are done with each cluster)

	out.write(str('peak_rank') + " / " +str('sequence_rank') + " / " + str('sequence') + " / " +  str('abundance') + " / "  + str('frequency') + " / " + str('distance to center') + " / " + str('status') + '\n')      # print the center seq in the file

	if keep_not_clustered == 'y':
		out2.write(str('sequence') + " / " + str('abundance') + " / " + str('frequency') + " / " + str('distance to center') + " / " + str('status') + '\n')      # print the center seq in the file

	for j in range (0,len(sorted_seq)):
		n_seq = 1
		printable=[]
		cands = []
		if sorted_seq[j]['abd'] >= min_abd_center:
			if sorted_seq[j]['status'] == 1 and sorted_seq[j]['len'] >= dist_cutoff :   # only sequences with st=1 and longer than the cutoff distance can be centers
				center = sorted_seq[j]     # define the center seq
				fr_center = float(center['abd'])/float(tot)
				center['freq'] = '%.4f' % fr_center
				sorted_seq.pop(j)     # remove the center from the list of seqs (so it cant be further used)
				sorted_seq.insert(j,{'seq': 'xxx', 'abd': 0, 'dist': 0, 'status': 0, 'freq': 0})      # this line is to ensure that the for loop doesnt skip a line after popping out the center
				for i, cand in enumerate(sorted_seq):
					if center['abd'] >= cand['abd'] >= min_abd:    # both the center and the candidate need to have abd > min_abd    
						if cand['len'] >= dist_cutoff:
							if cand['status'] == 1 or cand['status'] == rec_num:    # sequences can be clustered in more than one peak or not (depending of y/n asked for in the command line)
								cand['dist'] = Levenshtein.distance(center['seq'], cand['seq'])    # define distance between the center and the candidate
								if cand['dist'] <= dist_cutoff: 
									n_seq += 1
									fr_cand = float(cand['abd'])/float(tot)
									cand['freq'] = '%.4f' % fr_cand
									cands.append({'rank': str(n_seq), 'seq': str(cand['seq']), 'abd': str(cand['abd']), 'freq': str(cand['freq']), 'dist': str(cand['dist']), 'st': str(cand['status']) })	# store the candidate seqs to print them at the end (this is to make sure that there are at least min_seqs in the peak)
									sorted_seq[i]['status'] = 2
									sorted_seq.pop(i)     # remove the center from the list of seqs (so it cant be further used)
									sorted_seq.insert(i,{'seq': 'xxx', 'abd': 0, 'dist': 0, 'status': 0, 'freq': 0})      # this line is to ensure that the for loop doesnt skip a line after popping out the center
				if n_seq >= min_seqs:     # making sure that there are at least min_seqs in the peak
					out.write(str("---------- peak = ") + str(peak_number) + str("----------") + '\n')
					out.write('\n')
					out.write(str(peak_number).ljust(5) + str(seq_number).ljust(10) + str(center['seq']).ljust(len(center['seq'])+dist_cutoff +10)  +  str(center['abd']).ljust(20) + str(center['freq']).ljust(20) + str(0).ljust(20) + 'st=' + str(center['status']) + '\n')      # print the center seq in the file
					seq_number += 1
					for i in range (0, len(cands)):
						out.write(str(peak_number).ljust(5) + str(seq_number).ljust(10) + str(cands[i]['seq']).ljust(len(center['seq'])+dist_cutoff + 10) + str(cands[i]['abd']).ljust(20) + str(cands[i]['freq']).ljust(20) + str(cands[i]['dist']).ljust(20) + "st=" + str(cands[i]['st']) + '\n')     # print in the file the seqs stored in cands
						seq_number += 1     # add 1 to the number of seqs in the cluster
					out.write('\n')
					peak_number+=1    # count total peaks
					seq_number=1      # reset the seq number (so it is 1 when the next peak starts clustering)
				else:
					if keep_not_clustered == 'y':
						out2.write(str(center['seq']).ljust(len(center['seq'])+dist_cutoff +10)  +  str(center['abd']).ljust(20) + str(center['freq']).ljust(20) + str(0).ljust(20) + 'st=' + str(center['status']) + '\n')      # print the center seq in the file
						for i in range (0, len(cands)):
							out2.write(str(peak_number).ljust(5) + str(seq_number).ljust(10) + str(cands[i]['seq']).ljust(len(center['seq'])+dist_cutoff + 10) + str(cands[i]['abd']).ljust(20) + str(cands[i]['freq']).ljust(20) + str(cands[i]['dist']).ljust(20) + "st=" + str(cands[i]['st']) + '\n')	# print in the file the seqs stored in cands
							seq_number += 1 
	out.close()

	if keep_not_clustered == 'y':
		for element in sorted_seq:
			if element['seq'] != 'xxx':
				if element['status'] == 1:
					fr_element = float(element['abd'])/float(tot)
					element['freq'] = '%.4f' % fr_element
					out2.write(str(element['seq']) + str("\t") + str(element['abd'])+ str("\t") + str(element['freq'])+ str("\t") + str("-") + str("\t") + str("-") + '\n')
		out2.close()
	
	print("")
	print("Finished clustering sequences.")
	print("There are ", len(sequences), " different sequences in ",  f_name_in)
	print("Clustering resulted in ", peak_number-1, " peaks.")
	print("")
if __name__=='__main__':
	main()
	