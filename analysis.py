import math
try:
    import queue
except ImportError:
    import Queue as queue
from igraph import *

def main():
	#inputfile = sys.argv[1]
	inputfile = "neuraminidase.txt"
	seq_array = readinputasarray(inputfile)
	seq_dict = {}
	for seq in seq_array:
		seq_dict[seq[0]] = seq[1] 
	#seq_array = seq_dict.items()
	#print seq_dict["CY044149"]
	#print seq_array[0]
	#print checklengths(seq_dict)
	if checklengths(seq_dict):
		distmatrix = createdistmatrix(seq_array)
		distgraph = Adjacency(distmatrix)
		mst = spanning_tree(distgraph, return_tree = True)

		#ultra_distmatrix = ultrametrify(distmatrix)
	else:
		print "sequence not all the same length"


def createdistmatrix(seq_array):
	distarray = []
	for seq1_num in range(0, len(seq_array)):
		seq_dist = [0 for x in range(0, len(seq_array))]
		#P is transition fraction
		#Q is transversion fraction
		for seq2_num in range(0, len(seq_array)):
			seq1 = seq_array[seq1_num][1]
			seq2 = seq_array[seq2_num][1]
			seq_dist[seq2_num] = tamuradist(seq1, seq2)
			#seq_dist[seq2_num] = 5
		distarray.append(seq_dist)
	#print distarray[0]
	return distarray

#jukes cantor distance of two sequences
def jcdist(seq1, seq2):
	substitutions = 0
	for base in range(0, len(seq1)):
		if seq1[base] != seq2[base]:
			substitutions += 1
	P = float(substitutions)/float(len(seq1))
	dist = -0.75 * math.log(1 - 4.0 / 3 * P)
	return dist

#kimura 2 parameter distance of two sequences
def kimuradist(seq1, seq2):
	transitions = 0
	transversions = 0
	for base in range(0, len(seq1)):
		if seq1[base] != seq2[base]:
			if (seq1[base] == "a" and seq2[base] == "g") or (seq1[base] == "g" and seq2[base] == "a") or (seq1[base] == "c" and seq2[base] == "t") or (seq1[base] == "t" and seq2[base] == "c"):
				transitions += 1
			else:
				transversions += 1
	P = float(transitions)/float(len(seq1))
	Q = float(transversions)/float(len(seq1))
	dist = 0
	if (1 - 2 * P - Q) == 0:
		dist = 0
	else:
		dist = -0.5 * math.log((1 - 2* P - Q) * math.sqrt(1 - 2 * Q))
	return dist

#tamura distances accounts for different nucleotide compositions too
def tamuradist(seq1, seq2):
	transitions = 0
	transversions = 0
	seq1gc = 0
	seq2gc = 0
	for base in range(0, len(seq1)):
		if seq1[base] == "g" or seq1[base] == "c":
			seq1gc += 1
		if seq2[base] == "g" or seq2[base] == "c":
			seq2gc += 1
		if seq1[base] != seq2[base]:
			if (seq1[base] == "a" and seq2[base] == "g") or (seq1[base] == "g" and seq2[base] == "a") or (seq1[base] == "c" and seq2[base] == "t") or (seq1[base] == "t" and seq2[base] == "c"):
				transitions += 1
			else:
				transversions += 1
	P = float(transitions)/float(len(seq1))
	Q = float(transversions)/float(len(seq1))
	theta1 = float(seq1gc) / len(seq1)
	theta2 = float(seq2gc) / len(seq2)
	C = theta1 + theta2 - 2 * theta1 * theta2 
	dist = 0
	if C == 0 or (1 - (P)/(C) - Q) == 0 or 1- 2* Q == 0:
		dist = 0
	else:
		dist = -C * math.log(1 - (P)/(C) - Q) - 0.5 * (1 - C) * math.log(1 - 2* Q) 
	return dist




def ultrametrify(distmatrix):
	#for
	pass 


	


def checklengths(seq_dict):
	allsamelength = True
	length = len(seq_dict["GQ243758"])
	for sequence in seq_dict:
		if len(seq_dict[sequence]) != length:
			allsamelength = False
	#print length
	return allsamelength


def readinputasarray(file):
	seq_array = []
	seq_name = ""
	current_seq = ""
	count = 0
	for line in open(file,'r'):
		if line[0] == ">":
			if count != 0:
				seq_array.append([seq_name, current_seq])
				current_seq = ""
				seq_name = line[2:len(line)-1]
			else:
				seq_name = line[2:len(line) - 1]
		else:
			line = line.replace(" ", "")
			line = line.replace("\n", "")
			current_seq += line
		count += 1
	seq_array.append([seq_name, current_seq])
	#print seq_dict
	return seq_array

main()
