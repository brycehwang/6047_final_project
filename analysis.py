import math

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
		distmatrix = kimuradistmatrix(seq_array)
		#ultra_distmatrix = ultrametrify(distmatrix)
	else:
		print "sequence not all the same length"

def kimuradistmatrix(seq_array):
	distarray = []
	for seq1_num in range(0, len(seq_array)):
		seq_dist = [0 for x in range(0, len(seq_array))]
		#P is transition fraction
		#Q is transversion fraction
		for seq2_num in range(0, len(seq_array)):
			seq1 = seq_array[seq1_num][1]
			seq2 = seq_array[seq2_num][1]
			transitions = 0
			transversions = 0
			for base in range(0, len(seq1)):
				if seq1[base] != seq2[base]:
					if (seq1[base] == "A" and seq2[base] == "G") or (seq1[base] == "G" and seq2[base] == "A") or (seq1[base] == "C" and seq2[base] == "T") or (seq1[base] == "T" and seq2[base] == "C"):
						transitions += 1
					else:
						transversions += 1
			P = float(transitions)/float(len(seq1))
			Q = float(transversions)/float(len(seq1))
			if (1 - 2 * P - Q) == 0:
				seq_dist[seq2_num] = 0
			else:
				dist = -0.5 * math.log((1 - 2* P - Q) * math.sqrt(1 - 2 * Q))
				print dist
				seq_dist[seq2_num] = dist
			#seq_dist[seq2_num] = 5
		distarray.append(seq_dist)
	print distarray[0][0]
	return distarray

		

def ultrametrify(distmatrix):
	mst = prims(distmatrix)

def prims(distmatrix):
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
