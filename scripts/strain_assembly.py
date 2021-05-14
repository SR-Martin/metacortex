#!/usr/bin/env python

import sys, getopt, errno, array

try:
	opts, args = getopt.getopt(sys.argv[1:],"hg:b:")
except getopt.GetoptError:
	print("Option not recognised.")
	print("python strain_assembly.py -g <GFA2 file> -b <blast file>")
	print("python strain_assembly.py -h for further usage instructions.")
	sys.exit(2)
for opt, arg in opts:
	if opt == "-h":
		print("python strain_assembly.py -g <GFA2 file> -b <blast file> ")
		print("Assembles contigs for each reference from GFA2 file by choosing best path using blast mapping.")
		print("BLAST input must have been created with the option -outfmt '6 qseqid sseqid pident length qlen'")
		print("-g <GFA2 file>\t\t Name of GFA2 file output by MetaCortex")
		print("-b <blast file> Blast file containing hits to all sequences in GFA2 file against reference.")
		sys.exit()
	elif opt in ("-g"):
		gfa2Filename = arg
	elif opt in ("-b"):
		blastFilename = arg

references = []
sequences = dict()

class Alignment:
	def __init__(self, query, ref, queryCoverage, identity):
		self.query = query
		self.ref = ref
		self.score = queryCoverage * identity

class Sequence:
	def __init__(self, name, sequence):
		self.name = name
		self.sequence = sequence
		self.out_sequences = list()
		self.in_sequences = list()
		self.alignments = dict()

try:
	with open(gfa2Filename, 'r') as infile:
		for line in infile:
			fields = line.split()
			if fields[0] == "S":
				name = fields[1]
				sequence = fields[3]
				sequences[name] = Sequence(name, sequence)
			elif fields[0] == "E":
				from_sequence = fields[2].split("+")[0]
				to_sequence = fields[3].split("+")[0]
				sequences[from_sequence].out_sequences.append(sequences[to_sequence])
				sequences[to_sequence].in_sequences.append(sequences[from_sequence])

except (OSError, IOError) as e: 
	if getattr(e, 'errno', 0) == errno.ENOENT:
		print("Could not find file " + gfa2Filename)
	print(e)
	sys.exit(2)

try:
	with open(blastFilename, 'r') as infile:
		for line in infile:
			fields = line.split()
			query = fields[0]
			reference = fields[1]
			coverage = float(fields[3])/float(fields[4])
			identity = float(fields[2])/100

			if reference not in references:
				references.append(reference)
			sequences[query].alignments[reference] = Alignment(query, reference, coverage, identity)

except (OSError, IOError) as e: 
	if getattr(e, 'errno', 0) == errno.ENOENT:
		print("Could not find file " + blastFilename)
	print(e)
	sys.exit(2)

count = 0
for startingSequence in sequences.values():
	if len(startingSequence.in_sequences) == 0:
		for ref in references:
			current_sequence = startingSequence
			assembly = ""
			while len(current_sequence.out_sequences) > 0:
				assembly += current_sequence.sequence
				if len(current_sequence.out_sequences) == 1:
					current_sequence = current_sequence.out_sequences[0]
				else:
					best_score = 0
					best_seq = 0
					for next_seq in current_sequence.out_sequences:
						if ref in next_seq.alignments.keys() and next_seq.alignments[ref].score > best_score:
							best_score = next_seq.alignments[ref].score
							best_seq = next_seq
					if best_seq == 0:
						# in this case there are no mappings from ref to any of the out sequences. Pick the first which is guaranteed to exist.
						best_seq = current_sequence.out_sequences[0]
					current_sequence = best_seq
			assembly += current_sequence.sequence
			with open(ref + "_gfa_assembled.fa", 'w') as outfile:
				outfile.write(">" + ref + "_" + str(count) + "\n")
				outfile.write(assembly + "\n")
		count += 1

























