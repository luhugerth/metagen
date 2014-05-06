#!/usr/bin/env python

from Bio import SeqIO
import argparse
import csv

def main():
	parser = argparse.ArgumentParser(description='Output the sequences from each Concoct bin to a different file')
	parser.add_argument('concoct', help='csv associating contigs to bins')
	parser.add_argument('fasta', help='fasta file that was used as input for concoct')
	parser.add_argument('prefix', help='prefix for the output bin files: prefix-binNr.fasta')
	args = parser.parse_args()

	concoct_file=args.concoct
	contigbins = {}
	with open(concoct_file) as csvfile:
		binreader = csv.reader(csvfile)
		for row in binreader:
			thebin=row[1]
			contig=row[0]
			data=contig.split('__') ##remove extra information, if added
			contig=data[0]
			contigbins[contig] = thebin
				
	prefix=args.prefix
	infile=args.fasta
	records = SeqIO.parse(open(infile), "fasta")
	for seq in records:
		thisid=seq.id
		data=thisid.split('__') ##remove extra information, if added
		thisid=data[0]
		if thisid in contigbins:
			thisbin=contigbins[thisid]
			outname = prefix + "-bin" + str(thisbin) + ".fasta"
			handle = open(outname, "a")
			SeqIO.write(seq, handle, "fasta")
		
if __name__ == '__main__':
    main()
