#!/usr/bin/env python

from __future__ import division
import argparse
import csv
import sys

csv.field_size_limit(sys.maxsize)

def main():
	parser = argparse.ArgumentParser(description='Associate taxonomy to bins')
	parser.add_argument('concoct', help='csv associating contigs to bins')
	parser.add_argument('names', help='tsv associating tax id with its name')
	parser.add_argument('classification', help='kraken output file')
	parser.add_argument('taxonomy', help='tsv table giving ;-separated full phylogeny of kraken taxa')
	parser.add_argument('minlen', help='minimum length used in concoct', type=int)
	parser.add_argument('minprop', help='minimum support to keep a classification', type=float)
	args = parser.parse_args()

	## Read concoct file associating contigs to bins
	concoct_file=args.concoct
	contigbins = {}
	with open(concoct_file) as csvfile:
		binreader = csv.reader(csvfile)
		for row in binreader:
			thebin=row[1]
			contig=row[0]
			data=contig.split('__') ##remove extra information, if added
			contig=data[0]
			if thebin in contigbins:
				contigbins[thebin].append(contig)
			else:
				contigbins[thebin] = [contig]

	## Read kraken reference associating numbers to taxa
	kraken_tab=args.names
	taxids = {}
	taxids['0'] = 'Unclassified'
	with open(kraken_tab) as csvfile:
		taxreader = csv.reader(csvfile, delimiter="\t")
		for row in taxreader:
			taxids[row[0]]  = row[1]

	## Read kraken reference associating taxa to phylogeny
	kraken_tax=args.taxonomy
	fulltaxa = {}
	fulltaxa['0'] = '0'
	with open(kraken_tax) as csvfile:
		taxreader = csv.reader(csvfile, delimiter="\t")
		for row in taxreader:
			fulltaxa[row[0]]  = row[1]
	
	## Read kraken output files associating contigs to taxa numbers
	kraken=args.classification
	minlength=args.minlen
	contigclass = {}
	with open(kraken) as csvfile:
		krakreader = csv.reader(csvfile, delimiter="\t")

		for row in krakreader:		
			contig = row[1]
			length = int(row[3])
			taxa = row[4]
			alltaxa = taxa.split(' ')
			if length > minlength:	
#				contig = contigs[contig]

				for value in alltaxa:
					pair = value.split(':')
					thistax = pair[0]
					length = int(pair[1])
					fulltax = fulltaxa[thistax]
					eachtax = fulltax.split(';')
					depth = 0

					for taxid in eachtax:
						depth += 1
						taxname = taxids[taxid]
						if contig in contigclass:
							if depth in contigclass[contig]:
								if taxname in contigclass[contig][depth]:
									contigclass[contig][depth][taxname] += length
								else:
									contigclass[contig][depth][taxname] = length
							else:
									contigclass[contig][depth] =  {taxname:length}
						else:
							contigclass[contig] = {depth:{taxname:length}}

	## Output the per-kmer classification for each CONTIG in each bin
	for mybin in contigbins.keys():
		thesecontigs = contigbins[mybin]
		for mycontig in thesecontigs:
			totalkmer = 0
			outclass = ''
			for mydepth in contigclass[mycontig].keys():
				maxkmer = 0
				if mydepth == 1:
					for tax in contigclass[mycontig][mydepth].keys():
						kmer = contigclass[mycontig][mydepth][tax]
						totalkmer += kmer
						if kmer > maxkmer:
							toptax = tax
							maxkmer = kmer	
					maxprop = maxkmer/totalkmer
					outclass = outclass + toptax + ';' + str(maxprop) + ';'
				else:
					for tax in contigclass[mycontig][mydepth].keys():
						kmer = contigclass[mycontig][mydepth][tax]
						if kmer > maxkmer:
							toptax = tax
							maxkmer = kmer	
					maxprop = maxkmer/totalkmer
					outclass = outclass + toptax + ';' + str(maxprop) + ';'

			print "%s\t%s\t%s\t" %(mybin, mycontig, outclass)

if __name__ == '__main__':
    main()
