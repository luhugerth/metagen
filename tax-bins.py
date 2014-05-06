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
	parser.add_argument('minprop', help='minimum support to keep a classification' type=float)
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
			contigbins[contig]  = thebin

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
	
	## Read kraken output files associating contigs to taxa numbers and store PER BIN
	kraken=args.classification
	minlength=args.minlen
	binclass = {}
	with open(kraken) as csvfile:
		krakreader = csv.reader(csvfile, delimiter="\t")

		for row in krakreader:		
			contig = row[1]
			##upstream software may insert extra information which we strip here
			allinfo=contig.split('__')
			contig=allinfo[0]
			length = int(row[3])
			taxa = row[4]
			alltaxa = taxa.split(' ')
			if length > minlength:	
				contigbin = contigbins[contig]

				for value in alltaxa:
					pair = value.split(':')
					thistax = pair[0]
					length = int(pair[1])
					fulltax = fulltaxa[thistax]
					eachtax = fulltax.split(';')
					depth = 0

					for taxid in eachtax:
						depth += 1
						#print taxid
						taxname = taxids[taxid]
						if contigbin in binclass:
							if depth in binclass[contigbin]:
								if taxname in binclass[contigbin][depth]:
									binclass[contigbin][depth][taxname] += length
								else:
									binclass[contigbin][depth][taxname] = length
							else:
									binclass[contigbin][depth] =  {taxname:length}
						else:
							binclass[contigbin] = {depth:{taxname:length}}

	## Output the per-kmer classification for each BIN
	for mybin in binclass.keys():
		totalkmer = 0
		outclass = ''
		for mydepth in binclass[mybin].keys():
			maxkmer = 0
			if mydepth == 1:
				for tax in binclass[mybin][mydepth].keys():
					kmer = binclass[mybin][mydepth][tax]
					totalkmer += kmer
					if kmer > maxkmer:
						toptax = tax
						maxkmer = kmer	
										
				maxprop = maxkmer/totalkmer
				if maxprop >= minprop:
					outclass = outclass + toptax + ';' + str(maxprop) + ';'
			else:
				for tax in binclass[mybin][mydepth].keys():
					kmer = binclass[mybin][mydepth][tax]
					if kmer > maxkmer:
						toptax = tax
						maxkmer = kmer	
								
				maxprop = maxkmer/totalkmer
				outclass = outclass + toptax + ';' + str(maxprop) + ';'
		print "%s\t%s" %(mybin, outclass)

if __name__ == '__main__':
    main()
