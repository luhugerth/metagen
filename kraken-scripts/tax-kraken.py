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
	contigsupport = {}
	contignull = {}
	with open(kraken) as csvfile:
		krakreader = csv.reader(csvfile, delimiter="\t")

		for row in krakreader:		
			contig = row[1]
			totlength = int(row[3])

			winner = row[2]
			ancestors = fulltaxa[winner]
			family = ancestors.split(';')

			taxa = row[4]
			alltaxa = taxa.split(' ')

			if totlength > minlength:	

				unclass = 0
				support = 0

				for value in alltaxa:
					pair = value.split(':')
					thistax = pair[0]
					length = int(pair[1])
					
					if thistax == '0':
						unclass += length
					elif thistax in family:
						support += length
				
				contigclass[contig] = winner
				contigsupport[contig] = support/totlength
				contignull[contig] = unclass/totlength


	print "Bin\tContig\tClassification\tSupport\tUnclassified"
	for mybin in contigbins.keys():
		thesecontigs = contigbins[mybin]
		for mycontig in thesecontigs:
			winner = contigclass[mycontig]
			taxon = taxids[winner]
			support = contigsupport[mycontig]
			null = contignull[mycontig]

##add: print the whole taxonomy
			print "%s\t%s\t%s\t%f\t%f" %(mybin, mycontig, taxon, support, null)

if __name__ == '__main__':
    main()
