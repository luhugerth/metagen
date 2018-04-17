#!/usr/bin/env python

from Bio import Phylo
import argparse
import csv

def parse_taxonomy(taxfile, keeplist):
	keepers = keeplist.split(",")
	tax_dict = dict()
	with open(taxfile) as csvfile:
		reader = csv.reader(csvfile, delimiter="\t")
		for row in reader:
			ID = row[0]
			fulltax = row[1].split(";")
			taxa = []
			for keep in keepers:
				tax = fulltax[int(keep)]
				if tax is not None:
					taxa.append(tax)
				else:
					taxa.append('unknown')
			tax_dict[ID] = ";".join(taxa)
	return(tax_dict)

def convert_tree(taxa, treefile, outfile):
	tree = Phylo.read(treefile, 'newick')
	for leaf in tree.get_terminals(): 
		if leaf.name in taxa:
			leaf.name = taxa[leaf.name]
	Phylo.write(tree, outfile, 'newick')


def main(taxfile, keepers, tree, outfile):
	taxtab = parse_taxonomy(taxfile, keepers)
	convert_tree(taxtab, tree, outfile)

		
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Takes a tree in newick format and replaces node names with taxonomy, where available')
	parser.add_argument('-tax', '--taxonomy', help='tsv with taxon ID on the left and ;-separeted taxonomy on the right')
	parser.add_argument('-k', '--keep', default='1,5,6', help='comma-separated list of fields in the taxonomy to keep (0-based)')
	parser.add_argument('-tre', '--treefile',  help='Phylogenetic tree in newick format')
	parser.add_argument('-o', '--outfile', default='named_tree.nwk', help="Output file")

	args = parser.parse_args()

	main(args.taxonomy, args.keep, args.treefile, args.outfile)
