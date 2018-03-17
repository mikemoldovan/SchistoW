"""
Step 5: Merge and extend transcripts, parse TGICL results
Input:
1. TGICL generated *_cl_clusters file
2. TGICL generated asm_1/contigs file
3. Initial transcripts FASTA
4. Name of the output FASTA
Output:
FASTA with unique transcripts
"""

from optparse import OptionParser
from Bio import SeqIO

def clustfile_parse(clustfile):
	id_dict = dict()
	for s in open(clustfile):
		s = s.strip().split()
		if s[0][0] == '>':
			clust_id = s[0][1:]
			continue
		id_dict[clust_id] = s
	return id_dict

def make_new_fasta(clustfile, tgicl_contigs, init_fasta, outfasta_name):
	id_dict = clustfile_parse(clustfile)
	outhandle = open(outfasta_name, "w")
	discarded_ids = []
	for record in SeqIO.parse(tgicl_contigs, "fasta"):
		clust_id = record.id[:-7]
		discarded_ids += id_dict[clust_id]
		SeqIO.write(record, outhandle, "fasta")
	for record in SeqIO.parse(init_fasta, "fasta"):
		if record.id not in discarded_ids:
			SeqIO.write(record, outhandle, "fasta")
	outhandle.close()


parser = OptionParser()
parser.add_option("-c", "--clustfile", help="TGICL generated *_cl_clusters file")
parser.add_option("-a", "--tgicl_contigs", help="TGICL generated asm_1/contigs file")
parser.add_option("-i", "--init_fasta", help="Initial transcripts FASTA")
parser.add_option("-o", "--outfasta_name", help="Name of the output FASTA")
opt, args = parser.parse_args()

make_new_fasta(opt.clustfile, opt.tgicl_contigs, opt.init_fasta, opt.outfasta_name)
"""
def clustfile_parse(clustfile):
	id_dict = dict()
	for s in open(clustfile):
		if s[0] == '>':
			continue
		s = s.strip().split()
		for _id in s:
			id_dict[_id] = True
	return id_dict


def make_new_fasta(clustfile, tgicl_contigs, init_fasta, outfasta_name):
	id_dict = clustfile_parse(clustfile)
	outhandle = open(outfasta_name, "w")
	for record in SeqIO.parse(tgicl_contigs, "fasta"):
		SeqIO.write(record, outhandle, "fasta")
	for record in SeqIO.parse(init_fasta, "fasta"):
		try:
			temp = id_dict[record.id]
			print record.id
		except:
			SeqIO.write(record, outhandle, "fasta")
	outhandle.close()
"""