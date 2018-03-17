"""
Bactchtrog's pipeline -- step3
Map transcripts to male reference genome assembly
Discard transcripts that map
Takes info from the PSL BLAT output and selects FASTA Trinity sequences, that didn't map
Input:
1. PSL BLAT output
2. FASTA Trinity query
3. Mapping length threshold
4. Outfile name
Output:
1. Filtered FASTA with transcripts
"""

from Bio import SeqIO
from optparse import OptionParser

def psl_dict_build(psl_file, maplen_thresh):
	psl_dict = dict()
	with open(psl_file) as psl:
		for s in psl:
			s = s.strip().split()
			if len(s) < 1:
				continue
			if not s[0].isdigit():
				continue
			match = eval(s[0])
			q_size = eval(s[10])
			maplen = float(match)/q_size
			if maplen >= maplen_thresh:
				psl_dict[s[9]] = True
	return psl_dict


def fastaparse(trinity_fasta, psl_dict, outfile_name):
	outhandle = open(outfile_name)
	with open(trinity_fasta, "rU") as handle:
		for record in SeqIO.parse(handle, "fasta"):
			try:
				temp = psl_dict[record.id]
			except:
				SeqIO.write(record, outhandle, "fasta")
	outhandle.close()

parser = OptionParser()
parser.add_option("-p", "--psl_file", help="PSL BLAT output")
parser.add_option("-o", "--outfile_name", help="Outfile name")
parser.add_option("-i", "--trinity_fasta", help="FASTA Trinity query")
parser.add_option("-l", "--maplen_thresh", help="Mapping length threshold (default 0.9)", default="0.9")
opt, args = parser.parse_args()

psl_dict = psl_dict_build(opt.psl_file, eval(opt.maplen_thresh))
fastaparse(opt.trinity_fasta, psl_dict, opt.outfile_name)

