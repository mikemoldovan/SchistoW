"""
Bachtrog -- step2.2
given a coverage file, obtain sequences that didn't get any coverage
"""

from Bio import SeqIO
from optparse import OptionParser

def mk_covered_dict(covfile):
	covered_dict = dict()
	for s in open(covfile):
		s = s.strip().split()
		covered_dict[s[0]] = True
	return covered_dict

def fastaparse(fafile, covered_dict):
	singlets = dict()
	with open(fafile) as handle:
		for record in SeqIO.parse(handle, "fasta"):
			try:
				temp = covered_dict[record.id]
			except:
				singlets[record.id] = 0.0
	for k in singlets.keys():
		print k, 0.0 
	return singlets

def make_covdict(coverage_file):
	covdict = dict()
	with open(coverage_file) as covfile:
		for s in covfile:
			s = s.strip().split()
			try:
				covdict[s[0]][eval(s[1])] = eval(s[2])
			except:
				covdict[s[0]] = dict()
				covdict[s[0]][eval(s[1])] = eval(s[2])
	return covdict

def make_percent_covdict(covdict):
	percent_covdict = dict()
	for k in covdict.keys():
		cov = 0
		n = 0
		for i in covdict[k].keys():
			if covdict[k][i] != 0:
				cov += 1
			n += 1
		percent_covdict[k] = float(cov)/n
		print k, float(cov)/n
	return percent_covdict

def fastasample(infasta_name, outfasta_name, percent_covdict, threshold):
	outfasta = open(outfasta_name, "w")
	with open(infasta_name) as handle:
		for record in SeqIO.parse(handle, "fasta"):
			cov = percent_covdict[record.id]
			if cov <= threshold:
				SeqIO.write(record, outfasta, "fasta")
	outfasta.close()


def main(merged_covfile, threshold, infasta_name, outfasta_name):
#find sequences, onto which no male reads mapped
	covered_dict = mk_covered_dict(merged_covfile)
	singlets = fastaparse(infasta_name, covered_dict)
#Calculate coverage of the rest
	covdict = make_covdict(merged_covfile)
	percent_covdict = make_percent_covdict(covdict)
	percent_covdict.update(singlets)
#Filter fastafile
	fastasample(infasta_name, outfasta_name, percent_covdict, threshold)

parser = OptionParser()
parser.add_option("-c", "--merged_covfile", help="Coverage file of merged BAM files")
parser.add_option("-r", "--threshold", help="Coverage threshold")
parser.add_option("-i", "--infasta_name", help="Name of the scaffold fasta")
parser.add_option("-o", "--outfasta_name", help="Name of the output fasta")
opt, args = parser.parse_args()


main(opt.merged_covfile, eval(opt.threshold), opt.infasta_name, opt.outfasta_name)