"""
Step6 -- compare coverages
Takes two coverage files and returns the file with each scaffold's coverage in percentage
Input:
1. Covfile1 -- heterogametic coverage -- can't be 0
2. Covfile2
3. Threshold1 -- lower for heterogametic
4. Threshold2 -- upper for homogametic
5. Infasta name
6. Outfasta name
Output:
file with merged coverages
"""

from optparse import OptionParser
from Bio import SeqIO

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

def compare_covdicts(covdict1, covdict2, threshold1, threshold2):
	filtered_dict = dict()
	for k in covdict1.keys():
		if covdict1[k] < threshold1:
			continue
		try:
			val_hom = covdict2[k]
			if val_hom <= threshold2:
				filtered_dict[k] = {1:covdict1[k],2:covdict2[k]}
		except:
			filtered_dict[k] = {1:covdict1[k],2:0.0}
	for k in filtered_dict.keys():
		print k, filtered_dict[k][1], filtered_dict[k][2]
	return filtered_dict

def fastasample(infasta_name, outfasta_name, filtered_dict):
	outfasta = open(outfasta_name, "w")
	with open(infasta_name) as handle:
		for record in SeqIO.parse(handle, "fasta"):
			try:
				a = filtered_dict[record.id]
				SeqIO.write(record, outfasta, "fasta")
				print "sampled",record.id, a[1], a[2]
			except:
				pass
	outfasta.close()

def main(covfile1, covfile2, thresh1, thresh2, infasta_name, outfasta_name):
	covdict1 = make_percent_covdict(make_covdict(covfile1))
	covdict2 = make_percent_covdict(make_covdict(covfile2))
	filtered_dict = compare_covdicts(covdict1, covdict2, thresh1, thresh2)
	fastasample(infasta_name, outfasta_name, filtered_dict)

parser = OptionParser()
parser.add_option("-a", "--covfile1", help="Covfile1 -- heterogametic coverage -- can't be 0")
parser.add_option("-b", "--covfile2", help="Homogametic coverage")
parser.add_option("-1", "--thresh1", help="Threshold1 -- lower for heterogametic")
parser.add_option("-2", "--thresh2", help="Threshold2 -- upper for homogametic")
parser.add_option("-i", "--infasta_name", help="Name of the scaffold fasta")
parser.add_option("-o", "--outfasta_name", help="Name of the output fasta")
opt, args = parser.parse_args()

main(opt.covfile1, opt.covfile2, eval(opt.thresh1), eval(opt.thresh2), opt.infasta_name, opt.outfasta_name)
