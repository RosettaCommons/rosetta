#!/usr/bin/env python
import os,sys,re,string
#from AshworthUtil import get_oneletter

from optparse import OptionParser
p = OptionParser()
p.add_option('-m','--matrix',help='amino acid substitution matrix')
p.add_option('-p','--pdb',help='pdb file')
p.add_option('-w','--weight',type='float',default=1.0,help='weighting factor by which matrix scores are multiplied')
opt,args = p.parse_args()

oneletter = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G',
             'HIS':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N',
             'PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V',
             'TRP':'W','TYR':'Y',
						 'ADE':'a','CYT':'c','GUA':'g','THY':'t',
						 '  A':'a','  C':'c','  G':'g','  T':'t',
						 'A':'a','C':'c','G':'g','T':'t',
						 'DA':'a','DC':'c','DG':'g','DT':'t', # new pdb spec
						 ' DA':'a',' DC':'c',' DG':'g',' DT':'t'} # new pdb spec

def get_oneletter(aa3):
	if oneletter.has_key(aa3): return oneletter[aa3]
	return aa3.strip()[0]

# initialization of lookup matrix
# note: this list must be of length 20 and its order must match Rosetta aa ordering
aas = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
substitution_matrix = {}
for aa1 in aas:
	for aa2 in aas:
		substitution_matrix[ (aa1,aa2) ] = None

if opt.matrix:
	# read input matrix
	# assumes full 2D table format
	column_ids = []
	for line in open(opt.matrix):
		if line.startswith('#'): continue
		fields = line.strip().split()

		# first line assumed to identify columns
		if column_ids == []:
			column_ids = fields
			if len(column_ids) != 20:
				sys.stderr.write('error wrong # of columns in %s\n' %line); sys.exit()
			for letter in column_ids:
				if not letter in aas:
					sys.stderr.write('error unknown aa %s\n' %letter); sys.exit()
			continue

		if len(fields) != 21:
			sys.stderr.write('error wrong number of columns in %s\n' %line); sys.exit()
		aa = fields[0]
		if not aa in aas:
			sys.stderr.write('error unknown aa %s\n' %letter); sys.exit()
		for i in range(20):
			aa2 = column_ids[i]
			if substitution_matrix[ (aa,aa2) ] != None:
				sys.stderr.write('error, shouldn\'t be overwriting substitution score for (%s,%s)\n' %(aa,aa2)); sys.exit()
			substitution_matrix[ (aa,aa2) ] = float( fields[i+1] )

# verify/validate matrix
#print 'Verifying substitution matrix:'
for aa1 in aas:
	for aa2 in aas:
		if not substitution_matrix.has_key( (aa1,aa2) ):
			sys.stderr.write('error, substitution matrix missing key (%s,%s)\n' %(aa1,aa2)); sys.exit()
		score = substitution_matrix[(aa1,aa2)]
		if score == None:
			sys.stderr.write('error, missing score for %s->%s\n' %(aa1,aa2)); sys.exit()
#		print '%s->%s: %.0f' %(aa1,aa2,score)

# matches pdb c-alpha line
re_ca = re.compile( 'ATOM.+?CA  ([A-Z0-9 ]{3}) [A-Z]')

pdb = open(opt.pdb).read()

prof_file = open('sequence_profile','w')
constraint_file = open('constraints','w')

prof_file.write('made for pdb %s from matrix %s with factor %g. (note: SequenceProfile::read_from_checkpoint() skips this line)\n' %(opt.pdb,opt.matrix,opt.weight))
index = 1
for match in re_ca.finditer(pdb):
	constraint_file.write('SequenceProfile %i sequence_profile\n' %index)
	aa3 = match.groups()[0]
	aa1 = get_oneletter(aa3)
	prof_file.write(aa1)
	for aa2 in aas:
		prof_file.write( ' %g' %( opt.weight * substitution_matrix[(aa1,aa2)] ) )
	prof_file.write('\n')
	index += 1
