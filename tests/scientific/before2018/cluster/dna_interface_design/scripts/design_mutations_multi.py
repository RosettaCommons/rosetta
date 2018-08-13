#!/usr/bin/python
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

# ashworth
import sys,string
from AshworthUtil import *
from optparse import OptionParser

p = OptionParser()
p.add_option( '-f', '--filelist',help='file listing pdb paths' )
p.add_option( '-p', '--pdbpath', default='./', help='directory path where output pdbs are located' )
p.add_option( '-n', '--natives', help='directory path where input pdbs are located' )
p.add_option( '-r', '--rms', type='float', default=1.0, help='RMS cutoff value, for assesing wildtype rotamer recovery' )
p.add_option( '-o', '--out', default='design_mutations', help='output filename' )
opt,args = p.parse_args()

if not opt.pdbpath or not opt.natives:
	p.print_help(); sys.exit()

typeconversions = {
	' DA':'  A',
	' DC':'  C',
	' DG':'  G',
	' DT':'  T',
}

def officialtype(type):
	if typeconversions.has_key(type): return typeconversions[type]
	return type

################################################################################
# output distribution and recovery statistics

def statsfileheader():
	outstring = ''
	outstring += 'Statistics for residues that were allowed to mutate:\n'
	outstring += 'RMS cutoff for rotamer conservation (\'cnsrot\') in mutable residues: %4.2f\n' % opt.rms
	outstring += '-------------------------------------------------------------\n'
	outstring += 'AA  : %s %s\n' % ('--------Composition--------','----------Recovery---------')
	outstring += 'AA  : %6s %6s %6s %6s %6s %6s %6s %6s\n' % ('nat','f','des','f','cons','f','cnsrot','f')
	outstring += '-------------------------------------------------------------\n'
	return outstring

def statscsvfileheader():
	outstring = ''
	outstring += 'AA,counts_native,freq_native,count_design,freq_design,counts_recov,freq_recov,counts_rot,freq_rot\n'
	return outstring

def output_stats(aachoices,full=True):

	types = []

	typesets = {
		'protein' : ['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR'],
		'dna' : ['  A','  C','  G','  T'],
	}

	present = { 'protein' : False, 'dna' : False }
	for type,typeset in typesets.items():
		for t in typeset:
			if aachoices.has_key(t):
				present[type] = True
				break
		if present[type]: types.extend(typeset)

	total_counts = 0
	counts = {}
	for mode in ['wildtype','design','conserved','conserved_rotamer']:
		counts[mode] = {}
		for aa in types:
			counts[mode][aa] = 0

	submatrix = {}
	for aa1 in types:
		for aa2 in types:
			submatrix[ (aa1,aa2) ] = 0

	for aa in types:
		if not aachoices.has_key(aa): continue
		total_counts += len(aachoices[aa]) # will sum to total number of designable positions
		counts['wildtype'][aa] = len(aachoices[aa]) # for wildtype composition
		for choice in aachoices[aa]:
			des_aa = choice[0]
			submatrix[ (aa,des_aa) ] += 1
			counts['design'][des_aa] += 1 # for design composition
			if des_aa == aa:
				counts['conserved'][aa] += 1
				if choice[1]: counts['conserved_rotamer'][aa] += 1

	total_conserved = 0
	total_conserved_rotamers = 0

	outstrings = {'HR':'','csv':''} # HR: 'human-readable'

	for aa in types:
		n_wt  = counts['wildtype'][aa]
		n_des = counts['design'][aa]
		n_con = counts['conserved'][aa]
		n_rot = counts['conserved_rotamer'][aa]

		p_nat = '0'; p_des = '0'
		# composition of wildtype
		if total_counts != 0: p_nat = '%.3f' %( float(n_wt) / float(total_counts) )
		# composition of designs
		if total_counts != 0: p_des = '%.3f' %( float(n_des) / float(total_counts) )

		p_con = '0'; p_rot = '0'
		# conservation of wildtype by designs
		if n_wt != 0: p_con = '%.3f' %( float(n_con) / float(n_wt) )
		# conservation of wildtype rotamer by designs
		if n_wt != 0: p_rot = '%.3f' %( float(n_rot) / float(n_wt) )

		if full:
			outstrings['HR'] += '%s : %6i %6s %6i %6s %6i %6s %6i %6s\n' \
				%( aa, n_wt, p_nat, n_des, p_des, n_con, p_con, n_rot, p_rot )
			outstrings['csv'] += '%s,%i,%s,%i,%s,%i,%s,%i,%s\n' \
				%( aa, n_wt, p_nat, n_des, p_des, n_con, p_con, n_rot, p_rot )
		total_conserved += n_con
		total_conserved_rotamers += n_rot

	p_con = '0'; p_rot = '0'
	if total_counts != 0: p_con = '%.3f' % ( float(total_conserved) / float(total_counts) )
	if total_counts != 0: p_rot = '%.3f' % ( float(total_conserved_rotamers) / float(total_counts) )

	if full:
		cat = 'TOT'
		if       present['protein']     and present['dna']: cat = 'ALL'
		elif     present['protein'] and not present['dna']: cat = 'PRT'
		elif not present['protein']     and present['dna']: cat = 'DNA'

		outstrings['HR'] += '-------------------------------------------------------------\n'
		outstrings['HR'] += '%s : %6i %6.3f %6i %6.3f %6i %6s %6i %6s\n' \
			%(cat, total_counts, 1, total_counts, 1, total_conserved, p_con, total_conserved_rotamers, p_rot )
		outstrings['csv'] += '%s,%i,%.3f,%i,%.3f,%i,%s,%i,%s\n' \
			%(cat, total_counts, 1, total_counts, 1, total_conserved, p_con, total_conserved_rotamers, p_rot )

	else: outstrings['HR'] += p_con
	return outstrings

def output_submatrix(aachoices):
	outstring = ''
	submatrix = {}
	aas = aachoices.keys()
	aas.sort()
	for aa1 in aas:
		for aa2 in aas:
			submatrix[ (aa1,aa2) ] = 0

	for aa in aas:
		if not aachoices.has_key(aa): continue
		for choice in aachoices[aa]:
			des_aa = choice[0]
			if not submatrix.has_key( (aa,des_aa) ): submatrix[ (aa,des_aa) ] = 0
			submatrix[ (aa,des_aa) ] += 1

	outstring += '%s\n' %string.join( [ aa for aa in aas ], ',' )
	for aa1 in aas:
		outstring += aa1
		for aa2 in aas:
			outstring += ',%i' %submatrix[(aa1,aa2)]
		outstring += '\n'
	return outstring

#############
### START ###
#############

################################################################################
# load designs
designs = {}
if opt.filelist:
	if not os.path.exists(opt.filelist):
		print "error: %s not found" %opt.filelist
	for f in open(opt.filelist):
		filepath = f.strip()
		root = f.split('/')[-1][0:4]
		if not designs.has_key(root): designs[root] = []
		designs[root].append( RosettaPDBFile('f',filepath,True,['energy','molten']) )
else:
	for filepath in os.listdir(opt.pdbpath):
		if not re.search('\.pdb',filepath): continue
		root = filepath[0:4]
		if not designs.has_key(root): designs[root] = []
		designs[root].append( RosettaPDBFile('f','%s/%s'%(opt.pdbpath,filepath),True,['energy','molten']) )

################################################################################
# load reference pdbs
ref_pdbs = {}
for root in designs.keys():
	ref_path = '%s/%s.pdb' %(opt.natives,root)
	if not os.path.exists(ref_path):
		ref_path = '%s.gz' %ref_path
		if not os.path.exists(ref_path):
			print 'warning, can\'t find reference %s' %ref_path
	if not ref_pdbs.has_key(root):
		ref_pdbs[root] = RosettaPDBFile('f',ref_path,True,['energy'])

################################################################################
# despos contains the total set of designable residues for each given scaffold
despos = {}
for root,ref_pdb in ref_pdbs.items():
	if not despos.has_key(root): despos[root] = {}
	for chainpos,resn in ref_pdb.residues.items():
		vary = False
		for design in designs[root]:
			if chainpos in design.mutable():
				vary = True; break
		if vary: despos[root][chainpos] = []

################################################################################

mutationsfile = open( '%s.mutations' %opt.out, 'w' )
statsfile = open( '%s.stats' %opt.out, 'w' )
statscsvfile = open( '%s.stats.csv' %opt.out, 'w' )
fastafile = open( '%s.fasta' %opt.out, 'w' )
for ref in ref_pdbs.values():
	fastafile.write('%s\n' %ref.fasta())
for deslist in designs.values():
	for des in deslist:
		fastafile.write('%s\n' %des.fasta())
fastafile.close()
rmsfile = open('%s.rms' % opt.out, 'w' )

rmsfile.write('RMS from native sidechain if mutable residue remained wildtype.\n')
rmsfile.write('(for assessing wildtype recovery when designed against wildtype DNA sequence)\n')

################################################################################
#### main loop ####
# count aa types in native and designs for compositional statistics
choices = {}
last_root = ''
sorted_roots = designs.keys()
sorted_roots.sort()
for root in sorted_roots:
#	print root
	if not choices.has_key(root): choices[root] = {}

################################################################################
	# write reference lines
	if root != last_root:
#		print 'new root:',root
		if choices.has_key(last_root):
			thispdbstats = output_stats(choices[last_root],False)['HR']
#			print 'stats for old root:',root,thispdbstats
			statsfile.write('%s recovery: %s\n' %( last_root, thispdbstats ) )
		last_root = root

		# write reference lines for new root
		mutationsfile.write( '\n%50s' % 'pdb' )
		despos_thisroot = despos[root].keys()
		despos_thisroot.sort()
		for chainpos in despos_thisroot: mutationsfile.write( '%6s' %str(chainpos) )
		mutationsfile.write('\n')
		# write native info as header
		mutationsfile.write( '%50s' % opt.natives.split('/')[-1] )
		for chainpos in despos_thisroot:
			mutationsfile.write( '%6s' % ref_pdbs[root].residues[chainpos] )
		mutationsfile.write( '%9.2f\n' % ref_pdbs[root].E['energy'] )

	# loop over designs for this root
	for design in designs[root]:
		rmsfile.write( '%s\n' % design.filename )
		mutationsfile.write( '%50s' % design.filename )

	################################################################################
		# print conserved and mutated residues in despos
		despos_thisroot = despos[root].keys()
		despos_thisroot.sort()
		for chainpos in despos_thisroot:
			res = ' . '
			des_res = ' . '
			if design.residues.has_key(chainpos):
				des_res = officialtype( design.residues[chainpos].type )
				if des_res == 'UNK': des_res = ' . '
			nat_res = officialtype( ref_pdbs[root].residues[chainpos].type )
			if not des_res == nat_res: res = des_res
			mutationsfile.write( '%6s' % res )

		mutationsfile.write( '%9.2f\n' % design.E['energy'] )

	################################################################################
		# compile list of substitutions for distribution and recovery statistics
		# only count at mutable positions within each design--this could be a SUBSET of ref's despos!
		for chainpos in design.mutable():
			des_res = officialtype( design.residues[chainpos].type )
			if not ref_pdbs[root].residues.has_key(chainpos):
				print 'error: can\'t find chainpos %s in reference pdb %s for root %s' %(str(chainpos),ref_pdbs[root].filename,root)
			nat_res = officialtype( ref_pdbs[root].residues[chainpos].type )
			if not choices[root].has_key(nat_res): choices[root][nat_res] = []
			# check for rotamer conservation if same as wildtype
			rotcons = False
			# this may occur with sparse pdb output: pdb file indicated this was a design position, but no coordinates were output for it, because it remained wildtype (w/ wt rotamer)
			if des_res == 'UNK':
#				print 'warning, mutable residue %s \'UNK\' assumed to be wildtype' %chainpos
				des_res = nat_res
				rotcons = True
			if des_res == nat_res:
				if rotcons == False:
					rmsd = residue_rmsd( design.residues[chainpos], ref_pdbs[root].residues[chainpos] )
					rmsfile.write('\t%6s%4s%7.3f\n' %(str(chainpos),des_res,rmsd) )
					if rmsd <= opt.rms: rotcons = True
				despos[root][chainpos].append(1)
			else:
				despos[root][chainpos].append(0)
			choices[root][nat_res].append( (des_res,rotcons) )

if choices.has_key(last_root):
	statsfile.write('%s recovery: %s\n' %( last_root, output_stats(choices[last_root],False)['HR'] ) )

mutationsfile.close()
rmsfile.close()

desposfile = open( '%s.despos' %opt.out, 'w' )
desposfile.write('root,chainpos,recrate,n\n')
roots = despos.keys()
roots.sort()
for r in roots:
	cps = despos[r].keys()
	cps.sort()
	for cp in cps:
		dps = despos[r][cp]
		recrate = float(sum(dps))/float(len(dps))
		s = '%s,%s,%.3f,%i' %(r,cp,recrate,len(dps))
#		s = '%s,%s,%s' %(r,cp,string.join(['%i' %dp for dp in dps],','))
		desposfile.write('%s\n' %s)
desposfile.close()

combined_choice_lists = {}
for aachoices in choices.values():
	for aa,choices in aachoices.items():
		if not combined_choice_lists.has_key(aa): combined_choice_lists[aa] = []
		combined_choice_lists[aa].extend(choices)

results = output_stats(combined_choice_lists,True)

statsfile.write( statsfileheader() )
statsfile.write( results['HR'] )
statsfile.close()

statscsvfile.write( statscsvfileheader() )
statscsvfile.write( results['csv'] )
statscsvfile.close()

submatrixfile = open( '%s.matrix' %opt.out, 'w' )
submatrixfile.write( output_submatrix(combined_choice_lists) )
submatrixfile.close()
