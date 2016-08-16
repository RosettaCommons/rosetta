#!/usr/bin/env python

import string
import sys
import os
from sys import argv,stderr,stdout
from os import popen,system
from os.path import exists,basename
from amino_acids import longer_names
from amino_acids import modres
from string import uppercase
import gzip
import random


# remote host for downloading pdbs
remote_host = 'ws0.nrb'

shit_stat_insres = False
shit_stat_altpos = False
shit_stat_modres = False
shit_stat_misdns = False # missing density!

fastaseq = ""
pdbfile = ""

def download_pdb(pdb_id,dest_dir):
	#print "downloading %s" % ( pdb_id )
	url      = 'http://www.rcsb.org/pdb/files/%s.pdb.gz' % ( pdb_id.upper() )
	dest     = '%s/%s.pdb.gz' % ( os.path.abspath(dest_dir), ''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(8)) )
	wget_cmd = 'wget --quiet %s -O %s' % ( url, dest )
	#if remote_host:
		#wget_cmd = 'ssh %s %s' % ( remote_host, wget_cmd )

	system( wget_cmd )
	print wget_cmd
	if ( exists(dest) ):
		return dest
	else:
		print "Error: didn't download file!"

def check_and_print_pdb( count, residue_buffer, residue_letter ):
	global fastaseq
	global pdbfile

  ## Check that CA, N and C are present!def check_and_print_pdb( outid, residue_buffer )
	hasCA = False
	hasN = False
	hasC = False
	for line in residue_buffer:
		atomname = line[12:16]
		#Only add bb atoms if they have occupancy!
		occupancy=float(line[55:60])
		if atomname == " CA " and occupancy > 0.0: hasCA = True
		if atomname == " N  " and occupancy > 0.0:  hasN = True
		if atomname == " C  " and occupancy > 0.0:  hasC = True

  ## if all three backbone atoms are present withoccupancy proceed to print the residue
	if hasCA and hasN and hasC :
		for line in residue_buffer:
			## add linear residue count
			newnum = '%4d ' % count
			line_edit = line[0:22] + newnum + line[27:]
			## write the residue line
			pdbfile = pdbfile + line_edit

    ## finally print residue letter into fasta stream
		fastaseq = fastaseq + residue_letter

    ## count up residue number
		count = count + 1
		return True

	return False


files_to_unlink = []
assert( len(argv)>2)
pdbname = argv[1]
chainid = argv[2]

if (pdbname[-4:] != '.pdb' and pdbname[-8:] != '.pdb1.gz'):
	pdbname += '.pdb'

outfile = string.lower(pdbname[0:4]) + chainid + pdbname[4:]

nopdbout = 0
if argv.count('-nopdbout'):
	nopdbout = 1

removechain = 0
if argv.count('-nochain'):
	removechain = 1

ignorechain = 0
if argv.count('-ignorechain'):
	ignorechain = 1

netpdbname = '/net/wwpdb/' + pdbname[1:3] + '/' + pdbname
if not exists(netpdbname):
	netpdbname = pdbname

fixed_pdb = "/work/robetta/src/rosetta_server/fixed_pdbs/" + outfile
print  "Looking for: ", fixed_pdb
if os.path.isfile( fixed_pdb ):
		print "Found preoptimised or otherwise fixed PDB file. "
		netpdbname = fixed_pdb
else:
	print "File %s doesn't exist, downloading from internet." % ( netpdbname )
	netpdbname = download_pdb(pdbname[0:4],'.')
	files_to_unlink.append( netpdbname )

if netpdbname[-3:]=='.gz':
 	f = gzip.open( netpdbname , 'r' )
 	lines = f.readlines()
else:
 	lines = open(netpdbname,'r').readlines()


oldresnum = '   '
count = 1;
modifiedres = ''


residue_buffer = []
residue_letter = ''
residue_invalid = False

if chainid == '_':
	chainid = ' '

for i in range(len(lines)):
	line = lines[i]

	if len(line)>5 and line[:6]=='ENDMDL':break #Its an NMR model.

	if (chainid == line[21] or ignorechain):
		line_edit = line
		if line[0:3] == 'TER':
			continue
		elif (line[0:6] == 'HETATM'):
			ok = False

			## Is it a modified residue ?
			if modres.has_key( line[17:20] ):
			  ## if so replace it with its canonical equivalent !
				line_edit = 'ATOM  '+line[6:17]+modres[line[17:20]] +line[20:]
				modifiedres = modifiedres +  line[17:20]  + ',  '
				## dont count MSEs as modiied residues (cos they're so common and get_pdb deal with them previosuly)
				if line[17:20] != "MSE":
				  shit_stat_modres = True
				ok = True

			## other substitution (of atoms mainly)
			if (line[17:20]=='MSE'): #Selenomethionine
				if (line_edit[12:14] == 'SE'):
					line_edit = line_edit[0:12]+' S'+line_edit[14:]
				if len(line_edit)>75:
					if (line_edit[76:78] == 'SE'):
						line_edit = line_edit[0:76]+' S'+line_edit[78:]


			if not ok:
			  continue # skip this atom if we havnt found a conversion


		if line_edit[0:4] == 'ATOM': #or line_edit[0:6] == 'HETATM':

##			if line_edit[13:14]=='P': #Nucleic acid? Skip.
##				resnum = line_edit[23:26]
##				oldresnum = resnum
##				while (resnum == oldresnum):
##					print "HERE"
##					i += 1
##					line = lines[i]
##					resnum = line_edit[23:26]

			resnum = line_edit[22:27]

			insres = line[26]
			if insres != ' ': shit_stat_insres = True

			altpos = line[16]
			if altpos != ' ': shit_stat_altpos = True
			## Is thresidue_letter
			if not resnum == oldresnum:
				if residue_buffer != []:  ## is there a residue in the buffer ?
				  if not residue_invalid:
					if not check_and_print_pdb( count, residue_buffer, residue_letter ):
					  ## if unsuccessful
					  shit_stat_misdns = True
					else:
					  count = count + 1

				residue_buffer = []
				residue_letter = ''
				residue_invalid = False

				longname = line_edit[17:20]
				if longer_names.has_key(longname):
					residue_letter = longer_names[longname];
				else:
					residue_letter = 'X'
					residue_invalid = True

			oldresnum = resnum

			## What does this do ?
			if line_edit[16:17] == 'A':
				line_edit = line_edit[:16]+' '+line_edit[17:]

			if line_edit[16:17] != ' ':
				continue

			if removechain:
				line_edit = line_edit[0:21]+' '+line_edit[22:]


			residue_buffer.append( line_edit )


if not check_and_print_pdb( count, residue_buffer, residue_letter ):
	## if unsuccessful
	shit_stat_misdns = True
else:
	count = count + 1


flag_altpos = "---"
if shit_stat_altpos : flag_altpos = "ALT"
flag_insres = "---"
if shit_stat_insres : flag_insres = "INS"
flag_modres = "---"
if shit_stat_modres : flag_modres = "MOD"
flag_misdns = "---"
if shit_stat_misdns : flag_misdns = "DNS"

nres = len(fastaseq)

flag_successful = "OK"
if nres <= 0: flag_successful = "BAD"

print netpdbname, pdbname, chainid, "%5d"%nres, flag_altpos,  flag_insres,  flag_modres,  flag_misdns, flag_successful


if chainid == ' ':  chainid = '_'
if nres > 0:
	if ( nopdbout == 0 ):
		#outfile = string.lower( basename(outfile) )
		outfile = outfile.replace('.pdb1.gz','.pdb')
		outid = open( outfile, 'w')
		outid.write(pdbfile)
		outid.write("TER\n")
		outid.close()

	fastaid = stdout
	fastaid.write('>'+pdbname[0:4]+chainid+'\n');
	fastaid.write( fastaseq )
	fastaid.write('\n')

if len(files_to_unlink) > 0:
	for file in files_to_unlink:
		os.unlink(file)
