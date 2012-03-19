import Bio.PDB.Polypeptide
from rosettautil.util import fileutil


class pssm_map:
	def __init__(self,path,mode="score"):
		"""Parse a PSSM file from BLAST into a usable datastructure"""
		self.pssmmap = {}
		self.native_sequence = {}
		pssmfile = fileutil.universal_open(path,'r')

		
		pssmfile.readline()
		pssmfile.readline()
	
		header = pssmfile.readline()
		header = header.split()
		header = header[0:21]
		for line in pssmfile:
			#print line
			line = line.split()
			#self.native_sequence.append(
			if len(line) == 0:
				break
		
			res_num = int(line[0])
			res_id = line[1]
			self.native_sequence[res_num] = res_id
			line_map = {}
			
			if mode == "score":
				data = line[2:23]
			if mode == "percent":
				data = line[22:42]
				#print data
			for resname,score in zip(header,data):
				line_map[resname] = int(score)
			self.pssmmap[res_num] = line_map

		pssmfile.close()
	
	def get_score(self,seqpos,resname):
		"""get the score for a given sequence position and residue name"""
		if len(resname) == 3:
			resname =  Bio.PDB.Polypeptide.three_to_one(resname)
			linemap = self.pssmmap[seqpos]
			return linemap[resname]
			
		elif len(resname) == 1:
			linemap = self.pssmmap[seqpos]
			return linemap[resname]
		else:
			raise LookupError("this isn't a residue")
	def size(self):
		return len(self.pssmmap)

	def get_native_res(self,seqpos):
		return self.native_sequence[seqpos]
		
	
	def conserved(self,seqpos,resname):
		"""return true if the score of a given sequence position and residue name is greater than 0"""
		score = self.get_score(seqpos,resname)
		if(score >=0):
			return True
		else:
			return False
