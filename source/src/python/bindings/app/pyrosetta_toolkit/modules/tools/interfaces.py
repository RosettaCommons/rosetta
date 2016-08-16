#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   /GUIs/pyrosetta_toolkit/modules/tools/interfaces.py
## @brief  original interface functions, before I found them in Rosetta.  Still may be usefull...
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#Python Imports
from time import clock
from shutil import rmtree
import os

#Toolkit Imports
from app.pyrosetta_toolkit.modules.PythonPDB import PythonPDB
import app.pyrosetta_toolkit.modules.tools.loops
import general_tools

pwd = os.getcwd()

class around():
	"""
	Basic vacinity and interface class...
	"""

	def getVicinity(self, p, LisLoop, infile, cutoff = 5.0):
		"""
		Get vaccinity around each loop/residue/chain specified in LisLoop.
		Returns a VacDic: ([res:chain]=atomic contact #)
		Needs to be rewritten using Evan's way of finding neighbors.
		"""
		pdb_map = dict()
		pdb_map = PythonPDB(infile).get_pdb_map()
		#Get individual residues in LisLoop
		resList = [] #List of residues to compute vacinity.  Res:chain
		resList = loops.loopArea(p, LisLoop, 0)
		pdbResCoord = dict() #This is a dictionary with [res:chain] for each residue we need to find the vacinity for.  Dictionary is a list (x, y, z)
		for i in pdb_map:
			for res in resList:
				res = res.split(":")
				resNum = res[0]; resChain = res[1]
				if pdb_map[i]["resNum"]==resNum and pdb_map[i]["chain"]==resChain:
					#This is where we do get the coordinates of each residue.
					pdbResCoord[resNum+":"+resChain] = (float(pdb_map[i]["x"]), float(pdb_map[i]["y"]), float(pdb_map[i]["z"]))

		#Now we measure every atom in pdb_map.  Adding it if it is within x of the residues in the list.
		vacDic = dict() #This is a dictionary with [contact]=# seen (in case we need to use it later...)
		for res in pdbResCoord:
			for i in pdb_map:
				xyz2 = (float(pdb_map[i]["x"]), float(pdb_map[i]["y"]), float(pdb_map[i]["z"]))
				res_chain = pdb_map[i]["resNum"]+":"+pdb_map[i]["chain"]
				d = tools.general_tools.getDistGen(pdbResCoord[res], xyz2)
				if d <= cutoff:
					#print "Cutoff:"+repr(cutoff); print d
					#check_button_ck to make sure this is not BS due to python craptastical numbering.

					if vacDic.has_key(res_chain):
						vacDic[res_chain]+=1 #check_button_ck that this works
					else:
						vacDic[res_chain]=1
		for res in pdbResCoord:
			vacDic.pop(res)
		for key in vacDic:
			print "Residues in Vacinity:"
			print "Residue: " +key
			#print "AtomContacts: "+ repr(vacDic[key])
		print "Total Residues in Vicinity: " +repr(len(vacDic))
		return vacDic

	def getDimerInterface(self, p, chainA, chainB, infile, interDic, cutoff = 5.0):
		"""
		Gets residues common to your interface of chainA and chainB
		Returns a VacDic: ([res:chain]=atomic contact #)
		"""

		#Get pdb_map without Whitespace
		pdb_map = dict()
		pdb_map = PythonPDB(infile).get_pdb_map()
		#Get individual residues in LisLoop
		resList = [] #List of residues to compute vacinity.  Res:chain
		chain = "::"+chainA
		chainList = []; chainList.append(chain)
		print chain
		resList = loops.loopArea(p, chainList, 0)
		pdbResCoord = dict() #This is a dictionary with [res:chain] for each residue we need to find the vacinity for.  Dictionary is a list (x, y, z)
		for i in pdb_map:
			for res in resList:
				res = res.split(":")
				resNum = res[0]; resChain = res[1]
				if pdb_map[i]["resNum"]==resNum and pdb_map[i]["chain"]==resChain:
					#This is where we do get the coordinates of each residue.
					pdbResCoord[resNum+":"+resChain] = (float(pdb_map[i]["x"]), float(pdb_map[i]["y"]), float(pdb_map[i]["z"]))

		#Now we measure every atom in pdb_map.  Adding it if it is within x of the residues in the list, and part of the interface between only the two.
		print "Calculating Distances...."
		for res in pdbResCoord:
			for i in pdb_map:
				xyz2 = (float(pdb_map[i]["x"]), float(pdb_map[i]["y"]), float(pdb_map[i]["z"]))
				res_chain = pdb_map[i]["resNum"]+":"+pdb_map[i]["chain"]
				#res_chainB= res+":"+chainA
				d = tools.general_tools.getDistGen(pdbResCoord[res], xyz2)
				if (d <= cutoff) and pdb_map[i]["chain"]==chainB:
					#print "Cutoff:"+repr(cutoff); print d
					#check_button_ck to make sure this is not BS due to python craptastical numbering.

					if interDic.has_key(res_chain):
						interDic[res_chain]+=1 #check_button_ck that this works
					else:
						interDic[res_chain]=1

					if interDic.has_key(res):
						interDic[res]+=1
					else:
						interDic[res]=1
		#for res in pdbResCoord:
			#interDic.pop(res)
		for key in interDic:
			print "Residues in interface between chain"+chainA+"and chain "+chainB+":"
			print "Residue: " +key
			#print "AtomContacts: "+ repr(vacDic[key])
		print "Total Residues in half interface: " +repr(len(interDic))
		return interDic

	def getInterface(self, p, interface, cutoff = 5.0):
		interfaceDic = dict(); #Dictionary: [::chain] = (res:chain, res:chain, etc.).
				      #Represents interfaces between chain and other chains.  Can be parsed for more info later down the line...
		#An interface is defined like this:  'A:CD' or A:B:EF'

		interfaceSP = interface.split(":")
		if len(interfaceSP)==1:
			print "Please specify an actual interface..."
			return
		#Dumps PDB so we can read it and find out what is in the vaccinity through good old python.
		t = clock()
		tempdir = pwd + "/temp2/" + "_"+repr(t)
		os.mkdir(tempdir)
		temp = tempdir + "/temp.pdb"; p.dump_pdb(temp)

		#Here, we will specify all the combinations we need to make in order to return the interface.
		#Then, we parse the data as needed?
		comboInterface = []
		#This is basically the card sorting algorithm...
		tempinterface = interfaceSP
		interfaceDic = dict(); #This is the final interface dictionary used by everything else.  I know.  It doesn't look too good.
		interDic = dict(); #This is the dictionary with keys as res:chain pairs = counts.  Thats it.
		while len(tempinterface) >1:
			for chains in interfaceSP:
				tempinterface.remove(chains)
				for chainA in chains:
					for x in tempinterface:
						for chainB in x:
							print chainA+":"+chainB
							interDic = self.getDimerInterface(p, chainA, chainB, temp, interDic, cutoff)

		allchains = ""
		for chains in interfaceSP:
			allchains = chains+allchains
		LisLoop = []
		for chain in allchains:
			LisLoop.append("::"+chain)
		for item in LisLoop:
			interfaceDic[item] = []
			tempLoop = []
			tempLoop.append(item); #Spits the LisLoop into individual chains.
			itemSP = item.split(":")
			chain = itemSP[2]
			for key in interDic:
				keySP = key.split(":")
				if keySP[1] != chain:
					interfaceDic[item].append(key)

		#Removes temp files
		rmtree(tempdir)
		return interfaceDic; #Interface dictionary is a simple dictionary with each chain as key, and the corresponding residues that that chain makes as interface as value (not including itself)

class Antibodies:
	def __init__(self):
		self.cdrDic = dict()
		self.cdrDic["L1"]=(24, 42)
		self.cdrDic["L2"]=(57, 72)
		self.cdrDic["L3"]=(107, 138)

		self.cdrDic["H1"]=(24, 42)
		self.cdrDic["H2"]=(57, 69)
		self.cdrDic["H3"]=(107, 138)

	def getContactCDRs(self, p, cutoff):
		"""
		This function gets the residues making contact - any contact - with other residues not in H or L.
		It returns the residues within the CDR's for use in other things.
		Must be renumbered antibodies.
		Should be expanded to get more data and help in analysis of individual CDR residue contacts with other CDRs and the framework/dimer interface.
		"""
		#First - Find chain(s) making contact with the CDR's.
		#Make LisLoop, measure with vaccinity, then test whether the vaccinity has any non H or L chains.
		contactResidues = dict()
		for key in self.cdrDic:
			print key
			start = p.pdb_info().pdb2pose(key[0], self.cdrDic[key][0])
			end   = p.pdb_info().pdb2pose(key[0], self.cdrDic[key][1])
			for i in range(start, end+1):
				LisLoop = []
				x = (p.pdb_info().pose2pdb(i)).split()
				LisLoop.append(x[0]+":"+x[0]+":"+x[1])
				residue = x[0]+":"+x[0]+":"+x[1]
				vacDic = tools.input.load_vicinity(p, LisLoop, cutoff)
				for key in vacDic:
					keySP = key.split(":")
					if (keySP[1]!="H") and (keySP[1]!="L"):
						if contactResidues.has_key(residue):
							contactResidues[residue]+=1
						else:
							contactResidues[residue]=1
		print "CDR Residues in Contact with Antigen: "+repr(contactResidues)
		return contactResidues
