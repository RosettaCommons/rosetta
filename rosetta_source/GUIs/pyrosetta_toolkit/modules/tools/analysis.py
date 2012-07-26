#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/modules/tools/analysis.py
## @brief  Functions for analysis in the toolkit
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)




from __future__ import division
from rosetta import *
import math
import random
import general_tools
import os
import os.path
import re
from pdbs import pdbTools
import loops as loop_tools
from pdbs import pdbTools
from collections import defaultdict
from operator import itemgetter
import heapq
from shutil import copy


def RetFAEnergy(p):
    FaScore=create_score_function_ws_patch('standard', 'score12')
    #start = int(start); end = int(end); print start; print end;
    #start = p.pdb_info().pdb2pose(chain, start); end = p.pdb_info().pdb2pose(chain,end)
    #emap=rosetta.core.scoring.TwoBodyEMapVector()
    print FaScore(p)
    print FaScore.show(p)

    #TotLoop=0
    #for i in range(start, end+1):
        #print "Total Energy"
        #print i
        #print emap[rosetta.core.scoring.total_score]
        
def RetPhiPsi(p):
    tot=p.total_residue()
    print "Residue - Phi - Psi"
    for i in range(1, tot+1):
        print repr(p.pdb_info().number(i))+":    "+repr(p.phi(i))+"  "+repr(p.psi(i))

def loopPhiPsiOmega(p, start, end):
    '''
    posStart and posEnd numbers....not the pdb ones...
    '''
    DicRama = dict()
    for i in range(start, end+1):
        rama="%.5f"%(p.phi(i))
        rama=rama+":"+"%.5f"%(p.psi(i))
        rama=rama+":"+"%.5f"%(p.omega(i))
        print rama
        DicRama[i] = rama
    return DicRama
    
class cdrCluster():
    '''
    Deprecated?
    '''
    
    L1=(24, 42, "L"); L2=(57, 72, "L"); L3=(107, 138, "L")
    H1=(24, 42, "H"); H2=(57, 69, "H"); H3=(107, 138, "H")
        
    def cdrClusterFile(self, p, num, outName = 0, outPath = 0):
        loops = (self.L1, self.L2, self.L3, self.H1, self.H2, self.H3)
        ldef = ("L1", "L2", "L3", "H1", "H2", "H3")
        i =0
        file = 0
        if outName !=0:
            file= 1
            full = outPath+"/"+outName + ".txt"
            print "Outputing to file:"
            print full
            if os.path.exists(full):
                fileout = open(full, 'a')
            else:
                fileout = open(full, 'w')
        print "Phi:Psi:Omega for use in CDR clustering..."
        for loop in loops:
            seq = tools.sequence.get_sequence(p, int(loop[0])+":"+int(loop[1])+":"+loop[2])
            print ldef[i] +"_"+repr(len(seq))+'\t'+seq
            start = p.pdb_info().pdb2pose(loop[2], loop[0])
            end = p.pdb_info().pdb2pose(loop[2], loop[1])
            DicRama = loopPhiPsiOmega(p, start, end)
            if file ==1:
                name = outName+"_"+repr(num)+"_"+ldef[i]+"_"+repr(len(seq))
                name = name+"\t"+seq
                fileout.write(name); fileout.write('\t')
                for x in DicRama:
                    out = DicRama[x]+ ":"
                    fileout.write(out)
                fileout.write('\n')
            i+=1
        if ((outName !=0 and outPath==0) or (outName ==0 and outPath!=0)):
            print "Please Specify both OutName and OutPath if trying to write to a file..."
        if file ==1:
            fileout.close()
            
    def cdrClusterFolder(self, dir, outName, outPath):
        fileList = os.listdir(dir)
        i = 1
        p = Pose();
        for file in fileList:
            if re.search(".pdb", file):
                file = dir +"/"+file
                print file
                pose_from_pdb(p, file)
                self.cdrClusterFile(p, i, outName, outPath)
                print i
                i+=1
        print "cdr Clustering Complete"
    
    def clusterAnalyze(self, file, pdb, pdbPaths=0, outPath=0, outName=0):
        print pdb
        fileIn = open(file, 'r')
        print outPath
        if outPath != "0/":
            out = outPath + "/" + outName
            fileOut=open(out, 'w')
        else:
            outPath = 0
        #if pdb!=0:
            #out = outPath + "/" + outName
        DicData = dict()
        DicTime = dict()
        DicCounts = defaultdict(int)
        DicType = dict()
        DicSwitchCount = defaultdict(int)
        DicSwitch = defaultdict(int)
        LisCounts = list()
        i = 1
        start = 0
        cdrDic=dict()
        pdbObject = pdbTools()
        for line in fileIn:
            sp = line.split(",")
            LisCounts.append(int(sp[3]))
            DicTime[i]=int(sp[3])
            if start != int(sp[3]):
                sw = repr(start)+":"+sp[3]
                print sw
                DicSwitch[sw]+=1
            #else:
                #DicSwitchCount[sp[3]]+=1
            start = int(sp[3])
            
            if int(pdb)==1:
                sp0 =(sp[0]).split("_")
                print sp0; fName = sp0[0]+"_"+sp0[1]+"_"+"10ps"+"_"+sp0[3]+".pdb"
                inFile = pdbPaths+"/"+fName
                newPath = outPath + "/Cluster_" +(sp[3])+"/"                   
                print newPath
                if os.path.exists(newPath)==0:
                    os.mkdir(newPath)
                cdrDic={'H1':self.H1, 'H2':self.H2, 'H3':self.H3, 'L1':self.L1, 'L2': self.L2,'L3':self.L3}
                cdr = cdrDic[sp0[4]]
                newFile = newPath + sp0[0]+"_"+sp0[4]+"_Clus"+sp[3]+"_"+sp0[3]+"-ps.pdb"
                pdbObject.outputCDR(cdr[0], cdr[1], cdr[2], inFile, newFile)

                #p=Pose()
                #make_pose_from_sequence(p, sp[1], "fa_standard")
                #print p
                #a, b, c = 0, 1, 2
                #for i in range(1, len(sp[1])):
                #   print sp[4+a] +" : " +sp[5+b] +" : "+sp[6+c]
                #    las = (sp[6+c]).replace('\r\n', '')
                #    p.set_phi(i, float(sp[4+a])); p.set_psi(i, float(sp[5+b])); p.set_omega(i, float(las))
                #    a, b, c = a+3, b+3, c+3
                #outN = newPath+"/"+(sp[0])+".pdb"
                #print p
                #dump_pdb(p, outN)
            i +=1    
        for i in LisCounts:
            DicCounts[i]+=1
        if outPath !=0:
            fileOut.write('Cluster:Total Frames\n')
        for clus in DicCounts:
            print "Cluster: " + repr(clus) + " Total: "+ repr(DicCounts[clus])
            if outPath !=0:
                out = repr(clus)+ "\t" + repr(DicCounts[clus]) + '\n'
                fileOut.write(out)
        if outPath != 0:
            fileOut.write('Avg Time in Cluster B/4 switch\n')
        for sw in DicSwitch:
            s = sw.split(":")
            DicSwitchCount[int(s[1])]= 0       
        for sw in DicSwitch:
            s = sw.split(":")
            DicSwitchCount[int(s[1])]=DicSwitch[sw]+DicSwitchCount[int(s[1])]
        DicTspent= dict()
        for clus in DicCounts:
            print "Cluster: " + repr(clus)
            tSpent = ((DicCounts[clus])*1.0)/((DicSwitchCount[clus])*1.0)
            print "Cluster: " + repr(clus) +" Avg b/4 switch: " +  repr(tSpent)
            DicTspent[clus] = tSpent
        sort = sorted(DicTspent, key=DicTspent.get, reverse=True)
        for clus in sort:
            if outPath != 0:
                out = repr(clus) + '\t' + "%.2f"%(DicTspent[clus]) + '\n'
                fileOut.write(out)
        if outPath !=0:
            fileOut.write('Switch:Observed\n')
        sort2 = sorted(DicSwitch.iteritems(), key = itemgetter(1), reverse=True)
        for sw in sort2:
            print "Switch: " + sw[0] + " Observed: " +repr(sw[1])
            if outPath != 0:
                out = sw[0] + '\t' + repr(sw[1]) + '\n'
                fileOut.write(out)
        if outPath !=0:
            fileOut.write('Time:Cluster\n')
        for time in DicTime:
            #print "Time: " + repr(time) + " Cluster: " + repr(DicTime[time])
            if outPath != 0:
                out = repr(time)+ "\t" + repr(DicTime[time]) + '\n'
                fileOut.write(out)
        fileIn.close(); fileOut.close()
        print "Analysis of Clustering Data Complete."
        return
    def clusterOutput(self, dirOut, outName):
        pass
    
    
    
def RetFAEnergyAll(p):
    for i in range(0, p.total_residue()):
        p.energies().show(i)
        
def rmsd(rmsdP, p, loops_as_strings, ca_only = False, all_atom=False):
    '''
    Returns RMSD for Full Protein, as well as any loops in loops_as_strings.
    '''
    
    print "\n"
    if ca_only:
        print "\nCA RMSD %.3f"%CA_rmsd(rmsdP, p)
    elif all_atom:
        print "\nAll Atom RMSD %.3f"%all_atom_rmsd(rmsdP, p)
    else:
        print "\nBB RMSD %.3f"%bb_rmsd(rmsdP, p)
        
    #if start !=0 and end!=0:
    #    start = p.pdb_info().pdb2pose(chain, int(start)); end = p.pdb_info().pdb2pose(chain, int(end))
    #    cut = start+((end-start)/2); loo=Loop(start,end, cut)
    #    loops = Loops()
    #    loops.add_loop(loo)
        
    #    lrms = loop_rmsd(p, rmsdP, loops, False)
    #    print "Loop RMSD:" + str(lrms)
    
    if all_atom:bb_only = False
    else: bb_only=True
    if loops_as_strings:
        all_rosetta_loops = Loops()
        for loop_string in loops_as_strings:
            rosetta_loop = loop_tools.return_rosetta_Loop(p, loop_string)
            single_loops = Loops()
            
            single_loops.add_loop(rosetta_loop)
            all_rosetta_loops.add_loop(rosetta_loop)
            
            lrms = loop_rmsd(p, rmsdP, single_loops, ca_only, bb_only)
            print "\n"+loop_string+" RMSD %.3f"%lrms
        lrms = loop_rmsd(p, rmsdP, all_rosetta_loops, ca_only, bb_only)
        print "\nALL Loop RMSD:%.3f"%lrms
    return

def readFASC(fileName):
    File = open(fileName, 'r')
    fascData = dict()
    for line in File:
        lineSplit = line.split()
        if lineSplit[0]!="pdb":
            x = 3
            for i in range(3, (len(lineSplit)/2)+2):
                if not fascData.has_key(lineSplit[1]):
                    fascData[lineSplit[1]]=dict()
                fascData[lineSplit[1]][lineSplit[x-1]]=lineSplit[x]
                x = x+2
    return fascData
        
    
class Decoys():
    '''
    Depends on Fasc.  This is important.  JD2 only!!!  Needs to be rewritten if we are going to include it.
    '''
    def findLowDecoys(self, inDir, inName, settings, loops_as_strings, outDir=0):
        #Settings - [check_button_RMSD, check_button_LRMSD, check_button_Energy, check_button_Move, Percent, score]
        fileList = os.listdir(inDir)
        filenames = []
        topFileDic = dict()
        x = 0
        for file in fileList:
            file = inDir +"/"+file
            if re.search(inName, file) and re.search(".fasc", file):
                fasc = file
                x = 1
            elif re.search(inName, file) and re.search(".pdb", file):
                filenames.append(file)
        if x == 1:
            fascData = readFASC(fasc)
            total = len(filenames)
        else:
            print "Fasc File Not found..."
            return

        top = int(round(total*((int(settings[4].get()))/100)))
        if int(settings[4].get())==100:
            top = total
        print "TOP POSES: "+ repr(top)
        #Only Energies
        Escore=dict()
        if int(settings[3].get()) == 1:
            newPath = outDir + "/" +settings[5].get()
            if os.path.exists(newPath)==0:
                os.mkdir(newPath)
            file = newPath+"/"+ "pdbRecord.txt"
            pdbRecord = open(file, 'w')
            pdbRecord.write("Scores for: "+ settings[5].get()+" of "+ inName)
        topFiles = []
        if int(settings[2].get())==1:
            print "Returning Top Energies..."
            for key in fascData:
                Escore[key]=fascData[key][settings[5].get()+":"]
            topScores = heapq.nsmallest(total, Escore.iteritems(), itemgetter(1))
            print "Scores for: "+ settings[5].get()
            ave = 0
            for i in range(1, top+1):
                #Fixes the ordering of what heapq considers smallest (It's Dumb)
                if float(topScores[len(topScores)-1][1]) >= 0:
                    smallest = i-1
                else:
                    smallest = (len(topScores)-i)
		sp = topScores[smallest][0].split("/")
                name = sp[len(sp)-1]
                
                print name + " - "+topScores[smallest][1]
                topFiles.append(topScores[smallest][0])
                topFileDic[name]=topScores[smallest][1]
                ave = ave +float(topScores[smallest][1])
                if int(settings[3].get())==1:
                    scorefix = topScores[smallest][0].split("/")
                    newFile = inDir + "/" + scorefix[len(scorefix)-1]
                    copy(newFile, newPath)
		    
                    pdbRecord.write(name + "\t"+topScores[smallest][1]+ "\n")
                    
            #Global Helpful Information:
            if float(topScores[len(topScores)-1][1]) >= 0:
                max = float(topScores[len(topScores)-1][1])
                min = float(topScores[0][1])
            else:
                max = float(topScores[0][1])
                min = float(topScores[len(topScores)-1][1])
            ave = ave/(top); print "Average: " + repr(ave)
            print "Min: " + repr(min) +"Max: " +repr(max)
            
                
        
        #Additions - print Standard Deviation; Fix Negative/Positive crap
        #Additions - Select multiple, and get top from all!
        #Only RmsD
        elif int(settings[0].get())==1:
            for key in fascData:
                Rscore[key]=fascData[key]["rmsd:"]
            topScores = heapq.nsmallest(top, Escore.iteritems(), itemgetter(1))
            print topScores
        #Only Loop RMSD
        elif int(settings[1].get())==1:
            pass
        else:
            print "Please Specify what you would like to do"
            return
        if int(settings[3].get())==1:
            pdbRecord.close()
        #Then, we figure out how to do combos...
        #topFiles is a list of the top files, in order, while topFileDic has file:score
        return topFiles, topFileDic
    
    
    
    
    
    
    
    
    def seqRecovery(sefl, loops_as_strings, settings, inName, inDir, fileout):
        '''
        This is meant to be independent of the fasc/sc file.  Returns a sequence dictionary.
        Independant of score.  (Hopefully) DOES NOT WORK
        '''

        
        if re.search(inName, file) and re.search(".pdb", file):
                filenames.append(file)

    
    def AnalyzeSeq(self, loops_as_strings, settings, inName, inDir, fileout, PDBcompare, seqOut):
        '''
        This starts the analysis of decoy sequences and analyzes them, printing out info.
        inName is the name specified to look for decoys.
        '''
        topFiles, topFileDic = self.findLowDecoys(inDir, inName, settings, loops_as_strings)
        p = Pose()
        #This makes sure it works on any computer you do the analysis that is different from the generation of decoys
        fixedFiles = []
        for files in topFiles:
            files = files.split("/")
            fixedFile = inDir + "/" +files[(len(files)-1)]
            fixedFiles.append(fixedFile)
        topFiles = fixedFiles
        #Loads one file into rosetta to get the numbering of residues easily.
        #Each other file is not loaded into rosetta.
        print topFiles[0]
        pose_from_pdb(p, topFiles[0])
        PDBComSplit = PDBcompare.split("/")
        topFiles.append(PDBcompare)
        topFileDic[PDBComSplit[len(PDBComSplit)-1]] = "NA"
        seqDic = self.loadSeq(topFiles)
        result = self.compareSeq(p, seqDic, loops_as_strings, inName, fileout, topFiles, topFileDic, PDBcompare, seqOut)
        
        if result =="error":
            print "Error"
            return
        else:
            print "Analysis Complete...."
            return
        
    def loadSeq(self, filenames):
        '''
        Input is List of filenames, each sequence is loaded into a Sequence Dictionary and returned.
        '''
        seqDic = dict()
        p = Pose()
        for file in filenames:
            seqDic[file] = pdbTools().getSeq(file)
            
        return seqDic
    
    def compareSeq(self, p, seqDic, loops_as_strings, inName, fileoutReal, topFiles, topFileDic, PDBCompare, seqOut, native = 0, settings = 0):
        '''
        Input is Dictionary of Sequences, the function compares these and outputs stats for each residue.
        '''
        pCompare = Pose()
        pose_from_pdb(pCompare, PDBCompare)
        pwd = os.getcwd()
        fileout = pwd + "/Seqtemp"
        newLoops = []
        if len(loops_as_strings)!=0:
            for loo in loops_as_strings:
                looSp = loo.split(":")
                if loo[0] == ":" and loo[1]==":":
                    print "Comparing chain: "+ loo[2]
                    print "Currently not implemented..."
                    return "error"
                elif looSp[0]=="":
                    print "Comparing N-terminus of chain "+looSp[2] + " From" + looSp[1] + " to N-Term end of chain..."
                    print "Currently not implemented..."
                    return "error"
                elif looSp[1]=="":
                    print "Comparing C-terminus of chain "+looSp[2]+ " From" + looSp[0] + " to C-Term end of chain "
                    print "Currently not implemented..."
                    return "error"
                else:
                    print "Comparing Loops"
                    start = int(looSp[0]); end = int(looSp[1]); chain = looSp[2]
                    start = p.pdb_info().pdb2pose(chain, start); end = p.pdb_info().pdb2pose(chain, end)
                    #NewLoops are loops in rosetta numbering.
                    newLoops.append((repr(start)+":"+repr(end)))
                    
        else:
            print "Comparing Whole protein"
            start = 1; end = p.total_residue()
            newLoops.append((repr(start)+":"+repr(end)))
        resDic = dict()
        #Seq Dictionary Data Structure:
            #3-D Dictionary
            #"start-end:Sequence" -->Position-->AA = AA's seen at each position.
        realSeqList = []

        File = open(fileout, 'w')   ;#-1 accounts for native pose in counts...
        for c in range(0, (len(topFiles)-1)):
            fileIn = topFiles[c]
            seq = seqDic[fileIn]
            seqList = ""
            for loops in newLoops:
                f = fileIn.split("/")
                fLast = len(f)-1;
                f = f[fLast]
                looSp = loops.split(":")
                start = int(looSp[0]); end = int(looSp[1])
                LooSeq = seq[(start-1):(end)]
                print "LooSeq" + LooSeq
                seqList = seqList+loops+" : "+LooSeq +"\t"
                resDicKey =repr(start)+"-"+repr(end)
                if not resDic.has_key(resDicKey):
                    resDic[resDicKey]=dict()
                    for i in range(0, len(LooSeq)):
                        if not resDic[resDicKey].has_key(i):
                            resDic[resDicKey][i] = dict()
                            
                        if not resDic[resDicKey][i].has_key(LooSeq[i]):
                            resDic[resDicKey][i][LooSeq[i]]=1
                        else:
                            resDic[resDicKey][i][LooSeq[i]]+=1
                else:
                    for i in range(0, len(LooSeq)):
                        if resDic[resDicKey][i].has_key(LooSeq[i]):
                            resDic[resDicKey][i][LooSeq[i]]+=1
                        else:
                            resDic[resDicKey][i][LooSeq[i]]=1
	    fileInSp = fileIn.split("/")
	    f = fileInSp[len(fileInSp)-1]
            if seqOut =="1":
                print f + "\t"+seqList
                File.write(f+"\t"+seqList+topFileDic[f]+"\n")
            #if fileout:
                #File.write(fileIn + "\t"+seqList+"\n")
            
        #Heres the Analysis/Printout of the dictionary of mutations in each sequence in each loop.
        AA = ("GPAVLIKRQNEDFYWHMCST")
        tot = len(seqDic)
        #File.write("\t"+ inName + ":\n\n")
        #File.write("Native Sequence and Energy)
        #Writes to a temp file, all the data for each position and each loop specified
        for key in resDic:
            File.write(key+":\n")

            loo = key.split("-");
            start = p.pdb_info().pose2pdb(int(loo[0])); end = p.pdb_info().pose2pdb(int(loo[1]))
                
            StNum = start.split(); EnNum = end.split()
	    File.write("::")
            #May have to Fix for missing residue numbers....ugh...
            for i in range(int(StNum[0]), int(EnNum[0])+1):
                cur = repr(i) + StNum[1]
		cur = cur.rjust(5)
                File.write(cur)
	    seqLen = int(EnNum[0])-int(StNum[0])
            File.write("\n")
            File.write(":: ")
            x = pdbTools()
            olc = x.OneLetterCode
            #Writes Original Sequence for easy comparison
            for i in range(int(StNum[0]), int(EnNum[0])+1):
                res = p.pdb_info().pdb2pose(StNum[1], i)
                cur = pCompare.residue(res).name()
                cur = olc[cur]
                cur = cur.ljust(5)
                File.write(cur)
            
            File.write("\n")
            #Compares each position to the native, and prints
            File.write(":: ")
            
            v = 0
            for i in range(int(StNum[0]), int(EnNum[0])+1):
                res = p.pdb_info().pdb2pose(StNum[1], i)
                cur = pCompare.residue(res).name()
                cur = olc[cur]
                if resDic[key][v].has_key(cur):
                    count = (resDic[key][v][cur])
                    per = count/(tot-1)
                    if per == 1.00:
                        File.write("-----")
                    else:
                        per = "%.2f"%(per)+" "
                        File.write(per)
                else:
                    File.write("MUT  ")
                v+=1 
            File.write("\n")
            File.write(":: ")
            for i in range(int(StNum[0]), int(EnNum[0])+1):
                File.write(":::::")
            File.write("\n")
            #Creates the List AA found for each position and their numbers
            for i in range(0, len(AA)):
                File.write(AA[i]+"  ")
                for x in range(0, seqLen+1):
                    if resDic[key][x].has_key(AA[i]):
                        count = (resDic[key][x][AA[i]])
                        per = count/(tot-1)
                        per = "%.2f"%(per)+" "
                        File.write(per)
                    else:
                        File.write("0    ")
                File.write("\n")
            File.write("Total Poses:" +repr(tot-1)+"\n")
        File.close()
        
        #Here, we load the file into memory, and print it out line by line.  If a filename was specified, we save it to the correct place.
        lines = pdbTools().loadPDB(fileout)
        for line in lines:
            print line
        if fileoutReal!=0:
            FILEOutReal = open(fileoutReal, 'w')
            for line in lines:
                FILEOutReal.write(line)
            FILEOutReal.close()
        os.remove(fileout)
        
        #Now, we print out specific, sequence recovery info.
        #print "Sequence Recovery:"
        return

        
class rotamers():
    '''
    Does this work?
    '''
    
    def retEn(self, p, res, chain=0, type = fa_dun):
        #("fa_atr", "fa_rep", "fa_sol", "fa_intra_rep", "pro_close", "fa_pair", "hbond_sr_bb", "hbond_lr_bb", "hbond_bb_sc", "hbond_sc", "dslf_ss_dst", "dslf_cs_ang", "dslf_ca_dih", "fa_dun", "p_aa_pp")
        #Wish there was a better way for this...
        if type =="fa_atr":
            type = fa_atr
        elif type == "fa_sol":
            type = fa_sol
        elif type == "fa_intra_rep":
            type = fa_intra_rep
        elif type == "fa_pair":
            type = fa_pair
        elif type == "hbond_sr_bb":
            type = hbond_sr_bb
        elif type == "hbond_bb_sc":
            type = hbond_bb_sc
        elif type == "hbond_lr_bb":
            type =hbond_lr_bb
        elif type == "hbond_sc":
            type = hbond_sc
        elif type=="dslf_ss_dst":
            type =dslf_ss_dst
        elif type == "dslf_cs_ang":
            type = dslf_cs_ang
        elif type == "dslf_ca_dih":
            type = dslf_ca_dih
        elif type == "fa_dun":
            type=fa_dun
        elif type == "p_aa_pp":
            type=p_aa_pp
        elif type == fa_dun:
            ":)"
        else:
            print "We have a problem..."
            return
        if chain !=0:
            res = p.pdb_info().pdb2pose(chain, int(res))
        score = create_score_function_ws_patch('standard', 'score12')
        emap = core.scoring.EMapVector()
        score.eval_ci_1b(p.residue(res), p, emap)
        e = emap[type]
        return e
    
    
    
    
    def retProb(self, p, res, chain=0):
        e = self.retEn(p, res, chain)
        p = math.exp(-e)
        return p
    

class interface():
    '''
    Does this work?
    '''
    def interface_ddg(self, pose, resid, mut_res):
        '''
        Author: Sid Chaudhury
        '''
        
        dock_jump = 1
        wt_p = Pose()
        mut_p = Pose()

        wt_p.assign(pose)
        mut_p.assign(pose)
        mutate_residue(mut_p, resid, mut_res)

        scorefxn = ScoreFunction()
        scorefxn.set_weight(fa_atr, 0.44)
        scorefxn.set_weight(fa_rep, 0.07)
        scorefxn.set_weight(fa_sol, 1.0)
        scorefxn.set_weight(hbond_bb_sc, 0.5)
        scorefxn.set_weight(hbond_sc, 1.0)

        dock_prot = DockingProtocol(dock_jump)
        dock_prot.set_highres_scorefxn(scorefxn)

        #setting up packer options
        tf = standard_task_factory()
        tf.push_back(RestrictToRepacking())
        prevent_repacking = PreventRepacking()

        for i in range(1, pose.total_residue()+1):
            rsd1 = pose.residue(resid)
            rsd2 = pose.residue(i)
            if (rsd1.nbr_atom_xyz().distance_squared(rsd2.nbr_atom_xyz())>64.0):
                prevent_repacking.include_residue(i)	

        tf.push_back(prevent_repacking)

        #setup packer	
        pack_mover = PackRotamersMover(scorefxn)
        pack_mover.task_factory(tf)

        #applying side chain packer to wt and mutant
        pack_mover.apply(wt_p)
        pack_mover.apply(mut_p)

        wt_score = dock_prot.calc_interaction_energy(wt_p)
        mut_score = dock_prot.calc_interaction_energy(mut_p)
        ddg = mut_score - wt_score

        pdbname = p.pdb_info().chain(resid)+str(p.pdb_info().number(resid))+wt_p.sequence()[resid-1]+"to"
        #dump_pdb(wt_p,pdbname+ wt_p.sequence()[resid-1]+".pdb")
        #dump_pdb(mut_p,pdbname+ mut_p.sequence()[resid-1]+".pdb")
        return ddg
    
    def analyzeDDG(self, p):
        '''
        Author: Sid Chaudhury
        '''
        
        starting_p = Pose()

        #parameters
        dock_jump = 1
        interface_dist = 6.0 #angstroms
        
        #This part is not working, and must be changed to pyrosetta 2.0 stuff...
        DockingProtocol().set_autofoldtree(True)
    
        scorefxn = create_score_function('standard')
        scorefxn(p) #needed for proper Interface calculation
        
        interface = Interface(dock_jump)
        interface.distance(8.0)

        interface.calculate(p)

        #store starting structure
        starting_p.assign(p)

        ddg_scorefxn = ScoreFunction()
        ddg_scorefxn.set_weight(fa_atr, 0.44)
        ddg_scorefxn.set_weight(fa_rep, 0.07)
        ddg_scorefxn.set_weight(fa_sol, 1.0)
        ddg_scorefxn.set_weight(hbond_bb_sc, 0.5)
        ddg_scorefxn.set_weight(hbond_sc, 1.0)

        for i in range(1, p.total_residue()+1):
            if (interface.is_interface(i) == True):
                p.assign(starting_p)
                ddg = self.interface_ddg(p, i, 'A')
                #f = open("ala_scan_output.txt", 'a')
                mutname = p.pdb_info().chain(i)+str(p.pdb_info().number(i))+p.sequence()[i-1]+"to"+'A'+".pdb" 
                #f.write( mutname + "  ddG: " + str(round(ddg, 3)) +'\n')]
                print mutname + "  ddG: " + str(round(ddg, 3)) +'\n'

                
