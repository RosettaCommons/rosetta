import os
import re
import tkFileDialog

class pdbTools:
    def __init__(self):
        self.All =("Alanine:ALA:A", "Arginine:ARG:R", "Asparagine:ASN:N", "Aspartate:ASP:D", "Cysteine:CYS:C", "Glutamate:GLU:E", "Glutamine:GLN:Q", "Glycine:GLY:G", "Histidine:HIS:H", "Leucine:LEU:L", "Isoleucine:ILE:I", \
            "Lysine:LYS:K", "Methionine:MET:M", "Phenylalanine:PHE:F", "Proline:PRO:P", "Serine:SER:S", "Threonine:THR:T", "Tryptophan:TRP:W", "Tyrosine:TYR:Y", "Valine:VAL:V")
        self.OneLetterCode = dict()
        for  resFull in self.All:
            res = resFull.split(":")
            self.OneLetterCode[res[1]]=res[2]
        self.OneLetterCode["HIS_D"]="H"    
        
    def loadPDB(self, filename):
        '''
        Loads PDB/file into an array of lines
        '''
        
        FILE = open(filename)
        lines = FILE.readlines()
        FILE.close()
        return lines
    def cleanPDB(self, pdbDic):
        '''
        Removes HSD, Waters: Tries to fix atom name inconsistencies.
        Do not USE. Currently, untested for everything, and may not always work.
        '''
        if not pdbDic:
            print "Please Load PDB into Python to work on..."
            return
        #lines = self.loadPDB(filename)
        #pdbDic = dict()
        #pdbDic = self.pdbData(lines, pdbDic)
        His = 0; Water=0; ions=0#Keeps track of if a HSD or TIP residue is found (So that its not constantly printing)
        waters = []; #List of keys that have waters
        
        for key in pdbDic:
            if pdbDic[key]["ResName"] == "HSD ":
                His = 1
                pdbDic[key]["ResName"]="HIS "
            elif pdbDic[key]["ResName"]=="TIP3":
                Water=1
                waters.append(key)
            elif pdbDic[key]["chain"]=="I":
                ions = 1
                waters.append(key)
            #Changes Atom Names to correct ones: (Not sure why replace did not work)
            #alias = {}
            if pdbDic[key]["ResName"]== "SER ":
                dic = {"  HG1":"  HG "}

            elif pdbDic[key]["ResName"]=="ILE ":
                dic = {"  CD ":"  CD1"}
                pdbDic = self.pdbAlias(pdbDic, key, dic)    
            elif pdbDic[key]["ResName"]=="LEU ":
                dic = {"  OT1":"  O  ", "  OT2":"  OXT"}
                pdbDic = self.pdbAlias(pdbDic, key, dic)
            elif pdbDic[key]["ResName"]=="VAL ":
                dic = {"  OT1":"  O  ", "  OT2":"  OXT"}
                pdbDic = self.pdbAlias(pdbDic, key, dic)
            elif pdbDic[key]["ResName"]=="LYS ":
                dic = {"  HZ1":"  1HZ", "  HZ2":"  2HZ", "  HZ3":"  3HZ"}
                pdbDic = self.pdbAlias(pdbDic, key, dic)    
            elif pdbDic[key]["ResName"]=="ARG ":
                dic = {" HH11":" 1HH1", " HH12":" 2HH1", " HH21":" 1HH2", " HH22":" 2HH2"}
                pdbDic = self.pdbAlias(pdbDic, key, dic)
            elif pdbDic[key]["ResName"]=="ASN ":
                dic = {"HD21":"1HD2", "HD22":"2HD2"}
                pdbDic = self.pdbAlias(pdbDic, key, dic)
            elif pdbDic[key]["ResName"]=="PRO ":
                dic = {"  OT1":"  O  ", "  OT2":"  OXT", "  HD1":"  1HD", "  HD2":"  2HD", "  HB1":"  1HB", "  HG1":"  1HG", "  HG2":"  2HG"}
                pdbDic = self.pdbAlias(pdbDic, key, dic)    

        #Removes Waters+Ions
        for key in waters:
            pdbDic.pop(key)
        #Outputs what was found:
        if His ==1:
            print "HSD found...changed"
        if Water ==1:
            print "Water found...changed"
        if ions ==1:
            print "Ions found...changed"
        print "PDB file FIXED!!!"
        
 
        return pdbDic
    def pdbAlias(self, pdbDic, key, Dickey):
        '''
        Replaces atomNames with ones Rosetta is happy with.
        Dickey is all of the residues that need to be replaced start:end
        '''
        for start in Dickey:
            if pdbDic[key]["atomName"]==start:
                print pdbDic[key]["ResName"]+":"+pdbDic[key]["atomName"]+":"+ Dickey[start]
                pdbDic[key]["atomName"]=Dickey[start]
        return pdbDic
    '''
    def redoPDB(self, p, start, end, chain, inFile):
        pwd = os.getcwd(); temp=pwd +"/"+"tempDel"; temp2=pwd+"/"+"tempDel2.pdb"
        FILE2 = open(temp2, 'w')
        print p
        tools.output.dumpPDB(p, temp, inFile)
        temp = temp +"_1.pdb"
        lines = self.loadPDB(temp)
        i=1
        j=1
        pr=1
        for line in lines:
            if line.rfind("ATOM")!=1:
                for r in range(int(start), int(end)+1):
                    string = " "+repr(r)
                    if (line.rfind(string, 22, 26) !=-1) and (line.rfind(chain, 21) !=-1):
                        pr=0
                        print "Deletion Found"
                        
                        print line
            if pr ==1:
                newi = (str(i)).rjust(5)
                newLine = line[:6] + newi +line[11:]
                FILE2.write(newLine)
                i+=1
            pr=1
            j +=1
        p2 = Pose()
        pose_from_pdb(p2, temp2)
        print p2
        FILE2.close()
        os.remove(temp); os.remove(temp2)
        return p2
    '''
    def outputCDR(self, start, end, chain, inFile, outFile):
        #pwd = os.getcwd(); temp=pwd +"/"+"tempDel"; temp2=pwd+"/"+"tempDel2.pdb"
        FILE2 = open(outFile, 'w')
        #tools.output.dumpPDB(p, temp, inFile)
        #temp = temp +"_1.pdb"
        lines = self.loadPDB(inFile)
        i=1
        j=1
        pr=1
        for line in lines:
            if line.rfind("ATOM")!=-1:
                for r in range(int(start), int(end)+1):
                    string = " "+repr(r)
                    if (line.rfind(string, 22, 26) !=-1) and (line.rfind(chain, 21) !=-1):
                        pr=0
                        print "CDR Found"
                        FILE2.write(line)
            #if pr ==1:
                #newi = (str(i)).rjust(5)
                #newLine = line[:6] + newi +line[11:]
                #FILE2.write(newLine)
                #i+=1
            #pr=1
            #j +=1
        #p2 = Pose()
        #pose_from_pdb(p2, temp2)
        #print p2
        FILE2.close()
        #os.remove(temp); os.remove(temp2)
        #return p2
        
    def pdbData(self, lines, pdbDic, clean = 1, chains=0):
        '''
        Reads an array of lines and parses the contents of each line in accordance with pdb structure.
        Returns a dictionary of each line, starting from 1 with each column held exactly as is - With spaces.
        If clean=1, returns a dictionary without whitespace.
        
        Chains is list of chains you would like to extract from the PDB file.
        '''
        
        i = 1
        for line in lines:
            line = line.strip('\n')
            lineSP = line.split()
            if lineSP[0]=="ATOM":

                if chains==0:
                    if pdbDic.has_key(i):
                        pdbDic[i]["id"]="ATOM  "
                    else:
                        pdbDic[i] = dict()
                        pdbDic[i]["id"]="ATOM  "
                    if clean ==0:
                        pdbDic[i]["atomNumber"]=line[6:11]
                        pdbDic[i]["atomName"] = line[11:16]
                        pdbDic[i]["alternateLocation"]=line[16]
                        pdbDic[i]["ResName"] = line[17:21]
                        pdbDic[i]["chain"] = line[21]
                        pdbDic[i]["resNum"]= line[22:26]
                        pdbDic[i]["iCode"] = line[26]
                        pdbDic[i]["x"] = line[27:38]
                        pdbDic[i]["y"]= line[38:46]
                        pdbDic[i]["z"]= line[46:54]
                        pdbDic[i]["occupancy"] = line[54:60]
                        pdbDic[i]["Bfactor"]=line[60:66]
                        pdbDic[i]["element"]=line[66:78]
                        pdbDic[i]["charge"]=line[78:79]
                    elif clean ==1:
                        pdbDic[i]["atomNumber"]=line[6:11].strip()
                        pdbDic[i]["atomName"] = line[11:16].strip()
                        pdbDic[i]["alternateLocation"]=line[16].strip()
                        pdbDic[i]["ResName"] = line[17:21].strip()
                        pdbDic[i]["chain"] = line[21].strip()
                        pdbDic[i]["resNum"]= line[22:26].strip()
                        pdbDic[i]["iCode"] = line[26].strip()
                        pdbDic[i]["x"] = line[27:38].strip()
                        pdbDic[i]["y"]= line[38:46].strip()
                        pdbDic[i]["z"]= line[46:54].strip()
                        pdbDic[i]["occupancy"] = line[54:60].strip()
                        pdbDic[i]["Bfactor"]=line[60:66].strip()
                        pdbDic[i]["element"]=line[66:78].strip()
                        pdbDic[i]["charge"]=line[78:79].strip()
                else:
                    if re.search(line[21].strip(), chains):
                        if pdbDic.has_key(i):
                            pdbDic[i]["id"]="ATOM  "
                        else:
                            pdbDic[i] = dict()
                            pdbDic[i]["id"]="ATOM  "
                        pdbDic[i]["atomNumber"]=line[6:11]
                        pdbDic[i]["atomName"] = line[11:16]
                        pdbDic[i]["alternateLocation"]=line[16]
                        pdbDic[i]["ResName"] = line[17:21]
                        pdbDic[i]["chain"] = line[21]
                        pdbDic[i]["resNum"]= line[22:26]
                        pdbDic[i]["iCode"] = line[26]
                        pdbDic[i]["x"] = line[27:38]
                        pdbDic[i]["y"]= line[38:46]
                        pdbDic[i]["z"]= line[46:54]
                        pdbDic[i]["occupancy"] = line[54:60]
                        pdbDic[i]["Bfactor"]=line[60:66]
                        pdbDic[i]["element"]=line[66:78]
                        pdbDic[i]["charge"]=line[78:79]
                    else:
                        #print "Did not find what we needed...."
                        pass
                i +=1
        return pdbDic
        
        
    def parsePDB(self, filename, pdbDic, clean=0, chains=0):
        '''
        This Reads the PDB using 'loadPDB' and parses it using 'pdbData'
        returns a pdb Dictionary data structure.
        '''
        
        lines = self.loadPDB(filename)
        pdbDic = self.pdbData(lines, pdbDic, clean, chains)
        print "File Loaded into Memory..."
        return pdbDic
    
    def RemEle(self, pdbDic):
        '''
        Removes the extra stuff in the element column, but not the element itself.
        '''
        if not pdbDic:
            print "Please load PDB data into Python to work on...."
            return
        for i in range(1, len(pdbDic)+1):
            ele = pdbDic[i]["element"]
            e = ele[11]
            pdbDic[i]["element"]="           "+e
        print "Extra stuff in Element Columns Removed"
        return pdbDic
    def RemAlt(self, pdbDic):
        '''
        Removes any alternate residue codes and renumbers according to the start of each chain.  If the insertion is at the
        beginning of a chain, starts at residue number one.
        Returns a fixed pdb Dictionary data structure.
        NOT COMPLETELY TESTED.  USE WITH CAUTION.
        '''
        if not pdbDic:
            print "Please load PDB data into Python to work on..."
            return
        #pdbDic = self.pdbData(self.loadPDB(filename))
        inList = dict()
        switchFound = "off"
        chainStart = dict
        switchChain  = pdbDic[1]["chain"]
        code = 0
        addI = 0
        for i in range(1, len(pdbDic)+1):
            if pdbDic[i]["iCode"].rfind(" ")==-1:
                if switchChain == pdbDic[i]["chain"]:
                    switchFound = "on"
                    
                #Insertion = pdbDic[i]["chain"]+":"+repr(i) + pdbDic[i]["resNum"]
                #pdbDic[i]["resNum"] = "I"; pdbDic[i]["iCode"]=" "
                #inList[pdbDic[i]["chain"]] = "on"
                
                if code != pdbDic[i]["iCode"]:
                    if switchChain!= pdbDic[i]["chain"]:
                        print "Found residue insertion in begginning of chain..."
                        print "Starting chain at residue 1..."
                        print "Please Double check results, as this part has not been bug tested!!!"
                        one = 1
                        pdbDicp[i]["resNum"] = str(1).rjust(4)
                        addI = 1
                    else:
                        preNum = pdbDic[i-1]["resNum"]
                        preNum = int((preNum.split())[0])
                        pdbDic[i]["resNum"] = (str(preNum+1).rjust(4))
                        addI+=1
                    
                    print "Found residue insertion at " +pdbDic[i]["resNum"]
                    print ":"+ pdbDic[i]["iCode"]+":"    
                else:
                    pdbDic[i]["resNum"] = pdbDic[i-1]["resNum"]
                code =pdbDic[i]["iCode"]
                pdbDic[i]["iCode"]= " "
            else:
                code = 0
                
                if switchFound == "on" and switchChain ==pdbDic[i]["chain"]:
                    preNumOr = pdbDic[i]["resNum"]
                    preNumOr = int((preNumOr.split())[0])
                    newNum = preNumOr +addI
                    pdbDic[i]["resNum"] = str(newNum).rjust(4) 

            if switchChain != pdbDic[i]["chain"]:
                switchFound = "off"
                addI = 0
                

                

                
            

            switchChain = pdbDic[i]["chain"]
            
            
            
            
            
                
                
                
                
                
        print "Residue Renumbering Complete...."
        return pdbDic
    
    def chaOcc(self, pdbDic):
        '''
        Changes ALL occupancies in a PDB dictionary to 1.00
        Returns PDB Dictionary.
        '''
        if not pdbDic:
            print "Please load PDB file into Python to work on..."
            return
        check = 0
        for key in pdbDic:
            if pdbDic[key]["occupancy"].rfind("0.00")!=-1:
                print "Changing atom occupancy for residue " + pdbDic[key]["resNum"] + "To 1.00"
                check =1
            pdbDic[key]["occupancy"] = "  1.00"
        if check ==1:
            print "Occupancy Column OK for PyRosetta..."
        return pdbDic
    
    def savePDB(self, pdbDic, filename=0):
        '''
        Uses a PDB Dictionary to save the data as a PDB file.
        No END or TER is added, but should be in the future.
        '''
        if not pdbDic:
            print "No PDB file loaded...."
            return
        
        elements = ("id", "atomNumber", "atomName", "alternateLocation", "ResName", "chain", "resNum", "iCode", "x", "y", "z", "occupancy", "Bfactor", "element", "charge")
        dir = os.getcwd()
        if filename==0:
            FILE = tkFileDialog.asksaveasfile(mode = 'w', initialdir=dir,title='Save As...')
        else:
            FILE = open(filename, 'w')
        for key in pdbDic:
            for el in elements:
                FILE.write(pdbDic[key][el])
            FILE.write("\n")
        FILE.close()
        print "PDB File Written..."
        #if retFile !=0:
            #return 
        return
    
    def getSeq(self, file):
        pdbDic = dict()
        pdbDic = self.parsePDB(file, pdbDic)
        res = 0
        seqString = ""
        for line in pdbDic:
            resname = (pdbDic[line]["ResName"]).split()
            resname = resname[0]
            resnum = (pdbDic[line]["resNum"]).split()
            resnum = resnum[0]
            if res != resnum:
                seqString = seqString+ self.OneLetterCode[resname]
            res = resnum
        return seqString
            
            
        #Goes through each line, gets Seq