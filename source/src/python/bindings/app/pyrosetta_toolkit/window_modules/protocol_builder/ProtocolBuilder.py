#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/window_modules/protocol_builder/ProtocolBuilder.py
## @brief  Original protocol builder.  DEPRECATED.  being replaced by rosetta script creator.
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

from rosetta import *
from Tkinter import *
import glob
import tkFileDialog
import tkMessageBox
import tkSimpleDialog


#Protocol Variables (Main)
#Difinitions for PROTOCOLS - SHOULD CHANGE!!!
#THESE SHOULD BE IN SEPERATE FILE!!

ProtTypesDic = dict()
ProtTypes = ("Low-Res(Frag)", "High-Res", "Loop Minimization", "Full Minimization", "Vaccinity", "Interface", "ShearMoves", "Design(ResFile)", "CDR-Specific", "Tools")
ProtTypesDic["ShearMoves"]=("Shear Custom", "Shear CustomFile", "Shear 1", "Shear 2", "Shear 3", "Shear 4", "Shear 5", "Shear 10", "Shear 15", "Shear 20", "Shear 30", "Shear 50", "Shear 80", "Shear 100", "Shear 120", "Shear 150", "Shear 180")
ProtTypesDic["Low-Res(Frag)"] = ("Low: Pure Frag Anneal", "Low-Perturb CCD: Default", "Low-Perturb KIC")
ProtTypesDic["High-Res"] = ("Hi-Refine CCD: Default", "Hi-Refine KIC")
ProtTypesDic["Vaccinity"] = ("Pack Vaccinity", "Pack Vaccinity (SCWRL)", "Relax Vaccinity", "LowPureFragAnneal Vaccinity", "Minimize Vaccinity")
ProtTypesDic["Interface"] = ("Pack Interface", "Pack Interface (SCWRL)", "Relax Interface", "Minimize Interface")
ProtTypesDic["Tools"] = ("CCD Loop Closure", "Linearize Loops", "Randomize Loops")
ProtTypesDic["Full Minimization"]=("Fast Relax", "Fast Relax BB", "Classic Relax", "Backrub", "Minimize", "Optimize Rot", "SCWRL")
ProtTypesDic["Loop Minimization"]=("Classic Loop Relax", "Fast Loop Relax", "Fast Loop Relax BB", "Loop Backrub", "Loop Minimize", "Loop Optimize Rot", "SCWRL(Loops)", "SCWRL(Seq file)")
ProtTypesDic["CDR-Specific"]=("EMPTY", " ")
ProtTypesDic["Design(ResFile)"] = ("PackDesign", "Pack Vaccinity Design", "Pack Interface Design")

#Need to Set Prot and Set Frag upon startup....
#Have 

#Adds Data into protocol list.
def getProtocolLists(LisBox):
    VarItem = LisBox.get(LisBox.curselection())
    x = len(DicProtocol)
    x+=1
    #Creates Dictionaries to tie them all together
    p = "p" + repr(x)
    print p
    #OLD WAY: DicProtocol[p] = [VarItem, int(VarRounds.get()), int(Varcheck_button_Cen.get()), int(VarPym.get()), int(VarFragLen.get())]
    DicProtocol[p]=dict()
    #NEW WAY: LOAD INTO DICTofDICT
    #NEWNEW way is to have a class hold all of these in an object!!!
    #[name] = VarItem(name)
    #[rounds] = Rounds (integer)
    #[centroid]= Centroid (integer) 0 or 1
    #[obs] = PyMOL Observer (integer) 0 or 1

    #[shear] = Shear Settings
    #[lisloops] = Loop Settings
    #[score] = Score Function!!! (Literally)
        #Accessory Files
    #[frag] = Fragment Length (int) 9 or 3.  Should be context dependendant.
    #[fragfile] = Accessory File - Will be actual frag file soon.   
        #Vaccinity:
    #[neighborcuttoff] = Neighbor Distance for vaccinity measure
    #[onlyaround] = 0 or 1.  1 means only around vaccinity, not including residues
    #[fixloopbb] = 0 or 1.  1 indicates to keep bb of loop fixed during Relax Vaccinity
    #[fixvacbb] = 0 or 1.  1 indicates to keep bb of vaccinity fixed during Relax Vaccinity.    
        #Interface
    #[interfacecutoff] = Interface Cuttoff distance

    #[interfacefixed] = Chains to keep fixed in Interface Protocols: (A-B-C)
    #[interfacepack] = Chains to only allow repacking for Relax Protocols

    
    #Only Around Combinations: (onlyaround, fixloopbb, fixvacbb)
        #0, 0, 0: Relax All
        #0, 1, 0: Relax Vaccinity, Keep Loop chi open
        #0, 0, 1: Relax Loop, Keep vaccinity chi open
        #1, 0, 0: Relax Vaccinity, Keep Loop chi closed
        #Other : All default to Pack Vaccinity - So we change the name to Pack Vaccinity.
        
    DicProtocol[p] = dict()

    #Asks User to input the number of Rounds Desired.
    rounds = tkSimpleDialog.askinteger(title = 'rounds', prompt = 'Please enter the number of rounds desired: ', initialvalue = 1)
    VarRounds.set(repr(rounds))
    DicProtocol[p]['name']=VarItem
    #OLD: DicProtocol[p]['rounds'] = int(VarRounds.get())
    DicProtocol[p]['rounds'] = rounds
    DicProtocol[p]['centroid'] = int(Varcheck_button_Cen.get())
    DicProtocol[p]['obs']= int(VarPym.get())
    #DicProtocol[p]['frag'] = int(VarFraglen.get())
    #DicProtocol[p]['score] = Get Currently Set Score
    DicProtocol[p]['score'] = score
    
    #Adds the LisLoop - Python has references too, so be careful!
    DicProtocol[p]['lisloops'] = LisLoop[:]
    #This handles protocol specific data
        #Shear Data
    VarItemSplit = VarItem.split()
    
    #Peudofix for the rest of the code where it tests itemsplit[1] etc.
    if len(VarItemSplit) ==1:
        VarItemSplit[1:4] = "0"
        
    #if "Frag":
        #Prompt for FragLength and Choose FragFile
    #This Needs to be protocol specific prompting.
    #Sets ShearSettings from WINDOW Prompt
    if (VarItemSplit[0] == "Shear" and VarItemSplit[1] =="Custom"):
        print "Appending Custom Shear Protocol"
        ShearSettings=[int(Shear1.get()), int(Shear2.get()), int(Shear3.get()), int(Shear4.get()), ShearKT.get(), shearAnneal.get()]
        
        DicProtocol[p]['shearsettings']=ShearSettings[:]
    elif (VarItemSplit[0] == "Shear" and VarItemSplit[1]=="CustomFile"):
        print "Appending Accessory Shear File"
    elif (VarItemSplit[0] == "Shear"):
        print "Appending default shear settings..."
        ShearSettings = []
        for i in range(0, 4):
            ShearSettings.append(VarItemSplit[1])
        ShearSettings.append(ShearKT.get()); ShearSettings.append(shearAnneal.get())
        DicProtocol[p]['shearsettings'] = ShearSettings[:]
        #OLD: DicProtocol[p].append(ShearSettings)
    
    
    #Asks user for Cutoff Length and Vaccinity Settings for each protocol.  
    if VarItemSplit[1] == "Vaccinity":
        
        length = tkSimpleDialog.askfloat(title = "vaccinity", prompt = "entry_er Vaccinity Cuttoff (A)", initialvalue = 5.0)
        DicProtocol[p]['neighborcutoff'] = length
        if VarItemSplit[0]=="Pack" and VarItemSplit[1] == "Vaccinity":
            onlyaround = tkMessageBox.askquestion(title = "onlyaround", message = "Pack/Relax Given residues in addition to vaccinity?", default = tkMessageBox.YES)
            if onlyaround == "yes":
                print "Packing Vaccinity and residues chosen..."
                DicProtocol[p]['onlyaround']=0
            else:
                print "Only packing Vaccinity...."
                DicProtocol[p]['onlyaround']=1
            
        #Relax and Minimize Vaccinity Prompts
    if VarItem == "Relax Vaccinity" or VarItem == "Minimize Vaccinity":
            #OLD: fixBB = tkMessageBox.askquestion(title = "fix surrounding", message = "Fix Backbone of Vaccinity?", default = tkMessageBox.YES)
            #OLD: if fixBB == 'yes':
                #OLD: DicProtocol[p]['fixvacbb'] = 1
            #OLD: else:
                #OLD: DicProtocol[p]['fixvacbb'] = 0
            #OLD: if onlyaround == "yes":
                #OLD: fixBBLoop = tkMessageBox.askquestion(title = 'Loop', message = "Fix Backbone of Loop Residues?", default = tkMessageBox.NO)
                #OLD: if fixBBLoop == 'yes':
                    #OLD: DicProtocol[p]['fixloopbb'] = 1
                #OLD: else:
                    #OLD: DicProtocol[p]['fixloopbb'] = 0
            #OLD: if (onlyaround == "no" and fixBB == 'yes'):
                #OLD: tkMessageBox.showinfo(message = "BB vaccinity fixed, no packing of loop residues - Runnng Rosetta Pack Residues instead of relax.")
                #OLD: DicProtocol[p]['name']="Pack Vaccinity"
            #OLD: if onlyaround == 'yes' and fixBB =='yes' and fixBBLoop == 'yes':
                #OLD: tkMessageBox.showinfo(message = "BB vaccinity fixed for loop and surrounding residues - Running Rosetta Pack Residues instead of relax.")
                #OLD: DicProtocol[p]['name']="Pack Vaccinity"
        
        vacwindow= vaccinityrelaxwindow(c)
        #c.wait_window(vacwindow)
        
        if vacwindow.result == None:
            print "Canceled.  Using Default settings for protocol."
            FixTarget = "Open"; FixVaccinity = "Open"; FixBoth = "UnSet"
        else:
            (FixTarget, FixVaccinity, FixBoth) = vacwindow.result
        #Here, we translate what each of these mean (Hopefully):
        
        target = [0, 0]; #(BB, Chi) 0 means to keep; 1 means to fix
        vaccinity = [0, 0]
        
        #BOTH
        if FixBoth =="Vaccinity":
            vaccinity =[1, 1]
            print "Only Relaxing Target.  You Should run Loop Relax instead."
        elif FixBoth == "Target":
            target = [1, 1]
            print "Only Relaxing Vaccinity"
        
        #TARGET
        
        if FixTarget == "Open":
            print "Opening Target"
            target=[0,0]
        elif FixTarget == "Fix":
            print "Fixing all atoms of Target"
            target = [1, 1]
        elif FixTarget =="Fix BB":
            print "Fixing backbone of Target"
            target[0]=1
        elif FixTarget =="Fix Chi":
            print "Fixing rotamers of Target"
            target[1]=1
        
        
        #VACCINITY
        if FixVaccinity =="Open":
            print "Opening Vaccinity"
            vaccinity = [0, 0]
        elif FixVaccinity =="Fix":
            print "Fixing all atoms of Vaccinity"
            vaccinity = [1, 1]
        elif FixVaccinity =="Fix BB":
            print "Fixing backbone of Target"
            vaccinity[0]=1
        elif FixVaccinity =="Fix Chi":
            vaccinity[1]=1
            
            
        print vaccinity; print target
        print "Final Settings for Relax Protocol:\n"
        #TARGET
        if target[0] ==1 and target[1]==1:
            print "Fixing Target Completely"
        elif target[0]==0 and target[1]==0:
            print "Relaxing Target Completely"
        elif target[0]==0:
            print "Relaxing BackBone of Target"
        elif target[1]==0:
            print "Repacking Rotamers of Target"
        
        #VACCINITY
        if vaccinity[0] ==1 and vaccinity[1]==1:
            print "Fixing Vaccinity Completely"
        elif vaccinity[0]==0 and target[1]==0:
            print "Relaxing Vaccinity Completely"
        elif vaccinity[0]==0:
            print "Relaxing BackBone of Vaccinity"
        elif vaccinity[1]==0:
            print "Repacking Rotamers of Vaccinity"
            
        print "\n\n"
        print "Switching Protocols as Needed:"
        if vaccinity == [1, 1] and target == [1, 1]:
            tkMessageBox.showerror(message = "Everything is fixed!  Remove Protocol and Try Again!")
            return
        
        if vaccinity[0]==1 and target[0]==1:
            print "Only running Side Chain Packing Protocol..."
        
        DicProtocol[p]['vaccinity']=vaccinity
        DicProtocol[p]['target'] = target
        #print FixBB
        #print FixChi
        #print FixBoth
        #return
    #Interface Prompts
    print VarItem; print VarItemSplit[1]
    if VarItemSplit[1] == "Interface":
        #Chain Cutoff Distance
        length = tkSimpleDialog.askfloat(title = "interface", prompt = "entry_er Interface Cuttoff (A)", initialvalue = 6.0)
        interface = tkSimpleDialog.askstring(title = "interface", prompt = "entry_er Interface (A:CD / AB:CD:F)")
        #Options if you are specifically looking for interfaces between multiple chains.
            #Handles Fully Restricking some chains
        intopt = tkMessageBox.askquestion(title = "interface", message = "Fix any Interface Chains?", default = tkMessageBox.NO)
        if intopt =="yes":
            interfaceKeep = tkSimpleDialog.askstring(title = "interface", prompt = "entry_er chains to keep fixed (A-B-C) :", initialvalue = "0")
            DicProtocol[p]['interfacefixed'] = interfaceKeep
        else:
            DicProtocol[p]['interfacefixed']=0
        #Handles Restricting some chains to repacking
        if VarItem == "Relax Interface":    
            intopt2 = tkMessageBox.askquestion(title = 'interface', message = "Restrict any Interface Chains to rePacking?", default = tkMessageBox.NO)
            if intopt2 =="yes":
                interfacePack = tkSimpleDialog.askstring(title="interface", prompt = "Enter chains to only allow repacking (A:B:C)")
                DicProtocol[p]['interfacepack']=interfacePack
            else:
                DicProtocol[p]['interfacepack']=0
        DicProtocol[p]['interfacecutoff'] = length
        DicProtocol[p]['interface'] = interface
    #Asks for Tolerance and Type for minmimization
    if VarItemSplit[1] =="Minimize" or VarItemSplit[0] == "Minimize":
        tol = tkSimpleDialog.askfloat(title = "tolerance", prompt = "Please Enter Tolerance", initialvalue = .1)
        DicProtocol[p]['tolerance'] = tol
        
        
    #OLD: print DicProtocol
def getFragLists(LisBox):
    VarItem = LisBox.get(LisBox.curselection())
    VarItem = pwd + "/FragSets_Designs/" +VarItem
    x = len(DicProtocol)

    #Creates Dictionaries to tie them all together
    p = "p" + repr(x)
    print p
    #OLD: DicProtocol[p]['frag']=int(VarFragLen.get())
    #Will check_button_ck if this is a frag or design file.
    length = tkSimpleDialog.askinteger(title="frag", prompt="entry_er Fragment Length", initialvalue=3)
    DicProtocol[p]['frag'] = length
    DicProtocol[p]['fragfile'] = VarItem
    print DicProtocol
def clearLists(List):
    DicProtocol.clear()
    List.delete(0, END)
    print "Protocol Directions Reset"
    
    


#REMOVE PROTOCOL from DicProtocol and ListBox
def RemoveListItem(LisIn):
    #print DicProtocol

    tup = LisIn.curselection()
    num = tup[0]
    #print "Num" + num
    LisIn.delete(LisIn.curselection())
    num = int(num)+1
    key = "p"+repr(num)
    #print "Key" +key
    tot = len(DicProtocol)
    #print "total Protocols"
    end = tot
    start = num-1
    
    #Start ReArranging at Num.  Num represents deleted number +1.
    if num == tot:#Takes Care of End
        DicProtocol.pop(("p"+repr(tot)));
        print "Protocol Removed from Queue"
        return
    else:
        DicProtocol.pop(key)
        
        
    if num >= 1:
        for i in range(num, end):
            #print i
            #print num
            #print "\n"
            DicProtocol[("p"+repr(num))]=DicProtocol[("p"+repr(num+1))]
            num+=1
    DicProtocol.pop(("p"+repr(tot))); #Removes the last item, as this will be copied into the -1.
        
        
    #print DicProtocol
    print "Protocol Removed from Queue"
    


class MainProtocols():
    def __init__(self, main):
        self.main = main
        self.repeats = StringVar()
        self.repeats.set(0)
        self.DicProtocol = dict(); #Main Protocol Dictionary
        self.Centroid = StringVar(); #Centroid check_button_ckbox
        self.Centroid.set(0)
        self.FragLen = StringVar(); #Fragment Size
        self.RoundVar = StringVar(); #Number of Protocol Rounds
        

        
    #ADD Protocol: Places data into Protocol List using getProtocolLists
    def SetandUpdateListsP(self, LisIn, LisOut):
        
        LisItem = LisIn.get(LisIn.curselection())
        getProtocolLists(self.LisProt1)
        LisItem = self.RoundVar.get()+ "-"+LisItem 
        LisOut.insert(END, LisItem)
        
        
        
    #ADD Accessory: Places data into Protocol by concatinating pwd, folder, and item.
    def SetandUpdateListsF(self, LisIn, LisOut):
        LisItem = LisIn.get(LisIn.curselection())
        LisItem = pwd + "/FragSets_Designs/" +LisItem
        LisItem = self.FragLen.get() + "-"+ LisItem
        getFragLists(self.LisProt1Frag)
        
        
    #UPDATES Listboxes of protocols after choosing type:    
    def updateLisProt(self, ListTypes, ListTypesFull):
        ListTypesFull.delete(0, END)
        type = ListTypes.get(ListTypes.curselection())
        #print type
        for res in ProtTypesDic[type]:
            ListTypesFull.insert(END, res)
            
            
    #Shows ShearControl Window.  Should be chosen if/as Shear Control is being added. 
    def shoShearControl(self):
        WinShear = Toplevel(self.main)
        shearCon = ShearControl(WinShear)
        shearCon.setTk()
        shearCon.shoTk()


    def kickProtocol(self):
        """
        Handles kicking off the protocol with number of repeats.
        """
        repeat = int(self.repeats.get())
        count=1
        
        
        if p.total_residue() == 0:
            tkMessageBox.showerror(message = "No Pose Loaded!  ")
            return
        decoys.set(tkSimpleDialog.askstring(title = "decoys", prompt = "Please entry_er the Number of Decoys Desired:", initialvalue = decoys.get()))
        if int(decoys.get()) > 1:
            print dirnameout.get()
            if dirnameout.get() == "0":
                dirnameout.set(tkFileDialog.askdirectory(initialdir = pwd, title = "Pick directory for decoys: "))
                Outname.set(tkSimpleDialog.askstring(title = "decoy name", prompt = "Please enter a base name for the decoys: ", initialvalue = output.entry_Outname.get()))
                if (dirnameout.get()==None)| (Outname.get()==None):
                    return
            else:
                print "Saving decoys to: " + dirnameout.get()
                print "Base name for decoys: " + output.entry_Outname.get()
        if repeat >=1:
            print "Repeating protocol..."
            print repeat
            start = "1"+"-"+repr(len(self.DicProtocol))
            numbers = tkSimpleDialog.askstring(title = "Simple Repeats", prompt = "Please Specify which protocols to repeat (1-4)", initialvalue = start)
            numbers = numbers.split("-")
            start = int(numbers[0]); end = int(numbers[1])
            length = len(self.DicProtocol)
            for i in range(1, repeat):
                print i
                for x in range(start, end+1):
                    print x
                    newkey = count+length
                    print "newkey:"+repr(newkey)
                    newkey = "p"+repr(newkey)
                    oldkey = "p"+repr(x)
                    print "oldkey:"+repr(oldkey)
                    self.DicProtocol[newkey] = self.DicProtocol[oldkey]
                    count+=1
            print self.DicProtocol
            p.assign(general_tools.protocols().initLoopProtocols(p, self.DicProtocol, decoys.get(), dirnameout.get() + "/" + output.entry_Outname.get(), LisLoop))
        else:
            p.assign(general_tools.protocols().initLoopProtocols(p, self.DicProtocol, decoys.get(), dirnameout.get() + "/" + output.entry_Outname.get(), LisLoop))
            #This is the default kick protocol.
            
    #DELETES Accessory File
    def delAccessory(self):
        acc = self.LisProt1Frag.get(self.LisProt1Frag.curselection())
        file = pwd +"/FragSets_Designs/"+acc
        os.remove(file)
        self.LisProt1Frag.delete(self.LisProt1Frag.curselection())
        
    #Prints Protocol Information:
    def shoProt(self):
        tup = self.LisSeeProt.curselection()
        num = tup[0]
        num = int(num)+1
        key = "p"+repr(num)
        print key
        for settings in self.DicProtocol[key]:
            print settings + "  :  "+repr(self.DicProtocol[key][settings])
        Lmode.LisLoops.delete(0, END)
        #print "Start: "+repr(LisLoop)
        for i in range(0, len(LisLoop)):
            LisLoop.pop(i)
        #print "Deleted "+ repr(LisLoop)
        for loops in self.DicProtocol[key]['lisloops']:

            Lmode.LisLoops.insert(END, loops)
            LisLoop.append(loops)
        #print "Added " + repr(LisLoop)
        return LisLoop
    def editProt(self):
        """
        Does Not WORK!!!?
        """
        
        
        index = self.LisSeeProt.curselection()
        LisItem = self.LisSeeProt.get(index)
        tup = self.LisSeeProt.curselection()
        num = tup[0]
        num = int(num)+1
        key = "p"+repr(num)
        print key
        getProtocolLists(self.LisSeeProt)
        LisItem = self.RoundVar.get()+ "-"+LisItem
        self.LisSeeProt.insert(index, LisItem)
    def makeWindow(self, main):
        self.check_button_Cen = Checkbutton(self.main, text="Centroid Mode", variable=self.Centroid)
        self.label_Ed1 = Label(self.main, text="Protocol Builder", font="Arial")
        self.label_Dire= Label(self.main, text="--Add protocol then add Accessory File--")
        self.button_ProtAdd = Label(self.main, text="Rounds")
        
        #OLD: self.entry_Rounds= Entry(self.main, textvariable=self.RoundVar, justify=CENTER)
        self.RoundVar.set("0")
        self.entry_FragLen= Entry(self.main, textvariable=self.FragLen, justify=CENTER)
        self.FragLen.set("3")
        self.button_Sho = Button(self.main, text = "Show Protocol Settings", command = lambda: self.shoProt())
        self.button_Edi = Button(self.main, text = "Edit Protocol Settings", command = lambda: self.editProt())
        #OLD: self.button_FragAdd = Label(self.main, text="Frag Length")
        self.button_ProtRes = Button(self.main, text="Reset Protocol", command=lambda: clearLists(self.LisSeeProt))
        self.LisProt1 = Listbox(self.main); self.LisProt1Frag = Listbox(self.main)
        self.LisSeeProt = Listbox(self.main); self.LisProtTypes = Listbox(self.main)
        
        #Repetition of Protocol
        self.entry_Repeat = Entry(self.main, justify=CENTER, textvariable = self.repeats)
        self.label_Repeat = Label(self.main, text = "Repeats")
        
        #Sets Bindings
        #Choose Type:
        self.LisProtTypes.bind("<ButtonRelease-1>", lambda event: self.updateLisProt(self.LisProtTypes, self.LisProt1))
        #Choose Protocol + Frag
        self.LisProt1.bind("<Double-Button-1>", lambda event:self.SetandUpdateListsP(self.LisProt1, self.LisSeeProt))
        self.LisProt1Frag.bind("<Double-Button-1>", lambda event:self.SetandUpdateListsF(self.LisProt1Frag, self.LisSeeProt))
        
        
        #Final Listbox
        self.LisSeeProt.bind("<Double-Button-1>", lambda event:RemoveListItem(self.LisSeeProt))
        #OLD: self.LisSeeProt.bind('<ButtonRelease-1>', lambda event:self.shoProt())
        
        #Shear Moves:
        self.label_Shear = Label(self.main, text="Shear Moves")
        self.check_button_Shear = Checkbutton(self.main, text = "Anneal?", variable = shearAnneal)
        ShearKT.set(1.0)
        self.button_Shear = Button(self.main, text = "Custom", command=lambda: self.shoShearControl())
        #(self.LisProt1, self.LisSeeProt)
        self.button_DelAcc = Button(self.main, text = "Delete Accessory File", command = lambda: self.delAccessory())
        
        #Repeat of protocols
        self.label_Rep = Label(self.main, text = "Repeats")
        self.entry_Rep = Entry(self.main, justify = CENTER, textvariable = self.repeats)
        #Kicks Protocol
        self.kickprotocol=Button(self.main, text="Start Protocol", command=lambda: self.kickProtocol())
        
        self.setProt()
        self.setFrag()
        self.shoProt()
        
    def setProt(self):
        """
        Sets up protocol List
        """ 
        
        for type in ProtTypes:
            self.LisProtTypes.insert(END, type)
    def setFrag(self):
        """
        Sets Fragment Listbox up.  Should be able to customize this in the future.
        Should be a customize file that tk main looks for.  If found, it loads it.
        Can save configuration on the fly.  
        """
        files = os.listdir(pwd+"/FragSets_Designs")
        for file in files:
            file = file.split("/")
            self.LisProt1Frag.insert(END, file[len(file)-1])
        #self.LisProt1Frag.insert(END, "CDR_Specific")
        #self.LisProt1Frag.insert(END, "CDR_All")
        #self.LisProt1Frag.insert(END, "Loops_All")
        #self.LisProt1Frag.insert(END, "Sequence_Specific")
        #self.LisProt1Frag.insert(END, "Neighbor_Dependant")
    def shoTk(self):
        self.label_Ed1.grid(row=11, column=3, columnspan=2, pady=15)   
        self.LisProtTypes.grid(row=13, column=3, rowspan=6); self.LisProt1.grid(row=13, column=4, rowspan=6); self.LisProt1Frag.grid(row=15, column=5, rowspan=6)
        #self.button_ProtAdd.grid(row=14, column=2);
        #OLD: self.entry_Rounds.grid(row=13, column=2)
        #self.entry_FragLen.grid(row=15, column=2)
        #OLD: self.button_FragAdd.grid(row=16, column=2)
        self.button_ProtRes.grid(row=17, column=2)
        self.check_button_Cen.grid(row=18, column=2)
        self.label_Dire.grid(row=19, column=3, columnspan=2, sticky=W+E)
        self.LisSeeProt.grid(row=20, column=3, rowspan=6, columnspan=2);
        self.button_DelAcc.grid(row=21, column=5)
        self.button_Sho.grid(row=22, column=5)
        self.button_Edi.grid(row=23, column=5)
        self.label_Shear.grid(row=12, column=5)
        self.check_button_Shear.grid(row=13, column=5)
        self.button_Shear.grid(row=14, column=5)
        self.label_Rep.grid(row=27, column = 3, columnspan = 2)
        self.entry_Rep.grid(row=26, column = 3, columnspan = 2)
        self.kickprotocol.grid(row=24, column=2, sticky=W+E)
