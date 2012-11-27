from rosetta import *
from Tkinter import *
import tkMessageBox
import tkSimpleDialog

class Docking:
    """
    Basic frontend for interacting with rosetta docking protocols
    Though, only highres is here since lowres needs more decoys then it can handle, and the user should not try to use it in the GUI.
    """
    def __init__(self, pose, score_class, GUIOutput_class):
        self.pose = pose
        self.score_class = score_class
        self.output_class = GUIOutput_class
        
        
    def high_res_dock(self):
        """
        Uses DockMCMProtocol for high-res docking.  Checks scorefunction for docking patch.
        """
        
        to_dock = tkSimpleDialog.askstring(title = "Input", prompt="Please supply two partners (ex - A_B or A_BC )")
        to_dock = to_dock.upper()
        if self.score_class.ScorePatch.get()!="docking":
            result = tkMessageBox.askokcancel(title="Continue?", message="Docking patch not set in scorefunction. Continue?")
            print result
            if not result:return
        jd=PyJobDistributor(self.output_class.outdir.get() + "/" + self.output_class.outname.get(), int(self.output_class.decoys.get()), self.score_class.score);
        dock_mover = DockMCMProtocol()
        dock_mover.set_scorefxn(self.score_class.score)
        dock_mover.set_partners(to_dock)
        
        if self.output_class.auto_write.get():
            for i in range(1, self.output_class.decoys.get()+1):
                print "Decoy number "+repr(i)
                for x in range(1, self.output_class.rounds.get()+1):
                    print "Round number "+repr(i)
                    dock_mover.apply(self.pose)
                jd.output_decoy(self.pose)
        else:
            for x in range(1, self.output_class.rounds.get()+1):
                dock_mover.apply(self.pose)
        print "Complete..."    
        