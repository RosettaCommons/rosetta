from rosettautil.util import fileutil
import sys

class RosettaLoop:
    def __init__(self):
        """initialize an empty loop oject"""
        self.start = 0
        self.end = 0
        self.cutpoint = 0
        self.skip = 0.0
        self.extend = False
        
    def set_loop(self,start,end,cutpoint,skip,extend):
        """fill the loop with data"""
        self.start = start
        self.end = end
        self.cutpoint = cutpoint
        self.skip = skip
        self.extend = extend
        
    def set_loop_from_string(self,loopstring):
        """fill a loop from a line of a Rosetta3 loop file"""
        loop_array = loopstring.split()
        if len(loop_array) != 6:
            sys.exit("loop lines must be in this form: LOOP start end cutpoint skip extend")
        try:
            self.start = int(loop_array[1])
            self.end = int(loop_array[2])
            self.cutpoint = int(loop_array[3])
            self.skip = float(loop_array[4])
        except ValueError:
            sys.exit("start, end, cutpoint and skip all need to be numbers")
        if loop_array[5] =="0":  #this is the way it works in C++, we'll mimic it here for compatability
            self.extend = False
        else:
            self.extend = True
    
    def to_string(self):
        """output a loop as a line of a rosetta3 loop file"""
        if self.extend:
            return "LOOP "+str(self.start)+" "+str(self.end)+" "+str(self.cutpoint)+" "+str(self.skip)+" "+str(1)
        else:
            return "LOOP "+str(self.start)+" "+str(self.end)+" "+str(self.cutpoint)+" "+str(self.skip)+" "+str(0)

class RosettaLoopManager:
    def __init__(self): 
        self.looplist = []
    
    def __iter__(self):
        return iter(self.looplist)
        
    def read(self,filename,append=False):
        """read a rosetta 3 loop file into the loop manager.
        if append=True, add the contents of the loop file to the existing loops"""
        if not append:
            self.looplist = []
        loop_file = fileutil.universal_open(filename,"rU")
        for line in loop_file:
            fields = line.split()
            if line[0] == '#':
                continue #this is a comment
            if len(fields) <1 :
                continue #this is a blank line
            if fields[0] != "LOOP":
                continue #this is something that is not a loop line
            loop = RosettaLoop()
            loop.set_loop_from_string(line)
            self.looplist.append(loop)
        loop_file.close()
    
    def write(self,filename):
        """write the contents of the loop manager to a Rosetta3 loop file"""
        loop_file = fileutil.universal_open(filename,"w")
        for loop in self.looplist:
            loop_file.write(loop.to_string()+"\n")
        loop_file.close()
        
    def is_res_in_loop(self,resnum):
        """Return true if a given residue is inside a loop"""
        for loop in self.looplist:
            if loop.start < resnum and loop.end > resnum:
                return True
        return False
        
    def add_loop(self,start,end,cutpoint,skip,extend):
        """Add a new loop to the loop manager""" 
        loop = RosettaLoop()
        loop.set_loop(start,end,cutpoint,skip,extend)
        self.looplist.append(loop)
    
    
