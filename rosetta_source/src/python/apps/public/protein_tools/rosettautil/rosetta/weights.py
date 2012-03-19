from rosettautil.util import fileutil
aa_codes_in_order = ["ALA","CYS","ASP","GLU","PHE","GLY","HIS","ILE","LYS","LEU","MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR"]

class WeightFile:
    def __init__(self):
        self.ref_energies = {}
        self.weights = {}
    def read_file(self,filename):
        self.ref_energies = {}
        self.weights = {}
        
        in_file = fileutil.universal_open(filename,'rU')
        for line in in_file:
            line = line.split()
            if line[0] == "METHOD_WEIGHTS":
                #Reference energies are ordered by 1 letter name
                for aa, value in zip(aa_codes_in_order,line[2:len(line)]):
                    self.ref_energies[aa] = float(value)
            else:
                self.weights[line[0]] = float(line[1])
        in_file.close()
    
    def write_file(self,filename):
        out_file = fileutil.universal_open(filename,'w')
        
        #write reference energies
        out_file.write("METHOD_WEIGHTS\tref")
        for key in aa_codes_in_order:
            out_file.write("\t"+str(self.ref_energies[key]))
        out_file.write("\n")
        
        #write the other weights
        for key in self.weights:
            out_file.write(key+"\t"+str(self.weights[key])+"\n")
        out_file.close()
    
    def get_ref(self,aa):
        return self.ref_energies[aa]
    
    def get_weight(self,term):
        return self.weights[term]
    
    def set_ref(self,aa,value):
        self.ref_energies[aa] = value
    
    def set_weight(self,term,value):
        self.weights[term] = value
