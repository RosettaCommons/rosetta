from rosettautil.util import fileutil
import sys
class list_of_2D_vectors:
    def __init__(self):
        self.records = []

    def add_record(self, first_col,second_col):
        self.records.append((first_col,second_col))
    
    def write_bcl_file(self,path):
        out_file = fileutil.universal_open(path,'w')
        list_header ="bcl::storage::List<bcl::storage::VectorND2<bcl::math::Vector<double>>>"
        vector_header = "bcl::storage::VectorND2<bcl::math::Vector<double>>"
        double_header = "bcl::math::Vector<double>"
        
        out_file.write(list_header+"\n")
        out_file.write(str(len(self.records))+"\n")
        
        for first_col, second_col in self.records:
            out_file.write(vector_header+"\n")
            
            out_file.write(double_header+"\n")
            out_file.write(str(1)+"\n")
            out_file.write(str(first_col)+"\n")
            
            out_file.write(double_header+"\n")
            out_file.write(str(1)+"\n")
            out_file.write(str(second_col)+"\n")
        out_file.close()
    
    def read_bcl_file(self,path):
        print "This function doesn't work yet"
        sys.exit()
        out_file = fileutil.universal_open(path,'r')
        list_header ="bcl::storage::List<bcl::storage::VectorND2<bcl::math::Vector<double>>>"
        vector_header = "bcl::storage::VectorND2<bcl::math::Vector<double>>"
        double_header = "bcl::math::Vector<double>"
        list_scope = False
        vector_scope = False
        double_scope = False