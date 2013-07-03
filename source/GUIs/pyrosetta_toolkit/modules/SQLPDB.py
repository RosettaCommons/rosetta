#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/modules/PDB.py
## @brief  A PDB class for reading PDBs into an sqlite3 database.  
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#Python Imports
import os
import urllib2

try:
    import sqlite3
except ImportError:
    print "Please reinstall python with sqlite3 support. "
    exit()
    



class SQLPDB:
    def __init__(self, pdbID, modelID, structID, memory=False, path=False):
        """
        modelID is basically a label
        Path specifies the db to load.
        """
        
        if memory:
            self.db = sqlite3.connect(":memory:") #Can always switch this by using a.db = "adsfasdf" to connect to a different database and add information into it.  
        else:
            if path:
                if not os.path.exists(os.path.dirname(path)):
                    os.mkdir(os.path.dirname(path))
            else:
                path = "test_db.db"
            self.db = sqlite3.connect(path)
            
        self.db_util = PDB_database(self.db)
        self.set_basic_options(pdbID, modelID, structID)
        self.pdb_url = "http://www.rcsb.org/pdb/files"
        
        
    def __exit__(self):
        self.db.close()
        
    def set_basic_options(self,pdbID, modelID, structID):
        self.pdbID = pdbID
        self.modelID = modelID
        self.structID = structID
        
    def set_structID(self,structID):
        self.structID = structID
    def set_pdbID(self, pdbID):
        self.pdbID = pdbID
    def set_modelID(self, modelID):
        self.modelID = modelID
        
    def read_pdb_into_database_flat(self, filePath, specific_chain=False, read_header=False, header_only=False, ):
        """
        Reads the flat filepath specified into a database structure.
        This can then be parsed using the PDB_Database class.
        NOTE: Reading of header not implemented.
        if header_only is True, only loads the header.  Useful for just getting specific information.  More useful to D/L it from the pdb if possible.
        If Header only, reads the header into the database.  
        """
        if filePath == "PDB":
            print "Fetching "+self.pdbID+" from the PDB"
            FILE = urllib2.urlopen(self.pdb_url+'/'+self.pdbID.lower()+'.pdb')
        else:
            FILE = open(filePath)
        line_num = 1
        with self.db:
            cur = self.db.cursor()
            
            l = 1
            if read_header:
                cur.execute("CREATE TABLE IF NOT EXISTS header(id integer PRIMARY KEY, pdbID TEXT, modelID TEXT, method TEXT, resolution REAL, species TEXT, engineered TEXT, protein TEXT)")
            if not header_only:
                cur.execute("CREATE TABLE IF NOT EXISTS pdb(l integer PRIMARY KEY, pdbID TEXT, modelID TEXT, strucID INT, type TEXT, atomNum INT, atomName TEXT, altLoc TEXT, residue TEXT, chain TEXT, resNum INT, icode TEXT, x REAL, y REAL, z REAL, occupancy REAL, bfactor REAL)")
            print "Tables created.  Loading "+self.pdbID+" data into table."
            for line in FILE:
                line = line.strip()
                lineSP = line.split()
                if not header_only:
                    if (lineSP[0]=="ATOM" or lineSP[0]=="HETATM"):
                        
                        #Only copy a specific chain into the database.
                        if not specific_chain:
                            pass
                        else:
                            if specific_chain != line[21].strip():
                                continue
                            
                            
                        cur.execute("INSERT INTO pdb VALUES(NULL, ?, ?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", \
                        (self.pdbID, self.modelID, self.structID, lineSP[0], \
                         
                        
                        #atomNum INT             atomName TEXT             altLoc TEXT           residue TEXT
                        line[6:11].strip(),      line[12:16],               line[16],             line[17:21].strip(), \
                        #chain TEXT              resNum INT                icode TEXT            x REAL
                        line[21].strip(),        line[22:26].strip(),      line[26],             line[27:38].strip(), \
                        #y REAL                  z REAL                    occupancy REAL        bfactor REAL
                        line[38:46].strip(),     line[46:54].strip(),      line[54:60].strip(),  line[60:66].strip()))
                        l+=1
                        
                        
                        
                        #self.stripped_pdb[line_num]["element"]=line[66:78].strip();        self.stripped_pdb[line_num]["charge"]=line[78:79].strip())
                      
    
    def fetch_and_read_pdb_into_database(self, pdbID, read_header=False, header_only=False):
        """
        Uses the PDB file specified, grabs it from the PDB, and reads the data in.
        """
        self.pdbID = pdbID
        self.read_pdb_into_database_flat("PDB", read_header, header_only)
        
    




                    
class PDB_database:
    """
    This class is specifically for if we already have a database.
    Note:  This is not a ROSETTA database.  If you need to convert this, use ROSETTA (Which now works in PyRosetta!)
    Functions are to output the database as a PDB, output specific pieces of protein as a pdb and query the database.
    """
    
    def __init__(self, database):
        if type(database) is str:
            self.db = sqlite3.connect(database)
        else:
            self.db = database
        self.db.row_factory = sqlite3.Row #Sets sqlite3 to return values based on a nice list of dictionaries for each row.  Pretty awesome. 
        self._reset_cursor()
        self.occupancy_1 = False
        #MUST QUERY THE TABLE INTERESTED IN FIRST!!!
        
    def _set_db(self, database):
        self.db = database
    
    def _reset_cursor(self):
        self.cur = self.db.cursor()
        
    def _convert_rosetta_db_to_basic_db(self):
        """
        This MAY eventually convert a Rosetta PDB database into a basic database structure that I am using.  
        """
        pass
    def set_output_DIR(self, outDIR):
        self.outDIR = outDIR
        if not os.path.exists(self.outDIR):
            os.mkdir(self.outDIR)
    
    def set_output_occupancy_1(self, bool):
        """
        Output structures with 1.0 as occupancy.  Mainly for Rosetta use.
        """
        self.occupancy_1 = bool
        
#################Query on the Current Cursor of the PDB Database#############################################################
    
    def query_all(self, table="pdb"):
        self.cur.execute("SELECT * FROM "+table)
        
    def query_modelID(self, table, modelID):
        self.cur.execute("SELECT * FROM "+table+" WHERE modelID=?", (modelID,))
    
    def query_pdbID(self, table, pdbID):
        self.cur.execute("SELECT * FROM "+table+" WHERE pdbID=?", (pdbID,))
    
    def query_strucID(self, table, strucID):
        self.cur.execute("SELECT * FROM "+table+" WHERE strucID=?", (strucID,))
        
    def query_pdbID_and_strucID(self, table, pdbID, strucID):
        self.cur.execute("SELECT * FROM "+table+" WHERE pdbID=? AND strucID=?", (pdbID, strucID))
        
    def query_chain(self, table, chain):
        self.cur.execute("SELECT * FROM "+table+" WHERE chain=?", (chain,))
    
    def query_pdbID_and_chain(self, table, pdbID, chain):
        self.cur.execute("SELECT * FROM "+table+" WHERE pdbID=? AND chain=?", (pdbID, chain))
        
    def query_piece(self, table, start, end, chain):
        self.cur.execute("SELECT * FROM "+table+" WHERE chain=? AND resnum BETWEEN ? AND ?", ( chain, start, end))

    def query_piece_pdbID(self, table, pdbID, start, end, chain):
        self.cur.execute("SELECT * FROM "+table+" WHERE pdbID=? AND chain=? AND resnum BETWEEN ? AND ?", ( pdbID, chain, start, end))
    
    def query_piece_pdbID_and_strucID(self, table, pdbID, start, end, chain, strucID):
        self.cur.execute("SELECT * FROM "+table+" WHERE pdbID=? AND chain=? AND strucID=? and resnum BETWEEN ? AND ?", (pdbID, chain, strucID, start, end))
    
    def scrub(self, table_name):
        """
        This should help protect from sql injection.  Not that it's important now, but...
        Author:OrangeOctopus from stack overflow
        """
        
        return ''.join( chr for chr in table_name if chr.isalnum() )


################Update on current cursor of the Database######################################################################

    def update_modelID_CDRS(self, table):
        """
        Updates modelID to specify L1 through H3 and framework for possible future statistical analysis.
        """
        from modules.Structure import Antibody_Structure
        AB = Antibody_Structure()
        self.cur.execute("UPDATE pdb SET modelID='FRAMEWORK'")
        for cdr in AB.CDRS:   
            self.cur.execute("UPDATE pdb SET modelID=? WHERE chain=? AND resnum BETWEEN ? AND ?", (cdr.name, cdr.chain, cdr.Nter, cdr.Cter))
            
        self.db.commit()
        self._reset_cursor(); #Nessessary as after update, resultant cursor is NULL.
        
################Save PDB database as PDB or database!#########################################################################  
    
    def save_cur_as_pdb(self, outpath, supress_modelSep=False):
        """
        Saves the DB at the current cursor to a file. Make sure cursor is on the pdb table.
        """
        
        FILE = open(outpath, 'w')
        self._convert_and_save_to_FILE_pdbformat(self.cur, FILE, supress_modelSep)
        FILE.close()
    
    def save_whole_db_as_db(self, filename, seperate_structures=False):
        """
        Saves the whole database in MEMORY to a file....
        """
        
        if not seperate_structures:
            outPath = self.outDIR+'/'+filename
        new_db = sqlite3.connect(outPath)
        cur = new_db.cursor()
        cur.execute("ATTACH DATABASE ':memory:' AS pdb_db")
        new_db.close()
    
    def _convert_and_save_to_FILE_pdbformat(self, cur, FILE, supress_model_separation=False):
        """
        Saves all information of the current cur to PDB format.
        """
        

        rows = cur.fetchall()
        if not supress_model_separation:
            pdb_previous = rows[0]['pdbID']; struc_previous = rows[0]['strucID']
            FILE.write("MODEL\n")
            for row in rows:
                pdb=row['pdbID']; struc=row['strucID']
                if pdb!=pdb_previous or struc!=struc_previous:
                    FILE.write("ENDMDL\n")
                    FILE.write("MODEL\n")
                line = self._morph_db_row_to_pdb_line(row)
                FILE.write(line+"\n")
                pdb_previous = row['pdbID']; struc_previous = row['strucID']
            FILE.write('ENDMDL')
        else:
            FILE.write("MODEL\n")
            for row in rows:
                line = self._morph_db_row_to_pdb_line(row)
                FILE.write(line+"\n")
            FILE.write('ENDMDL')
        FILE.close()
        
    def _morph_db_row_to_pdb_line(self, row):
        """
        Oh What fun. ;)
        (6,5,4,3,1,4,8,8,8,4,5);
        atomNum INT, atomName TEXT, altLoc TEXT, residue TEXT, chain TEXT, resNum INT, icode TEXT, x REAL, y REAL, z REAL, occupancy REAL, bfactor REAL")
        """
        occupancy = 0
        if self.occupancy_1:
            occupancy = 1.0
        else:
            occupancy = row['occupancy']
        #Create the PDB line.
        line = str(row['type']).ljust(6)+     str(row['atomNum']).rjust(5)+" "+str(row['atomName'])+ \
               str(row['altLoc'])+            (str(row['residue']).rjust(3)).ljust(4)+ str(row['chain'])+             \
               str(row['resNum']).rjust(4)+   str(row['icode']) +                                              \
               ("%.3f"%row['x']).rjust(11)+   ("%.3f"%row['y']).rjust(8)+       ("%.3f"%row['z']).rjust(8) +   \
               str(occupancy).rjust(6)+str(row['bfactor']).rjust(6)

        
        return line
    

    
    def _add_new_struct_to_existing_database(self, db, filename):
        """
        Concatonates PDBs in databases. The database should have a unique pdbID+modelID, and not be in Rosetta format
        """
        
        pass
    
if __name__ == '__main__':
    """
    Specifically for testing ATM.
    """
    pdbID = '1hil'
    modelID = 'FULL'
    outDIR = os.path.split(os.path.abspath(__file__))[0]
    s = PDB(pdbID, modelID, "x", False, outDIR+'/testing.db')
    s.fetch_and_read_pdb_into_database()
    
    table =s.db_util.scrub('pdb')
    s.db_util.query_all(table)
    s.db_util.set_output_DIR(outDIR)
    #s.db_util.save_whole_db_as_db("Whole_db.db")
    
    #s.db_util.save_cur_as_pdb("Test_Full_PDB.pdb")
    #s.db_util.save_cur_as_db("Test_Full_DB.db")
    
    #Test Queries.
    
    #s.db_util.query_chain(table, 'H')
    #s.db_util.save_cur_as_pdb("Test_H.pdb")
    #s.db_util.save_cur_as_db("Test_H.db")
    
    #s.db_util._reset_cursor()
    #s.db_util.query_pdbID(table, "2j88")
    #s.db_util.query_modelID(table, "FULL")
    #s.db_util.save_cur_as_pdb("model_query.pdb")
    
    s.db_util._reset_cursor()
    s.db_util.query_chain(table, 'A')
    s.db_util.query_piece(table, 24, 48, 'A')
    s.db_util.save_cur_as_pdb("H1A_CDR.pdb")
    
    #Loading Local Files.
    
    