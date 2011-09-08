import subprocess, os, glob, sys

class AddAPSA:

    def __init__( self ):
        self.apsa_path = "/Users/momeara/workspace_momeara/apsa/apsa"
        self.output_fname = "apsa_output.csv"
        self.tmp_fname = "tmp_output.apsa"

    def run_apsa( self, struct_fname, apsa_output_fname ):

        command = "%s %s %s" % \
           ( self.apsa_path,
             struct_fname,
             apsa_output_fname )

        print command
        p = subprocess.Popen(
            command,
            close_fds=True,
            shell=True,
            stdout=subprocess.PIPE)
        output = p.stdout.read()

    def parse_apsa_output( self, struct_fname, apsa_output_fname ):
        data = []
        struct_tag = ".".join(struct_fname.split(os.sep)[-1].split(".")[0:-1])+"_0001.pdb"

        f = open(apsa_output_fname)
        l = f.readline()
        while not l.startswith("# guide points"):
            l = f.readline()

        #move to start of data
        l = f.readline(); l = f.readline()
        
        while l != "\n":
            pos = l[61:64].strip()
            peak = l[67:70].strip()
            sec = l[73:76].strip()
            comment = l[78:].strip()

            data.append( (struct_tag, pos, peak, sec, comment) )
            l = f.readline()

        f.close()
        return data

    def add_to_database( self, output_fname, database_fname ):
        sql = """
DROP TABLE IF EXISTS apsa_raw;
CREATE TABLE apsa_raw(
  struct_fileName TEXT,
  resNum INTEGER,
  primary_struct TEXT,
  secondary_struct TEXT,
  comment TEXT);

CREATE INDEX IF NOT EXISTS apsa_raw_index ON
  apsa_raw( struct_fileName, resNum );

DROP TABLE IF EXISTS apsa;
CREATE TABLE apsa(
  struct_id TEXT,
  resNum INTEGER,
  primary_struct TEXT,
  secondary_struct TEXT,
  comment TEXT,

  FOREIGN KEY         ( struct_id, resNum )
  REFERENCES residues ( struct_id, resNum )
  DEFERRABLE INITIALLY DEFERRED,
  PRIMARY KEY ( struct_id, resNum ) );

.separator ','
.import %s apsa_raw

INSERT OR IGNORE INTO apsa
SELECT
  struct.struct_id,
  apsa_raw.resNum,
  apsa_raw.primary_struct,
  apsa_raw.secondary_struct,
  apsa_raw.comment
FROM
  structures as struct,
  apsa_raw
WHERE
  struct.fileName = apsa_raw.struct_fileName;

DROP TABLE apsa_raw;

""" % output_fname

        command = "sqlite3 %s <<HERE_BLOCK\n%s\nHERE_BLOCK" % ( database_fname, sql )
        p = subprocess.Popen(
            command,
            close_fds=True,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        print p.stderr.read()
        print p.stdout.read()

    def run( self, structure_directory, database_fname ):
        struct_fnames = glob.glob( structure_directory + os.sep + "*.pdb" )
        output_f = open( self.output_fname, 'w' )


        for struct_fname in struct_fnames:
            self.run_apsa( struct_fname, self.tmp_fname )
            data = self.parse_apsa_output( struct_fname, self.tmp_fname )
            [output_f.write( ",".join( row ) + "\n" ) for row in data ]

        output_f.close()

        self.add_to_database( self.output_fname, database_fname )




if __name__ == "__main__":
    
    if len( sys.argv ) != 3:
        print """
This runs the apsa secondary structure prediction program on all
structures in the specified directory and puts all the data in one file

To run:
    python add_apsa_sec_struct.py <directory_of_structures>

Output File Format of apsa.output

   <struct_tag>,<sequence position>,<primary structure>,<secondary structure> <database_fname>

Primary Structure categories
XFIRST:  X1
XLAST:  Xn
THREE: _3 
THREEPLUS:  3+
AFIRST: _A 
APLUS:  A+
AMINUS:  A-
HPLUS:  H+
BMINUS:  B-
BPLUS:  B+
JMINUS:  J-
JPLUS:  J+
WMINUS:  W-
WPLUS:  W+
UMINUS:  U-
UPLUS:  U+
NMINUS:  N-
NPLUS:  N+

Secondary Structure labels
HELIXHEAD: <A 
HELIXTAIL:  A>
HELIXREG:  A 
HELIXDIST:  AD
BETAHEAD: <B 
BETABODY:  B 
BETATAIL:  B>
betaHEAD: <b 
betaTAIL:  b>
betaBODY:  b 
TURNHEAD: <T 
TURNBODY:  T 
TURNTAIL:  T>
TURNSINGLE: <T>



S. Raganathen, D. Izotov, E. Kraka, and D. Cremer, "Automated and
accurate protein structure description: Distribution of ideal 
structureal  units in natural proteins", J. Phys. Chem. B, 2008, 
in press.

"""
        exit(0)

    unit = AddAPSA()
    unit.run( sys.argv[1], sys.argv[2] )
