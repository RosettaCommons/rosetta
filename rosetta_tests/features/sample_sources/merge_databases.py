

import sys, shutil, subprocess, os
import heapq

class MergeDatabases:

    def __init__(self, argv):
        if len(argv) < 3:
            print self.help_str()
            return 1

	self.merge_stmt_prototype = '''
ATTACH DATABASE "%(input_fname)s" AS input_db;
BEGIN TRANSACTION;
%(insert_stmts)s
COMMIT;
'''

        output_fname = argv[1]
        input_fnames = argv[2:]
        self.table_names = self.get_table_names(input_fnames[1])

        self.check_database_parts_are_mergable(input_fnames)

        self.agg_merge(output_fname, input_fnames)
        


    def help_str(self):
        return '''
NAME
     merge_database.py -- Merge multiple SQLite3 database with same schema

SYNOPSIS
     pythos merge_databases.py output input1 [input2 [input3 ...]]

DESCRIPTION
     merge_databases.py puts all data in each input database into output database

     The first argument should be the file name of the output
     database.  The remaining arguments should be the file name of the
     input databases.

WARNING
     This script does very little error detection.  So be careful when
     merging databases that may have different schemas etc.
'''

    def get_table_names(self, fname):
        p = subprocess.Popen(
            args=['sqlite3', fname, "SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%';"],
            stdout=subprocess.PIPE)
        p.wait()
        stdout, stderr = p.communicate()
        table_names = stdout.split("\n")
        return table_names[:-1]

    def apply_sqlite_stmt(self, db_fname, stmt):
        p = subprocess.Popen(
            args=['sqlite3', db_fname, stmt])
        p.wait()

    def check_database_parts_are_mergable(self, input_fnames):
        assert(len(input_fnames) >= 1)
        print "Checking if all database parts have the same tables..."
        self.check_consistent_schema(input_fnames[2:])
        print "Checking that all the structures are identifiable..."
        self.check_structures_are_identifiable(input_fnames) 

    def check_consistent_schema(self, input_fnames):
        for input_fname in input_fnames:
            if self.table_names != self.get_table_names(input_fname):
                schema1 = "\n\t".join(self.table_names)
                schema2 = "\n\t".join(self.get_table_names(input_fname))
                print """ERROR: Database parts have the different tables:
ERROR: Database schema for %s:
ERROR: 	%s
ERROR: 
ERROR: Database schema for %s:
ERROR: 	%s
ERROR:
ERROR: This may mean that something has gone wrong. To draw your attention to
ERROR: this fact, this script is exiting. However, if you know what you are
ERROR: doing, then feel free to remove this check and continue on your way.
""" % (input_fnames[1], schema1, input_fname, schema2)
                
                exit(1)

        
    def check_structures_are_identifiable(self, input_fnames):
        # Check that in the different database parts:
        #   No two structures have the same struct_id
        #   No two structures have the same (protocol_id, tag) 


        # fill these dict objects with info about the structures found so far
        struct_ids = {}
        protocol_tags = {}

        for input_fname in input_fnames:
            p = subprocess.Popen(
                args=['sqlite3', input_fname, "SELECT protocol_id, struct_id, tag FROM structures;"],
            	stdout=subprocess.PIPE)
            p.wait()
            stdout, stderr = p.communicate()
            rows = stdout.split("\n")
            for row in rows[:-1]:
                protocol_id, struct_id, tag = row.split("|")
                if struct_id in struct_ids:
                    print """ERROR: The database parts have duplicate 'struct_id' field in the structures table
ERROR: 
ERROR:     database part '%s' and '%s' both have struct_id '%s'.
ERROR: 
ERROR: When generating feature database parts that are to be merged using the
ERROR: ReportToDB mover in the RosettaScripts and not using MPI, you
ERROR: need to manually specify the initial starting struct_id. For
ERROR: example If you have 2000 structures split into 100 runs
ERROR: (indexed from 0 to 99 by run_index in your script) of 20
ERROR: structures each set the first_struct_id field like this:
ERROR: 
ERROR:     <ReportToDB ... first_struct_id=run_index*20 ...>
""" % (struct_ids[struct_id], input_fname, struct_id)
                    exit(1)

                else:
                    struct_ids[struct_id] = input_fname

                if (protocol_id, tag) in protocol_tags:
                    print """ERROR: The database parts have duplicate "(protocol_id, tag) fields in the structures table

ERROR:   database part '%s' and '%s' both have protocol_id: '%s' and struct_id: '%s' 
ERROR:
ERROR: When generating feature database parts that are to be merged using the
ERROR: ReportToDB mover in the RosettaScripts where the same tag may be used
ERROR: more than once, you should specify the protocol_id field to keep the
ERROR: structures separate:
ERROR: 
ERROR:     <ReportToDB ... protocol_id=(%string) ...>
""" % (protocol_tags[(protocol_id, tag)], input_fname, protocol_id, tag)
                    exit(1)
                else:
                    protocol_tags[(protocol_id, tag)] = input_fname



    def merge_database(self, source_fname, target_fname):
        for table_name in self.table_names:
            insert_stmt = "INSERT OR IGNORE INTO %(table_name)s SELECT * FROM input_db.%(table_name)s;"
            insert_stmt = insert_stmt % {"table_name" : table_name}

            stmt = self.merge_stmt_prototype % {
                "input_fname" : source_fname,
                "insert_stmts" : insert_stmt}

            self.apply_sqlite_stmt(target_fname, stmt)
        
    def get_nstruct(self, database_fname):
        p = subprocess.Popen(
            args=['sqlite3', database_fname, "SELECT count(*) FROM structures;"],
            stdout=subprocess.PIPE)
        p.wait()
        stdout, stderr = p.communicate()
        print "nstruct:", stdout
        return int(stdout[:-1])


    def agg_merge(self, output_fname, input_fnames):
        # For performance, merging smaller databases seems to be faster

        clusters = []
        [heapq.heappush(clusters, (self.get_nstruct(input_fname), input_fname)) for input_fname in input_fnames]
        for i in range(len(input_fnames) - 1):
            n1, input_fname1 = heapq.heappop(clusters)
            n2, input_fname2 = heapq.heappop(clusters)
            print "merging '%s' -> '%s', giving %s structs. %s dbs remaining..."%\
                (input_fname1, input_fname2, n1 + n2, len(input_fnames)-1-i)
            self.merge_database(input_fname1, input_fname2)
            heapq.heappush(clusters, (n1 + n2, input_fname2))
            os.remove(input_fname1)
        nstruct, last_fname = heapq.heappop(clusters)
        print "Moving last db '%s' to output db '%s'." % (last_fname, output_fname)
        shutil.move(last_fname, output_fname)

            
if __name__ == '__main__':
    MergeDatabases(sys.argv)
