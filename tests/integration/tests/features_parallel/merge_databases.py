

import sys, shutil, subprocess

class MergeDatabases:

    def __init__(self, argv):
        if len(argv) < 3:
            print self.help_str()
        
        self.merge_databases(argv[1], argv[2:])
        

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

    def build_copy_stmt(self, input_fname, table_names):
        stmt = '''
PRAGMA foreign_keys=ON;
ATTACH DATABASE "%(input_fname)s" AS input_db;
BEGIN TRANSACTION;
%(insert_stmts)s
COMMIT;
'''
        insert_stmts = []
        for table_name in table_names:
            insert_stmt = "INSERT INTO %(table_name)s SELECT * FROM input_db.%(table_name)s;"
            insert_stmt = insert_stmt % {"table_name" : table_name}
            insert_stmts.append(insert_stmt)

        stmt = stmt % {
            "input_fname" : input_fname,
            "insert_stmts" : "\n".join(insert_stmts)}
        return stmt

    def apply_sqlite_stmt(self, db_fname, stmt):
        p = subprocess.Popen(
            args=['sqlite3', db_fname, stmt])
        p.wait()


    def merge_databases(self, output_fname, input_fnames):
        shutil.copy(input_fnames[0], output_fname)
        table_names = self.get_table_names(output_fname)

        for input_fname in input_fnames[1:]:
            stmt = self.build_copy_stmt(input_fname, table_names)
            self.apply_sqlite_stmt(output_fname, stmt)
            
if __name__ == '__main__':
    MergeDatabases(sys.argv)
