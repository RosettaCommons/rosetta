# database_md5

This is a simple test which collects the MD5 signatures of all the files in the database.

This is mainly to be a diagnostic aid in debugging odd test errors on the test server,
by leting us know if there's some sort of undesired state difference between different backend
servers on the test server.

If you expect the contents of a database file to change with your commit, then don't worry about
the test diffs here. 
