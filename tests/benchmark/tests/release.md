# release test suite
This is a special test suite. In contrast to other suites each sub-test in this suite is designed to create various Rosetta release packages during run.
If run is successful then test package will became available for download in Rosetta Release system.

-----
### source
Create Rosetta **source** release package. On test completion package will became available from http://graylab.jhu.edu/download/

-----
### binary
Create Rosetta binary release package. On test completion package will became available from http://graylab.jhu.edu/download/

-----
### PyRosetta4 Debug, Release. MinSizeRel, RelWithDebInfo
Create corresponding PyRosetta4 release package. On test completion package will became available from http://graylab.jhu.edu/download/

-----
### PyRosetta4.debug
Build PyRosetta4 and generate API documentaton for it. If this is successful upload new documentation so it became available from http://graylab.jhu.edu/PyRosetta.documentation/
