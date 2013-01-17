APBS Release Procedures
-----------------------

 * Change Version Number
   - Navigate to the apbs root folder
     $ cd ~/apbs
   - Edit CMakeLists.txt
   - Go to the "Set up basic project stuff" section (around line 45)
   - Change the value for the APBS_VERSION variable:
     set(APBS_VERSION "X.X")
     
     
     
 * Update the ChangeLog
   - Navigate to the apbs root folder
     $ cd ~/apbs
   - Edit ChangeLog
   - Docmument major changes for this release
   
   
   
 * Update License info
   - Update license dates and information in source files
   - In apbs/src edit all .c source files and all .h header files, update dates
   
   
   
 * Test release
 
   - Set up seperate machines or virtual machines for target deploy platforms:
     Ubuntu x86_64
     Ubuntu i386
     Redhat
     Mac OSX
     Windows 7 x86_64
     Windows 7 i386
     Windows 7 with Cygwin
     Windows 7 with Mingw
     
   - On testing platforms install or verify presence of required tools:
     Essential compile toolchain 
     python
     git
     CMake
     Doxygen
     LaTeX builder like texlive, tetex (usually already available in linux)
     
   - Clone apbs git repository from sourceforge to testing machines:
     $ mkdir ~/apbs
     $ cd ~/apbs
     $ git init .
     $ git remote add apbs \
     > ssh://tuckerbeck@apbs.git.sourceforge.net/gitroot/apbs/apbs 
     $ git pull apbs master
       
   - Build testing
     > On each platform test machine
     > Navigate to apbs root folder
       $ cd ~/apbs
     > Run CMake selecting an install target with open permissions
       $ cmake -DCMAKE_INSTALL_PREFIX=~/apbs-install
     > Make and install to target apbs location
       $ make install
     > Ensure that build progresses without any errors.
     
   - Run testing
     > Navigate to apbs testing folder
       $ cd ~/apbs/tests
     > run standard testing suite on test platforms
       $ python apbs_tester.py
     > Ensure that all tests run without seg fault and results are acceptable
     
     
     
 * Build and upload Binary Packages
 
   - On the following test plotforms:
     Ubuntu x86_64
     Ubuntu i386
     Mac OSX
     Windows 7 x86_64
     Windows 7 i386
     
   - Add entire install structure to archive file named in the following way
     APBS-X.X-PLATFORM-ARCHITECTURE.ARCHIVE_EXTENSION
     $ tar -czvf APBS-1.4-linux-x86_64.tar.gz ~/apbs-install
     
   - Upload the archive to apbs project on sourceforge
     https://sourceforge.net/projects/apbs/files/apbs/
     
     
     
 * Upload Source Package
 
   - On development machine, Navigate to apbs directory
     $ cd ~/apbs
   
   - Use git to remove all non-versioned files and directories
     Use -dfn flags first for a dry run and make sure right files are rm'ed
     $ git clean -dfq
     
   - Add the whole apbs directory to an arcive
     $ tar -czvf APBS-1.4-source.tar.gz ~/apbs/
     
   - Upload the archive to apbs project on sourceforge
     https://sourceforge.net/projects/apbs/files/apbs/
     
     
     
 * Upload Package Programmer Guide Package
 
   - Build documentation
     > Navigate to apbs build folder
       $ cd ~/apbs/build
     > Run CMake with documentation building option
       $ cmake .. -DBUILD_DOC=YES
     > Build everything
       $ make
     > LaTeX and html output are generated in apbs doc folder
     
   - Navigate to the apbs documentation guide folder
     $ cd ~/apbs/doc
     
   - Add the whole programmer guide directory to an archive
     $ tar -czvf APBS-1.4-programmer_guide.tar.gz
     
   - Upload the archive to apbs project on sourceforge
     https://sourceforge.net/projects/apbs/files/apbs/

 * Update http://www.poissonboltzmann.org/apbs/release-history with new release information.
