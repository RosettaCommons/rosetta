README

How to use these scripts
========================

This code generates an updated Xcode project from the Xcode project template and Rosetta source files.
Do this whenever moving/adding/deleting files or updating from significant master changes. 

1) cd into this directory
2) python make_project.py all


That should be it! 

Please DO NOT commit changes to the project template unless you know what you are doing. 


Cmake And Classic projects
==========================

Two make projects exist.  The classic is the make_project.py script, which organized Rosetta by library, allows better indexing, searching (shift-cmd-o) and allows Xcode to flip between header and source (ctr-cmd-up). However, currently, updating is a pain.  We are working to automate this.  Please see below. 

The make_project_cmake.py script allows automatic project creation, and a directory-focused code organization.  

make_project_cmake.py 
Author: Andy Watkins 

make_project.py
Author: Unknown, perhaps through FoldIt
Updated/Documented: Jared Adolf-Bryfogle


Static Versions
===============

The rosetta_static and rosetta_scripts_static targets are a proof of principle showing how one could make static binaries using the same object files as those produced by another project. This would allow one to switch all the other libraries to dynamically loaded libriaries, while still being able to compile static executables. The empty rosetta_static.cc file is part of the rosetta_static build.


HOW TO UPDATE XCODE PROJECT WITH NEW LIBRARIES
==============================================

NOTE: This is ONLY for adding new Rosetta libraries - aka protocols.4.src.settings. 
If you add, move, delete files and directories, no change needs to happen for you and re-running make_project.py will 
update xcode with your new code. 

You will need to first add the new library to your Xcode project.  Automation will happen at some point, but is not trivial.

 Steps:
        1) Double Click on Rosetta.  See the Target list.  Right click on the closest library to yours and hit duplicate
        2) Click on the target like you do on Mac finder window and rename.  Drag the target to where it should be.
        3) Double click on your target.  To the right select build phases.  Add your dependencies that are missing to target dependancies
             and link binary to library
        4) Right click sources, hit create group.  Name the group the library.  Drag the group to the appropriate place.
        5) Now, you have to open the src.setting file and find a directory within it.  Right click on your new library in sources.
            Add the file into the group.  Select options (lower left).  Make sure the target is selected and new group is checked
        6) Finally, meander back to targets. Go to build phases.  If you see a compile sources section, click the upper right x to clear it.
           Next, select the + button on the upper left.  Select 'new compile sources' tag.  Under the section, click the + button.
           Select a .cc file corresponding in the directory named the same as this target.
        7) Run make_project.py - make sure everything works.  Copy your .pbxproj file into .pbxproj.template.
        8) Use the grep commands below to get the keys.  Source key is right before it says /* source * /

 group keys grep:
 grep -A 50 'BE8A17150CA8365000D67A6F \/\* Sources \*\/ =' Rosetta.xcodeproj/project.pbxproj
 
 protocols sources keys grep:
 grep -A 2000 '^\/\* Begin PBXNativeTarget' Rosetta.xcodeproj/project.pbxproj | grep -A 6 '\/\* protocols' | grep -E '(\/\* Sources|\/\* protocols)'

