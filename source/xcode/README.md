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
  1. Click on Rosetta. Select from the 'Targets' list the closest library to yours, right click and hit duplicate.
  2. Click on your newly duplicated target and rename it to the match the new library. Drag the new target in the list to where it should be.
  3. Click on your newly duplicated target and look at the top of xcode for four tabs. Click the third "Build Phases" tab.
    * In the "Dependencies" panel, add missing dependencies (e.g. if you just created core.6 from core.5, add core.5).
      * Repeat this dependencies process for any higher level libraries.
    * In the "Link Binary With Libraries" link the appropriate binaries to the new library.
      * Repeat this linking process for any higher level libraries.
    * Delete the "Compile Sources" by clicking the "x".
  4. Now to the folder navigation on the right panel. Right click sources, hit "Create Group" (don't worry about the newly create folder). Rename the group with the library name. Drag the group to the appropriate place.
  5. Outside of Xcode, find the appropriate src.settings file in Rosetta/main/source/src (e.g. core.6.src.settings). From that file, select a *.cc file and note its path.
  6. Back to Xcode, Right click the newly created group under "Sources" and select "Add files to "Rosetta"". Then add the noted file to the group. When adding, make sure the corret target is selected and unchecked any additional targets.
  7. Select the Rosetta again and navigate to "Build Phases". Click "+" and add a "New Compile Sources Phase". In the newly added section, click "+" and add your *.cc file.
  8. Copy your .pbxproj file into .pbxproj.template. This is important this will add all your hard work creating the new library to the template!
  9. Grep some identifying values. The below commands did not work for me, but I looked in the corresponding sections of Rosetta.xcodeproj/project.pbxproj for my library (core.6).
  10. Use the grep commands below to get the keys.  Source key is right before it says /* source * /

 group keys grep:
 grep -A 50 'BE8A17150CA8365000D67A6F \/\* Sources \*\/ =' Rosetta.xcodeproj/project.pbxproj

 protocols sources keys grep:
 grep -A 2000 '^\/\* Begin PBXNativeTarget' Rosetta.xcodeproj/project.pbxproj | grep -A 6 '\/\* protocols' | grep -E '(\/\* Sources|\/\* protocols)'
