Setup:

1) Make sure paths for PyRosetta are set; Version = 2.1

2) Run program,
	-cd into Toolkit directory and run using "Python pyrosetta_toolkit.py"
	-If sqlite3 does not work, reinstall python




NOTE FOR DEVELOPERS:  More thorough information at: https://wiki.rosettacommons.org/index.php/PyRosetta_Toolkit



Usage:

Many things are available to do within the program.  You can design, model, analyze, and work with pdb files.
Many menus are available to complete these tasks.  To start modeling, load a PDB by clicking the find button.
This loads the pose into rosetta.  Once there, you can quickly do minimization and get the sequence of regions.
run protocols using the right side of the screen, setup design by launching the Design File Toolbox (Located in the upper left of the menu) or look at energies, rmsd, and backbone Phi/Psi by selecting the analysis options on the right of the screen.  Each instance of the program runs on one processor; so if you have multiple protocols you need to run and multiple processors, launch as many instances you have processors!


PyMol Integration:

The ability to output any modeling step and your pose at any time is provided.  This is built into the Rosetta code, thanks to the Gray Lab.  First - Open Pymol and cd into your PyRosetta directory.  Next type "run PyMOLPyRosettaServer.py".  You can click "Show In Pymol" to see your pose at any point after it is returned to you.  If you click the "Observe in Pymol?" checkbox, you will be watching each protocol that is available to be watched.  If the protocol has no intrinsic pymol output, the pose is output each time the score of the pose is checked.
Note that this severly limits the speed of the protocol, so only use sparingly.  Especially the Relax protocol.  Though, it is extremely useful to quickly understand what the protocol is doing.


Scores:

Scores can be chosen and configured using the Score Menu.  This gives you access to all of the in-built Rosetta Scores, as well as the ability to create new ones.  


Regions:

Loops can be added by specifying the Start and End of the loop in the left part of the main screen, as well as the chain and pushing the "Add Loop" button.  Any number of loops or regions can be created.  If you would like to work on a specific chain, simply specify the chain and nothing else.  Most protocols support individual chain modeling.  If you would like to work on the N or C term, specify where the tail begins or ends.  Example: You would like to model from residue 1 to residue 15 on chain H.  Here, you would only specify 15 as 'end of loop' and chain H.  This will do Rosettas Tail Modeling if possible, otherwise it will write a movemap to minimize or disregard any movement from other atoms in your protein.

Design:

Use the Design file toolbox to build a resfile.  Each residue needs to be specified.  You can do multiple residues by specifying "start:end" instead of one residue. 

FullControl:

FullControl Toolbox allows you to check the energies of any residue in the pose, make point mutations, and change Phi/Psi manually.


PDB Tools:

Some useful PDB tools are provided.  These include: Chain Character change, Occupancy change to 1.0, and fixing pdb files with iCodes for residue numbers.  You may also use the Meiler label_ Clean_PDB.py script that is floating around.  More useful tools will follow.

Options:

Default PYRosetta options can be configured using the menu.



Rosetta Command-line creator:

For users to create a config file and look at all available options for an app.



Please See Wiki for more details.

Notes:

Rosetta requires all occupancies to be 1.0, or they will not be loaded into the pose.
This can be fixed by specifying a filename, selecting the pdb menu, loading the pdb, and then selecting 'change occupancies to 1.0'.  

If you have names of atoms or residues that are uncommon (Like that from MD, you will need to clean the pdb option to clean up your pdb before use in PyRosetta.  You can do this through the menu, or through the Meiler label_'s 'clean_pdb' script.  Clicking the menu option allows you to save the new file any where you would like and name it.







Created by Jared Adolf-Bryfogle, see the wiki for more details.


