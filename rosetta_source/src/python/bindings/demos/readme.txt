this folder is intended to hold demo scripts

as they are currently, these scripts are in test/
they live there because they are useful for testing PyRosetta functionality

we want to move those scripts here, mainly since users often ignore test/
even though demonstrative scripts live there

after the scripts are moved, make sure their default input and output
are in the test/ directory AND that every python script in the demos/
directory is run during the testing process (i.e. test ALL in demos/ and test/
directories)

please move the following scripts from test/ to demos/

test/D010_Pose_structure.py
test/D020_Pose_scoring.py
test/D030_Fold_tree.py
test/D040_Movemap.py
test/D050_Packer_task.py
test/D060_Folding.py
test/D070_Refinement.py
test/D080_Loop_modeling.py
test/D090_Ala_scan.py
test/D100_Docking.py
test/D110_DNA_interface.py
test/D120_Ligand_interface.py


Thanks!
