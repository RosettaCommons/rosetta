USE_INPUT_POSE

CLAIMER SequenceClaimer
CMD_FLAG #take sequence from -in:file:fasta
END_CLAIMER


#CLAIMER CoordConstraintClaimer
#PDB_FILE start_autobuild_4_trim_cterm_0001.pdb
# put root into one of the helices
#CST_FROM_INPUT_POSE
#ROOT 33 
#CST_FILE coord_in.cst

#REGION 
#LOOP 1 42 0 0 0
#END_REGION
#POTENTIAL BOUNDED 0.0 4 1 xyz
#END_CLAIMER

#CLAIMER StartStructClaimer
#NO_USE_INPUT_POSE
#FILE start_autobuild_4_trim_cterm_0001.pdb
#END_CLAIMER

CLAIMER RigidChunkClaimer
# defines a chunk
USE_THREADING_LOOPS
#default is 0 
#MIN_LOOP_SIZE 5
END_CLAIMER
