USE_INPUT_POSE

CLAIMER SequenceClaimer
CMD_FLAG #take sequence from -in:file:fasta
LABEL main
END_CLAIMER


CLAIMER CoordConstraintClaimer
#PDB_FILE start_autobuild_4_trim_cterm_0001.pdb
# put root into one of the helices
CST_FROM_INPUT_POSE
#ROOT 33 
POTENTIAL HARMONIC 0 1 
ASK_FOR_ROOT ALL
CST_FILE input/coord.cst
#USE_XYZ_FROM_CSTFILE
END_CLAIMER

CLAIMER StartStructClaimer
#NO_USE_INPUT_POSE
END_CLAIMER

CLAIMER RigidChunkClaimer
# defines a chunk
#NO_USE_INPUT_POSE
REGION_FILE input/core.rigid
END_CLAIMER
