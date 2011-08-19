CLAIMER SequenceClaimer
#either FILE or DEF.. END_DEF
LABEL DEFAULT
FILE input/1dsvA.fasta
END_CLAIMER

CLAIMER MetalloClaimer
ligand ZN 
anchor 5
aa CYS ZN
cst_atoms SG CB CA ZN V1 V2
jump_atoms N CA C ZN V1 V2
disAB 2.20
angleA 68.0
angleB 70.5
dihedralA -150.0 -120.0 -90.0 -60.0 -30.0 0.0 30.0 60.0 90.0 120.0 150.0 180.0
dihedralAB -150.0 -120.0 -90.0 -60.0 -30.0 0.0 30.0 60.0 90.0 120.0 150.0 180.0
dihedralB 120.0
END_CLAIMER

#CLAIMER ConstraintClaimer
#FILE 
#FILE input/t286_.pdb.distances.csts.harm
#END_CLAIMER

ABINITIO_FRAGS
#this will initialize fragments mover as in classic-abinitio
#small frags and large frags   smooth moves in stage4
LARGE input/aa1dsvA09_05.200_v1_3
SMALL input/aa1dsvA03_05.200_v1_3
END_ABINITO

