#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http:#www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   torsion_scoring.py
## @brief  no clue
## @author Andrew Watkins


if __name__ == '__main__':
	#devel.init(argc, argv)

	score_fxn = scoring.ScoreFunction()
	score_fxn.set_weight(core.scoring.rna_torsion, 1 )
	residue_name = sys.argv[1]

	pose = Pose()
	for filename in sys.argv[2:]:

		import_pose.pose_from_file( pose, filename, core.import_pose.PDB_file)
		print "Importing pose from ", filename

		Real score = ( *score_fxn )( pose )
		print "Total torsion score: ", score, std.endl
		# TODO: chainbreaks
		for iim1 in xrange(1,pose.size()-1):
			ii = iim1 + 1
			if not pose.residue_type( ii - 1 ).is_RNA(): continue
			if not pose.residue_type( ii ).is_RNA(): continue
			if not pose.residue_type( ii + 1 ).is_RNA(): continue

			# Skip non res
			#TR, "residue_name ", residue_name, " ", pose.residue_type( ii ).name3(), std.endl
			if residue_name != "ALL" && pose.residue_type( ii ).name3() != residue_name: continue

			print ii, " ", pose.residue_type( ii ).name3(), " ",
				pose.residue_type( ii-1 ).name3(), " ",  pose.residue_type( ii+1 ).name3(), # flanking
				" ", pose.torsion( TorsionID( ii, id.BB, BETA  ) ), # intra
				" ", pose.torsion( TorsionID( ii, id.BB, GAMMA ) ), # intra
				" ", pose.torsion( TorsionID( ii, id.BB, DELTA ) ), # intra
				" ", pose.torsion( TorsionID( ii, id.CHI, 1 ) ), # intra, chi1
				" ", pose.torsion( TorsionID( ii, id.CHI, 2 ) ), # intra, nu2
				" ", pose.torsion( TorsionID( ii, id.CHI, 3 ) ), # intra, nu1
				" ", pose.torsion( TorsionID( ii, id.CHI, 4 ) ), # intra, proton chi
				" ", pose.torsion( TorsionID( ii, id.BB, ALPHA ) ), # inter
				" ", pose.torsion( TorsionID( ii - 1, id.BB, ZETA    ) ), # inter
				" ", pose.torsion( TorsionID( ii - 1, id.BB, EPSILON ) ), # inter
				" ", pose.torsion( TorsionID( ii - 1, id.BB, DELTA   ) ), # inter
				#, " ", pose.torsion( TorsionID( ii, id.BB, DELTA ) ) # inter
				" ", pose.torsion( TorsionID( ii, id.BB, EPSILON ) ), # inter
				" ", pose.torsion( TorsionID( ii, id.BB, ZETA ) ), # inter
				" ", pose.torsion( TorsionID( ii + 1, id.BB, ALPHA ) ), # inter
				# LATER: make this just report one term.
				" ", pose.energies().residue_total_energy( ii )


