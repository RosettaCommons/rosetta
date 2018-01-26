#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http:#www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   test_mm_lj.py
## @brief  ??
## @author Andrew Watkins

def examine_interval( i, j, path_distance, min, max, step, mmljet ):

	for d in xrange(min, max, step):

		# arbitrary factor, just so we don't spend SO much of our time zooming in
		Size mult = 2
		Real rep1, atr1, drep1, datr1
		Real rep2, atr2, drep2, datr2

		mmljet.score(i, j, path_distance, d, rep1, atr1)
		mmljet.deriv_score(i, j, path_distance, d, drep1, datr1)

		Real d2 = d + step
		mmljet.score(i, j, path_distance, d2, rep2, atr2)
		mmljet.deriv_score(i, j, path_distance, d2, drep2, datr2)

		Real ndrep = (rep2 - rep1) / step
		Real ndatr = (atr2 - atr1) / step
		Real drep = (drep1 + drep2) / 2
		Real datr = (datr1 + datr2) / 2
		if ndrep - drep > step * mult or ndrep - drep < -1 * step * mult:
			print "Confirming possible issue (", ndrep, "vs", drep, ") between", d, d2, "at", step, "level."
			examine_interval(i, j, path_distance, d - step / mult, d + step, step / mult, mmljet)

		if ndatr - datr > step*mult or ndatr - datr < -1*step*mult:
			print "Confirming possible issue (", ndrep, "vs", drep, ") between ", d, d2, "at", step, "level."
			examine_interval(i, j, path_distance, d - step / mult, d + step, step / mult, mmljet)


if __name__ == '__main__':
	#devel.init(argc, argv)

	mmljet = core.scoring.mm.MMLJEnergyTable()
	mmljs = core.scoring.mm.MMLJScore()

	for i in xrange(1, 104+1):
		for j in xrange(1, 104+1):
			for path_distance in [3, 4]:
				Real step = 1e-8
				#Real orig_step = 1e-8
				#Real orig_d
				#std.cout, "Testing ", i, " ", j, " ", path_distance, " switch at ", ( 0.6*mmljs.min_dist( i, j, path_distance ) ), " min at ", mmljs.min_dist( i, j, path_distance ), std.endl
				#std.cout, "(in dist_squared that's ", (0.6*mmljs.min_dist( i, j, path_distance )*0.6*mmljs.min_dist( i, j, path_distance )), " and  ", mmljs.min_dist( i, j, path_distance )*mmljs.min_dist( i, j, path_distance )<< std.endl
				for d in xrange(mmljs.min_dist(i, j, path_distance) * mmljs.min_dist(i, j, path_distance) * .6, 64, step*1e5):
					rep1, atr1, drep1, datr1 = 0, 0, 0, 0
					rep2, atr2, drep2, datr2 = 0, 0, 0, 0
					
					mmljet.score(i, j, path_distance, d, rep1, atr1)
					mmljet.score(i, j, path_distance, d+step, rep2, atr2)
					mmljet.deriv_score(i, j, path_distance, d, drep1, datr1)
					mmljet.deriv_score(i, j, path_distance, d+step, drep2, datr2)
					ndrep = (rep2 - rep1) / step
					ndatr = (atr2 - atr1) / step
					drep = (drep1 + drep2) / 2
					datr = (datr1 + datr2) / 2
					if ndrep - drep > 1e-6 or ndrep - drep < -1e-6 or ndatr - datr > 1e-6 or ndatr - datr < -1e-6:
						if ndrep - drep > 1e-6 or ndrep - drep < -1e-6:
							print i, " ", j, " rep", "\t", mmljs.min_dist(i, j, path_distance), "\t", d, "\t", ndrep, "\t", drep
						if ndatr - datr > 1e-6 or ndatr - datr < -1e-6:
							print i, " ", j, " atr", "\t", mmljs.min_dist(i, j, path_distance), "\t", d, "\t", ndatr, "\t", datr
					elif ndrep - drep > 1e-7 or ndrep - drep < -1e-7 or ndatr - datr > 1e-7 or ndatr - datr < -1e-7:
						mmljet.score(i, j, path_distance, d + step / 10, rep2, atr2)
						mmljet.deriv_score(i, j, path_distance, d + step / 10, drep2, datr2)
						ndrep = (rep2 - rep1) / step * 10
						ndatr = (atr2 - atr1) / step * 10
						drep = (drep1 + drep2) / 2
						datr = (datr1 + datr2) / 2
						if ndrep - drep > 1e-7 or ndrep - drep < -1e-7:
							print i, " ", j, " rep", "\t", mmljs.min_dist(i, j, path_distance), "\t", d, "\t", ndrep, "\t", drep
						if ndatr - datr > 1e-7 or ndatr - datr < -1e-7:
							print i, " ", j, " atr", "\t", mmljs.min_dist(i, j, path_distance), "\t", d, "\t", ndatr, "\t", datr
						

					#examine_interval( i, j, path_distance, d, d+step, step, mmljet )

