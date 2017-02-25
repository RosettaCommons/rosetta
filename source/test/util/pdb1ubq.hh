// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/util/pdb1rpb.hh
/// @brief  Create a nine-residue pose from a fragment of the ubiquitin structure.
/// @author Andrew Leaver-Fay

#ifndef INCLUDED_util_pdb1ubq_HH
#define INCLUDED_util_pdb1ubq_HH

#include <test/util/pose_funcs.hh>

std::string pdb1ubq5to13();

inline
core::pose::Pose
pdb1ubq5to13_pose()
{
	return fullatom_pose_from_string( pdb1ubq5to13() );
}

inline
core::pose::PoseOP
pdb1ubq5to13_poseop()
{
	return fullatom_poseop_from_string( pdb1ubq5to13() );
}

inline
std::string
pdb1ubq5to13()
{
	return
		"ATOM     37  N   VAL     5      28.260  33.943  11.096  1.00  4.44           N  \n"
		"ATOM     38  CA  VAL     5      28.605  33.965  12.503  1.00  3.87           C  \n"
		"ATOM     39  C   VAL     5      28.638  35.461  12.900  1.00  4.93           C  \n"
		"ATOM     40  O   VAL     5      29.522  36.103  12.320  1.00  6.84           O  \n"
		"ATOM     41  CB  VAL     5      29.963  33.317  12.814  1.00  2.99           C  \n"
		"ATOM     42  CG1 VAL     5      30.211  33.394  14.304  1.00  5.28           C  \n"
		"ATOM     43  CG2 VAL     5      29.957  31.838  12.352  1.00  9.13           C  \n"
		"ATOM     44  N   LYS     6      27.751  35.867  13.740  1.00  6.04           N  \n"
		"ATOM     45  CA  LYS     6      27.691  37.315  14.143  1.00  6.12           C  \n"
		"ATOM     46  C   LYS     6      28.469  37.475  15.420  1.00  6.57           C  \n"
		"ATOM     47  O   LYS     6      28.213  36.753  16.411  1.00  5.76           O  \n"
		"ATOM     48  CB  LYS     6      26.219  37.684  14.307  1.00  7.45           C  \n"
		"ATOM     49  CG  LYS     6      25.884  39.139  14.615  1.00 11.12           C  \n"
		"ATOM     50  CD  LYS     6      24.348  39.296  14.642  1.00 14.54           C  \n"
		"ATOM     51  CE  LYS     6      23.865  40.723  14.749  1.00 18.84           C  \n"
		"ATOM     52  NZ  LYS     6      22.375  40.720  14.907  1.00 20.55           N  \n"
		"ATOM     53  N   THR     7      29.426  38.430  15.446  1.00  7.41           N  \n"
		"ATOM     54  CA  THR     7      30.225  38.643  16.662  1.00  7.48           C  \n"
		"ATOM     55  C   THR     7      29.664  39.839  17.434  1.00  8.75           C  \n"
		"ATOM     56  O   THR     7      28.850  40.565  16.859  1.00  8.58           O  \n"
		"ATOM     57  CB  THR     7      31.744  38.879  16.299  1.00  9.61           C  \n"
		"ATOM     58  OG1 THR     7      31.737  40.257  15.824  1.00 11.78           O  \n"
		"ATOM     59  CG2 THR     7      32.260  37.969  15.171  1.00  9.17           C  \n"
		"ATOM     60  N   LEU     8      30.132  40.069  18.642  1.00  9.84           N  \n"
		"ATOM     61  CA  LEU     8      29.607  41.180  19.467  1.00 14.15           C  \n"
		"ATOM     62  C   LEU     8      30.075  42.538  18.984  1.00 17.37           C  \n"
		"ATOM     63  O   LEU     8      29.586  43.570  19.483  1.00 17.01           O  \n"
		"ATOM     64  CB  LEU     8      29.919  40.890  20.938  1.00 16.63           C  \n"
		"ATOM     65  CG  LEU     8      29.183  39.722  21.581  1.00 18.88           C  \n"
		"ATOM     66  CD1 LEU     8      29.308  39.750  23.095  1.00 19.31           C  \n"
		"ATOM     67  CD2 LEU     8      27.700  39.721  21.228  1.00 18.59           C  \n"
		"ATOM     68  N   THR     9      30.991  42.571  17.998  1.00 18.33           N  \n"
		"ATOM     69  CA  THR     9      31.422  43.940  17.553  1.00 19.24           C  \n"
		"ATOM     70  C   THR     9      30.755  44.351  16.277  1.00 19.48           C  \n"
		"ATOM     71  O   THR     9      31.207  45.268  15.566  1.00 23.14           O  \n"
		"ATOM     72  CB  THR     9      32.979  43.918  17.445  1.00 18.97           C  \n"
		"ATOM     73  OG1 THR     9      33.174  43.067  16.265  1.00 20.24           O  \n"
		"ATOM     74  CG2 THR     9      33.657  43.319  18.672  1.00 19.70           C  \n"
		"ATOM     75  N   GLY    10      29.721  43.673  15.885  1.00 19.43           N  \n"
		"ATOM     76  CA  GLY    10      28.978  43.960  14.678  1.00 18.74           C  \n"
		"ATOM     77  C   GLY    10      29.604  43.507  13.393  1.00 17.62           C  \n"
		"ATOM     78  O   GLY    10      29.219  43.981  12.301  1.00 19.74           O  \n"
		"ATOM     79  N   LYS    11      30.563  42.623  13.495  1.00 13.56           N  \n"
		"ATOM     80  CA  LYS    11      31.191  42.012  12.331  1.00 11.91           C  \n"
		"ATOM     81  C   LYS    11      30.459  40.666  12.130  1.00 10.18           C  \n"
		"ATOM     82  O   LYS    11      30.253  39.991  13.133  1.00  9.10           O  \n"
		"ATOM     83  CB  LYS    11      32.672  41.717  12.505  1.00 13.43           C  \n"
		"ATOM     84  CG  LYS    11      33.280  41.086  11.227  1.00 16.69           C  \n"
		"ATOM     85  CD  LYS    11      34.762  40.799  11.470  1.00 17.92           C  \n"
		"ATOM     86  CE  LYS    11      35.614  40.847  10.240  1.00 20.81           C  \n"
		"ATOM     87  NZ  LYS    11      35.100  40.073   9.101  1.00 21.93           N  \n"
		"ATOM     88  N   THR    12      30.163  40.338  10.886  1.00  9.63           N  \n"
		"ATOM     89  CA  THR    12      29.542  39.020  10.653  1.00  9.85           C  \n"
		"ATOM     90  C   THR    12      30.494  38.261   9.729  1.00 11.66           C  \n"
		"ATOM     91  O   THR    12      30.849  38.850   8.706  1.00 12.33           O  \n"
		"ATOM     92  CB  THR    12      28.113  39.049  10.015  1.00 10.85           C  \n"
		"ATOM     93  OG1 THR    12      27.280  39.722  10.996  1.00 10.91           O  \n"
		"ATOM     94  CG2 THR    12      27.588  37.635   9.715  1.00  9.63           C  \n"
		"ATOM     95  N   ILE    13      30.795  37.015  10.095  1.00 10.42           N  \n"
		"ATOM     96  CA  ILE    13      31.720  36.289   9.176  1.00 11.84           C  \n"
		"ATOM     97  C   ILE    13      30.955  35.211   8.459  1.00 10.55           C  \n"
		"ATOM     98  O   ILE    13      30.025  34.618   9.040  1.00 11.92           O  \n"
		"ATOM     99  CB  ILE    13      32.995  35.883   9.934  1.00 14.86           C  \n"
		"ATOM    100  CG1 ILE    13      33.306  34.381   9.840  1.00 14.87           C  \n"
		"ATOM    101  CG2 ILE    13      33.109  36.381  11.435  1.00 17.08           C  \n"
		"ATOM    102  CD1 ILE    13      34.535  34.028  10.720  1.00 16.46           C  \n";
}

#endif
