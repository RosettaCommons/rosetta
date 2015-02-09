// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_make_rot_lib_makerotlib_HH
#define INCLUDED_protocols_make_rot_lib_makerotlib_HH

// core headers
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// unit headers
#include <protocols/make_rot_lib/RotData.fwd.hh>

// c++ headers
// AUTO-REMOVED #include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace make_rot_lib {

// hack hack hack
void
asp_corrections( RotVec & );

void
glu_corrections( RotVec & );

void
phe_tyr_corrections( RotVec & );

// minimize side chain dihedral angles of each rotamer
void
min_rotamers (RotVec &, core::scoring::ScoreFunctionOP, std::string);

// fill in centroid with cmd-line opts
void
init_rotamers_centroids(RotVec &, RotVec &, core::Size &, std::string, std::string &, bool, core::Real, core::Real );

// itterates over distances to centroids, assigns to closest centroid
// assigns rots to clusters and determins if clusters changed
bool
calc_rotamer_clusters (RotVec &);

// find new centroids
bool
calc_centroids (RotVec &, RotVec &);

//calculates distance b/w 2 RotData objects
core::Real
calc_dist (RotData &, RotData &);

//calculates avg distance between centroid and all points in its cluster
core::Real
avg_cluster_cen_dist (RotVec &, core::Size &);

//finds all distances between all rotamers and all centroids
void
calc_all_dist (RotVec &, RotVec &);

// pull out best rots from rotamers and add them to final_rotamers
void
get_final_rots(RotVec &, RotVec &, core::Size &);

// calc probabilites for final rots & normilize
void
get_final_rot_probs( RotVec &);

void
calc_std_dev (RotVec &, core::scoring::ScoreFunctionOP, std::string );

// prints out rotamers
void
pretty_print_rd( RotData & );

// print out rotamers acording to the Dunbrack format
void
dunbrack_print( RotVec &, RotVec &, std::string );

} // namespace make_rot_lib
} // namespace protocols

#endif // INCLUDED_protocols_makerotlib_makerotlib_HH

