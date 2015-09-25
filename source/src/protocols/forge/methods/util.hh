// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/forge/methods/util.hh
/// @brief  miscellaneous utility functions for forge
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_protocols_forge_methods_util_hh
#define INCLUDED_protocols_forge_methods_util_hh

// type headers
#include <core/types.hh>

// package headers
#include <protocols/forge/build/Interval.fwd.hh>

// project headers
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/graph/DisjointSets.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>

// C++ headers
#include <set>
#include <list>
#include <iostream>

#include <utility/vector1.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>


namespace protocols {
namespace forge {
namespace methods {


/// @brief return a set containing values in the closed interval [left, right]
///  filled by looping from left -> right with the given increment
template< typename T >
std::set< T >
closed_range(
	T const left,
	T const right,
	T const increment = 1
)
{
	std::set< T > s;

	for ( T i = left; i <= right; i += increment ) {
		s.insert( i );
	}

	return s;
}


/// @brief add the values in the closed interval [left, right] to a set
///  by looping from left -> right with the given increment
template< typename T >
void
insert_closed_range(
	T const left,
	T const right,
	std::set< T > & s,
	T const increment = 1
)
{
	for ( T i = left; i <= right; i += increment ) {
		s.insert( i );
	}
}


/// @brief return a set containing values in the half-open interval [left, right)
///  filled by looping from left -> right with the given increment
template< typename T >
std::set< T >
half_open_range(
	T const left,
	T const right,
	T const increment = 1
)
{
	std::set< T > s;

	for ( T i = left; i < right; i += increment ) {
		s.insert( i );
	}

	return s;
}


/// @brief add the values in the half-open interval [left, right) to a set
///  by looping from left -> right with the given increment
template< typename T >
void
insert_half_open_range(
	T const left,
	T const right,
	std::set< T > & s,
	T const increment = 1
)
{
	for ( T i = left; i < right; i += increment ) {
		s.insert( i );
	}
}


/// @brief perform union( root, i ) for all 'i' within the closed interval
///  [left, right]
/// @param[in] root position to union with; can be any number, does not have
///  to be a true root of a set
/// @param[in] left start of the interval
/// @param[in] right end of the interval
/// @param[in,out] uf
void
union_interval(
	core::Size const root,
	core::Size const left,
	core::Size const right,
	core::graph::DisjointSets & uf
);


/// @brief moving left to right, find the first true cutpoint within the specified
///  extent
/// @return the cutpoint position, otherwise 0 if not found
core::Size
find_cutpoint(
	core::pose::Pose const & pose,
	core::Size left,
	core::Size right
);


/// @brief moving left to right, count the number of true cutpoints within the
///  specified extent
core::Size
count_cutpoints(
	core::pose::Pose const & pose,
	core::Size left,
	core::Size right
);


/// @brief set omega to 180 for a range of residues [left, right]
void
trans_omega(
	core::Size const left,
	core::Size const right,
	core::pose::Pose & pose
);


/// @brief create Loop object w/ random cutpoint from an Interval
protocols::loops::Loop
interval_to_loop( protocols::forge::build::Interval const & interval );


/// @brief create Loops object w/ random cutpoints from a collection of
///  Intervals
template< typename IntervalIterator >
protocols::loops::Loops
intervals_to_loops(
	IntervalIterator begin,
	IntervalIterator end
)
{
	using protocols::loops::Loop;
	using protocols::loops::Loops;

	Loops loops;

	for ( IntervalIterator i = begin; i != end; ++i ) {
		loops.add_loop( interval_to_loop( *i ) );
	}

	return loops;
}

/// @brief create Loops object w/ random cutpoints from a collection of
///  Intervals for KIC confirmation purpose
template< typename IntervalIterator >
protocols::loops::Loops
intervals_to_confirmation_loops(
	IntervalIterator begin,
	IntervalIterator end,
	core::Size nres
)
{
	using protocols::loops::Loop;
	using protocols::loops::Loops;

	Loops loops;

	for ( IntervalIterator i = begin; i != end; ++i ) {
		core::Size const cut = (i->right - i->left +1)/2 + i->left; // just try cutting somewhere in the middle
		std::cout << "left: " << i->left << " right: " << i->right << " cut: " << cut << std::endl;

		protocols::loops::Loop loop;
		if ( i->left <= 3 ) { //N-term
			Loop loop( 1, i->right+2, cut);
			loops.add_loop(loop);
			//std::cout << "added loop " << "1" << ":" << i->right+2 << ":" << cut << std::endl;
		} else if ( i->right >= nres-2 ) { //cterm
			Loop loop(i->left-2, nres, cut);
			loops.add_loop(loop);
			//std::cout << "added loop " << i->left-2 << ":" << nres << ":" << cut << std::endl;
		} else {
			Loop loop(i->left-2, i->right+2, cut);
			loops.add_loop(loop);
			//std::cout << "added loop " << i->left-2 << ":" << i->right+2 << ":" << cut << std::endl;
		}
	}
	return loops;
}

/// @brief create fold tree from loops
/// @remarks This is a generic replacement function for the one in protocols::loops
///  and will be moved there in the near future.
core::kinematics::FoldTree
fold_tree_from_loops(
	core::pose::Pose const & pose,
	protocols::loops::Loops const & loops
);


/// @brief set a single loop fold tree
/// @remarks This is a generic replacement function for the one in protocols::loops
///  and will be moved there in the near future.
void
set_single_loop_fold_tree(
	core::pose::Pose & pose,
	protocols::loops::Loop const & loop
);

utility::vector1<bool>
parse_resfile_string_with_no_lockdown( core::pose::Pose const & pose, core::pack::task::PackerTask & the_task, std::string const & resfile_string
);

core::pack::task::TaskFactoryOP
remodel_generic_taskfactory();

void
fill_non_loop_cst_set(
	core::pose::Pose & pose,
	protocols::loops::Loops loops);

utility::vector1< core::Real > const
calc_rsd_sasa( core::pose::Pose const & pose );

void
apply_transformation(
	core::pose::Pose & mod_pose,
	std::list <core::Size> const & residue_list,
	numeric::xyzMatrix< core::Real > const & R, numeric::xyzVector< core::Real > const & preT, numeric::xyzVector< core::Real > const & postT
);

} // namespace methods
} // namespace forge
} // namespace protocols


#endif /* INCLUDED_protocols_forge_methods_util_HH */
