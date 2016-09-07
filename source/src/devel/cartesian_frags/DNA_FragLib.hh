// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief

#ifndef INCLUDED_devel_cartesian_frags_DNA_FragLib_hh
#define INCLUDED_devel_cartesian_frags_DNA_FragLib_hh

#include <devel/cartesian_frags/DNA_FragLib.fwd.hh>
#include <devel/cartesian_frags/CartesianFragment.hh>


// libRosetta headers


#include <core/pose/Pose.hh>


#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>


// // // C++ headers
#include <string>
#include <map>

#include <utility/vector1.hh>


// //silly using/typedef

namespace devel {
namespace cartesian_frags {

/// @brief  A container class for cartesian fragments of DNA structure

class DNA_FragLib : public utility::pointer::ReferenceCount {
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	~DNA_FragLib() override;

	utility::vector1< CartesianFragment > const &
	base_pairs( std::string const & bp ) const
	{
		assert( base_pairs_.count(bp) );
		return base_pairs_.find( bp )->second;
	}

	utility::vector1< CartesianFragment > const &
	base_steps( std::string const & bp ) const
	{
		assert( base_steps_.count(bp) );
		return base_steps_.find( bp )->second;
	}

	void
	add_base_pair(
		std::string const & bp,
		CartesianFragment const & frag
	)
	{
		if ( !base_pairs_.count( bp ) ) base_pairs_[ bp ];
		base_pairs_.find( bp )->second.push_back( frag );
	}

	void
	add_base_step(
		std::string const & bs,
		CartesianFragment const & frag
	)
	{
		if ( !base_steps_.count( bs ) ) base_steps_[ bs ];
		base_steps_.find( bs )->second.push_back( frag );
	}

	core::pose::Pose &
	suite_pose() const;

public:
	utility::vector1< CartesianFragment > sugars;
	utility::vector1< CartesianFragment > forward_suites;
	utility::vector1< CartesianFragment > backward_suites;


private:
	std::map< std::string, utility::vector1< CartesianFragment > > base_pairs_;
	std::map< std::string, utility::vector1< CartesianFragment > > base_steps_;

	mutable core::pose::Pose suite_pose_;
};


typedef utility::pointer::shared_ptr< DNA_FragLib > DNA_FragLibOP;


void
build_frag_libraries(
	utility::vector1< std::string > const & files,
	DNA_FragLib & lib
);


void
setup_suite_pose( core::pose::Pose & suite_pose );


} // ns cartesian_frags
} // ns devel

#endif
