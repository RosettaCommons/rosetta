// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/scoring/InterfaceInfo.cc
/// @brief  Statistically derived rotamer pair potentials
/// @detailed For docking (or between chains) only those residues at the interface
///						and between the two interfaces need to be evaluated
/// @author Monica Berrondo


// Unit headers
#include <protocols/scoring/InterchainPotential.hh>

// AUTO-REMOVED #include <core/scoring/AtomVDW.hh>
#include <core/scoring/EnvPairPotential.hh>
// AUTO-REMOVED #include <core/scoring/ScoringManager.hh>

// Package headers

// AUTO-REMOVED #include <core/scoring/EnergyGraph.hh>

// Project headers
#include <core/chemical/AA.hh>
// AUTO-REMOVED #include <core/chemical/VariantType.hh>
// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/kinematics/FoldTree.hh>
// AUTO-REMOVED #include <basic/database/open.hh>
// AUTO-REMOVED #include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/pose/datacache/CacheableDataType.hh>
// AUTO-REMOVED #include <basic/datacache/BasicDataCache.hh>

// Utility headers
// AUTO-REMOVED #include <utility/io/izstream.hh>
// AUTO-REMOVED #include <utility/utility.functions.hh>

// just for debugging
//#include <ObjexxFCL/format.hh>

#include <basic/Tracer.hh>

#include <protocols/scoring/InterfaceInfo.hh>
#include <utility/vector1.hh>

using basic::T;
using basic::Error;
using basic::Warning;


// C++

namespace protocols {
namespace scoring {

InterfaceInfo::InterfaceInfo( InterfaceInfo const & src ) :
	CacheableData(), calculated_(false)
{
	rb_jump_ = src.rb_jump_;
	num_jump_ = src.num_jump_;
	distance_ = src.distance_;

	interface_list_ = src.interface_list_;
}

void
//InterfaceInfo::initialize( pose::Pose const & pose )
InterfaceInfo::initialize()
{

	//get the number of jumps in the fold tree
	num_jump_ = rb_jump_.size();

	interface_list_.resize(num_jump_);

	//initialize interface objects for each interface in pose
	for (core::Size i = 1; i <= num_jump_; i++){
		interface_list_[i] = 	new protocols::scoring::Interface( rb_jump_[i] );
		interface_list_[i]->distance(6.0);
		}

}

bool
InterfaceInfo::is_interface(
	core::conformation::Residue rsd
	) const
{

	bool is_interface(false);
	for (core::Size i = 1; i <= num_jump_; i++){
		if (interface_list_[i]->is_interface( rsd )) is_interface = true;
		}

	return is_interface;
}

bool
InterfaceInfo::is_pair(
	core::conformation::Residue rsd1,
	core::conformation::Residue rsd2
	) const
{

	bool is_pair(false);

	for(core::Size i = 1; i<= num_jump_; i++){
		if (interface_list_[i]->is_pair(rsd1, rsd2)) is_pair = true;
	}

	return is_pair;
}

core::Size
InterfaceInfo::interface_nres( core::Size jump_num ) const
{
	return interface_list_[jump_num]->interface_nres();
}

void
InterfaceInfo::calculate( core::pose::Pose const & pose )
{
	for (core::Size i = 1; i <= num_jump_; i++){
		interface_list_[i]->calculate( pose );
		}
}

core::Size
InterfaceInfo::closest_interface_residue( core::pose::Pose const & pose, core::Size src_rsd, core::Real & distance ) const
{
	core::Size ret_rsd (0);
	core::Real min_distance (100000.0), temp_dist;
	for ( core::Size i = 1; i <= num_jump_; ++i ){
		core::Size temp_rsd = interface_list_[i]->closest_interface_residue( pose, src_rsd, temp_dist );
		if ( temp_dist < min_distance ) {
			ret_rsd = temp_rsd;
			min_distance = temp_dist;
		}
	}
	distance = temp_dist;
	return ret_rsd;
}

}
}
