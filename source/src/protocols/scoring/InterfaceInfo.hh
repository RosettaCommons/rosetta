// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/InterchainEnergy.cc
/// @brief  Statistically derived rotamer pair potentials
/// @details For docking (or between chains) only those residues at the interface
///						and between the two interfaces need to be evaluated
/// @author Monica Berrondo

#ifndef INCLUDED_protocols_scoring_InterfaceInfo_hh
#define INCLUDED_protocols_scoring_InterfaceInfo_hh

#include <core/types.hh>

// Unit headers
#include <protocols/scoring/InterfaceInfo.fwd.hh>

// Package headers
#include <protocols/scoring/Interface.hh>
#include <core/conformation/Residue.fwd.hh>

#include <basic/datacache/CacheableData.hh>

#include <ObjexxFCL/FArray1D.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace scoring {

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Keep track of the interface information

class InterfaceInfo : public basic::datacache::CacheableData {
public:

	/// @brief Default constructor with rb jump = 1
	InterfaceInfo(): calculated_(false)
	{
		rb_jump_.push_back( 1 );
		distance_ = 6.0; //default
		initialize();
	}

	/// @brief Constructor with arguments for non-default rb jump
	InterfaceInfo( core::Size rb_jump_in ): calculated_(false), distance_ (6.0)
	{
		rb_jump_.push_back( rb_jump_in );
		distance_ = 6.0; //default
		initialize();
	}

	/// @brief Constructor with arguments for multiple jumps
	InterfaceInfo( utility::vector1_size rb_jump_in ): calculated_(false)
	{
		rb_jump_ = rb_jump_in;
		distance_ = 6.0; //default
		initialize();
	}

	InterfaceInfo( InterfaceInfo const & src );

	basic::datacache::CacheableDataOP
	clone() const
	{
		return basic::datacache::CacheableDataOP( new InterfaceInfo( *this ) );
	}

	/// @brief Removes all jumps from the interface calculation
	void clear_jumps(){ rb_jump_.clear(); }

	/// @brief Adds another jump to the interface calculation, for example
	///for multi-body docking
	void add_jump( core::Size jump_in ){ rb_jump_.push_back( jump_in ); }

	/// @brief Sets the distance cutoff for interface calculations
	void distance( core::Real distance_in ){ distance_ = distance_in; }

	/// @brief Returns if interface calculation is up to date
	bool calculated() const { return calculated_; }

	/// @brief Returns if interface calculation is up to date
	bool & calculated() { return calculated_; }

	/// @brief Returns the number of jumps that are being used
	/// in interface calculations
	core::Size num_jump() const { return num_jump_; }

	/// @brief Sets up InterfaceInfo members such as interface_list_
	///based on variables from construction
	void
	initialize();

	/// @brief Returns whether a residue is at any of the interfaces
	bool
	is_interface( core::conformation::Residue rsd ) const ;

	/// @brief Returns whether the two residues are considered a
	///residue pair at any of the interfaces
	bool
	is_pair(
		core::conformation::Residue rsd1,
		core::conformation::Residue rsd2
	) const;

	/// @brief Returns the number of resides at the interface defined
	///by jump_num
	core::Size
	interface_nres( core::Size jump_num ) const;

	/// @brief Calculates the interface for all jumps specified in
	///rb_jump_
	void
	calculate( core::pose::Pose const & pose);

	core::Size
	closest_interface_residue( core::pose::Pose const & pose, core::Size src_rsd, core::Real & distance ) const;

	InterfaceCOP interface( core::Size interface_num ) const;

private:
	bool calculated_;
	core::Size num_jump_;
	core::Real distance_;

	utility::vector1_size rb_jump_;

	utility::vector1< InterfaceOP > interface_list_;
};

} // ns scoring
} // ns protocols

#endif
