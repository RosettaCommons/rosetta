// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/BetaTurnDetection.hh
/// @brief  determine the presence and type of beta turn at a specific postion in a pose
/// @author Brian D. Weitzner (brian.weitzner@gmail.com)

#ifndef INCLUDED_protocols_features_BetaTurnDetection_HH
#define INCLUDED_protocols_features_BetaTurnDetection_HH

// Unit Headers
#include <utility/pointer/ReferenceCount.hh>
#include <protocols/features/BetaTurnDetection.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

// Utility Headers
#include <utility/vector1.fwd.hh>

// C++ Headers
#include <map>
#include <string>


namespace protocols{
namespace features{

enum RamachandranHash {
	A = 1,
	B,
	L,
	E,
	number_of_ramachandran_hashes = E
};

class BetaTurnDetection : public utility::pointer::ReferenceCount {
public:
	BetaTurnDetection();

	BetaTurnDetection( BetaTurnDetection const & );

	virtual ~BetaTurnDetection();

	/// @brief return string with class name
	virtual std::string
	type_name() const;

public:
	core::Size beta_turn_length() const { return beta_turn_length_; }
	core::Real beta_turn_distance_cutoff() const { return beta_turn_distance_cutoff_; }

public:
	static std::map< std::string, std::string > const & get_conformation_to_turn_type_map();
	static utility::vector1< std::string > const & get_valid_ramachandran_hashes();

	bool all_turn_residues_are_on_the_same_chain( core::pose::Pose const & pose, Size first_residue ) const;

	bool residue_range_is_protein( core::pose::Pose const & pose, Size range_begin, Size range_end ) const;

	bool beta_turn_present( core::pose::Pose const & pose, Size first_residue ) const;

	std::string const & beta_turn_type( core::pose::Pose const & pose, Size first_residue ) const;
	
	std::string determine_ramachandran_hash( core::pose::Pose const & pose, core::Size first_residue ) const;

	std::string determine_ramachandran_hash_for_residue_with_dihedrals( core::Real phi, core::Real psi, core::Real omega ) const;

	void validate_ramachandran_hash( std::string & rama_hash ) const;

private:
	core::Size const beta_turn_length_;
	core::Real const beta_turn_distance_cutoff_;

};

} // features namespace
} // protocols namespace

#endif //INCLUDED_protocols_features_BetaTurnDetection_HH
