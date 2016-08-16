// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/UnfoldedStatePotential.hh
/// @brief  Unfolded state energies based on energies of residues in fragments, declaration (header) file
/// @author Ron Jacak (ronj@email.unc.edu)

#ifndef INCLUDED_core_scoring_UnfoldedStatePotential_hh
#define INCLUDED_core_scoring_UnfoldedStatePotential_hh

// Unit Headers
#include <core/scoring/UnfoldedStatePotential.fwd.hh>

// Package headers

// Project headers

#include <core/scoring/EnergyMap.hh>
#include <core/pose/Pose.fwd.hh>

// c++ headers
#include <map>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {


/// @remarks
/// making this a separate class because it relies on a database file. rather than putting the code to read the database file
/// in the energy method class, it seems like the design we're following is to have a class like this one that reads the file
/// and provides the lookup into the data structure holding the database information and a separate class for the energy method
/// implementation.
///
class UnfoldedStatePotential : public utility::pointer::ReferenceCount {

public:
	/// @brief ctor - calls the function which reads in the database file
	UnfoldedStatePotential( std::string const & filename );
	virtual ~UnfoldedStatePotential();

	/// @brief returns the database values for an aa in the unfolded state - these are unweighted values!
	void
	raw_unfolded_state_energymap( std::string const & aa_name3, scoring::EnergyMap & e ) const;

	/// @brief returns the unweighted unfolded state energy for the whole pose as an emap (i.e. broken up by score type)
	void
	pose_raw_unfolded_state_energymap( pose::Pose const & pose, scoring::EnergyMap & e ) const;

	/// @brief returns an emap of the energy method weights specfied in the unfolded energy file
	scoring::EnergyMap
	get_unfoled_potential_file_weights() const;

private:
	/// @brief Read the amino acid energy file
	void read_database_file( std::string const & filename );


private:
	/// @brief Unfolded state energies by residue
	std::map< std::string, scoring::EnergyMap > unfolded_energy_;

	/// @brief energy method weights listed in the energies file
	scoring::EnergyMap unfolded_potential_file_weights_;

};

} // namespace scoring
} // namespace core


#endif // INCLUDED_core_scoring_UnfoldedStatePotential_HH
