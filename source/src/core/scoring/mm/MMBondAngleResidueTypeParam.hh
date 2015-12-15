// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/mm/MMBondAngleResidueTypeParam.hh
/// @brief  Class to store bond angle parameters for a particular ResidueType
/// @author Colin A. Smith (colin.smith@ucsf.edu)


#ifndef INCLUDED_core_scoring_mm_MMBondAngleResidueTypeParam_hh
#define INCLUDED_core_scoring_mm_MMBondAngleResidueTypeParam_hh

// Unit headers
#include <core/scoring/mm/MMBondAngleResidueTypeParam.fwd.hh>

// Project headers
#include <core/chemical/ResidueType.fwd.hh>
#include <core/scoring/mm/MMBondAngleLibrary.fwd.hh>
#include <core/types.hh>

// Utility header
#include <utility/keys/Key2Tuple.hh>
#include <utility/keys/Key3Tuple.hh>

// C++ headers
#include <map>

#include <utility/vector1_bool.hh>


namespace core {
namespace scoring {
namespace mm {

typedef utility::keys::Key2Tuple< Size, Size > two_atom_set;
typedef utility::keys::Key3Tuple< Size, Size, Size > three_atom_set;

class MMBondAngleResidueTypeParam
{

public:

	/// @brief Default ctor
	MMBondAngleResidueTypeParam();

	/// @brief initialize the parameters
	void
	init(
		core::chemical::ResidueType const & residue_type,
		MMBondAngleLibrary const & mm_bondangle_library,
		bool use_residue_type_theta0,
		utility::vector1<std::string> const & central_atoms_to_score
	);

	/// @brief get number of intraresidue bond angles
	Size
	num_bondangles() const
	{
		return bondangle_atom_sets_.size();
	}

	/// @brief Return the indices for the set of atoms that define a particular
	/// intraresidue angle
	three_atom_set const &
	bondangle( Size const bondang ) const
	{
		return bondangle_atom_sets_[ bondang ];
	}

	/// @brief get Ktheta for a particular intraresidue angle
	core::Real
	Ktheta( Size const bondang ) const
	{
		return Ktheta_[ bondang ];
	}

	/// @brief get Ktheta for a particular intraresidue angle
	core::Real
	theta0( Size const bondang ) const
	{
		return theta0_[ bondang ];
	}

	/// @brief Returns the list of all of the indices of all the intraresidue
	/// bond angles a particular atom is involved in.
	/// Useful for calculating the derivatives for an atom.
	utility::vector1< Size > const &
	bondangles_for_atom( Size atomno ) const
	{
		return bondangles_for_atom_[ atomno ];
	}

	/// @brief Returns the index of the intraresidue bond angle, 0 if not found
	core::Size
	bondangle_index( three_atom_set const & atom_set ) const
	{
		std::map< three_atom_set, core::Size >::const_iterator iter(bondangle_index_.find(atom_set));

		if ( iter == bondangle_index_.end() ) return 0;

		return iter->second;
	}

	/// @brief number of ResidueConnections, counting polymeric residue connections
	Size
	n_residue_connections() const
	{
		return connection_atom_sets_.size();
	}

	/// @brief number of ResidueConnections, counting polymeric residue connections
	Size
	n_connection_pairs( Size const connection ) const
	{
		return connection_atom_sets_[ connection ].size();
	}

	/// @brief Return the indices for the set of two atoms that form part of a interresidue bond angle
	two_atom_set const &
	connection_pair( Size const connection, Size const bondang ) const
	{
		return connection_atom_sets_[ connection ][ bondang ];
	}

	/// @brief Return ResidueType derived theta0 for a interresidue bond angle
	core::Real
	connection_theta0( Size const connection, Size const bondang ) const
	{
		return connection_theta0_[ connection ][ bondang ];
	}

	/// @brief Return whether to use ResidueType derived theta0 for a interresidue bond angle
	bool
	connection_use_theta0( Size const connection, Size const bondang ) const
	{
		return connection_use_theta0_[ connection ][ bondang ];
	}

	/// @brief Returns the index of the interresidue bond angle, 0 if not found
	core::Size
	connection_index( Size const connection, two_atom_set const & atom_set ) const
	{
		std::map< two_atom_set, core::Size >::const_iterator iter(connection_index_[ connection ].find(atom_set));

		if ( iter == connection_index_[ connection ].end() ) return 0;

		return iter->second;
	}

	/// @brief stream << MMBondAngleResidueTypeParam
	friend
	std::ostream &
	operator <<(
		std::ostream & os,
		MMBondAngleResidueTypeParam const & residue_type_param
	);

private:

	/// @brief vector of sets of atoms that make up bond angles in the residue
	utility::vector1< three_atom_set > bondangle_atom_sets_;
	/// @brief vector of Ktheta values for each set
	utility::vector1< core::Real > Ktheta_;
	/// @brief vector of theta0 values for each set
	utility::vector1< core::Real > theta0_;
	/// @brief all intra-residue bond angles that each atom "participates" in.
	utility::vector1< utility::vector1< Size > > bondangles_for_atom_;
	/// @brief map to lookup a bondangle given its three atoms
	std::map< three_atom_set, core::Size > bondangle_index_;

	/// @brief vector of vector of pairs of atoms that take part in interresidue bond angles
	utility::vector1< utility::vector1< two_atom_set > > connection_atom_sets_;
	/// @brief vector of vector of theta0 values for pairs of atoms that take part in interresidue bond angles
	utility::vector1< utility::vector1< core::Real > > connection_theta0_;
	/// @brief vector of vector of booleans indicating whether to use theta0_ or lookup on the fly
	utility::vector1< utility::vector1< bool > > connection_use_theta0_;
	/// @brief vector of maps for lookup of interresidue theta0 parameters
	utility::vector1< std::map< two_atom_set, core::Size > > connection_index_;

	// The case in which a bond spanse two residue connections simultaneously should be able to be handled. One would
	// need to write code to detect multiple residue connections to the same atom, build both connected atoms
	// simultaneously, then calculate the theta0 value. Some additional instance variable data structures would need to
	// be incorporated to store/lookup the information.

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // namespace mm
} // namespace scoring
} // namespace core


#endif // INCLUDED_core_mm_MMBondAngleResidueTypeParam_HH
